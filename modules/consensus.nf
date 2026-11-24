#!/usr/bin/env nextflow
process Consensus {
    cpus 1
    tag "${name}" 

    input:
    tuple val(name), path(jaffal_csv), path(longgf_log), path(ctat_tsv), path(ctat_html)

    output:
    tuple val(name), path("jaffal_consensus_NOnonmalignant.csv"), path("jaffal_consensus_read_and_classification_filtered_NOnonmalignant.csv"), path("jaffal_breakpoint_consensus_NOnonmalignant.csv"), path("longgf_consensus_NOnonmalignant.csv"), path("longgf_consensus_read_and_classification_filtered_NOnonmalignant.csv"), path("longgf_breakpoint_consensus_NOnonmalignant.csv"), path("ctat_consensus_NOnonmalignant.csv"), path("ctat_consensus_read_and_classification_filtered_NOnonmalignant.csv"), path("ctat_breakpoint_consensus_NOnonmalignant.csv")

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    #defining functions for the script to normalise gene pairs alphabetically (for comparison of fusions) and to match fusions within specified distance between breakpoints called by 2 tool
    condense_ig <- function(g) {
        case_when(
            str_detect(g, "^IGH[VDJ][0-9I]|^IGH(M|D|E|G[1-4P]|A[12]|EP[12])\$") ~ "IGH",
            str_detect(g, "^IGK[VJ][0-9]|^IGK(C|DEL)\$")                        ~ "IGK",
            str_detect(g, "^IGL[VJ][0-9]|^IGLC[0-9]\$|^IGLL[0-9]\$")             ~ "IGL",
            TRUE ~ g                                                            # non-IG genes unchanged
        )
    }
    #this function sorts the genes alphabetically (normalizes pair for comparison)
    condense_ig_pair <- function(x, sep = ":+") {
        map_chr(str_split(x, sep), function(g) {
            if (all(is.na(g))) return(NA_character_)
            paste(sort(condense_ig(g)), collapse = ":")
        })
    }

    #this function will create a new column fusion_norm which will include alphabetically sorted fusion gene names with IG genes condensed to IGH/K/L. It also filters fusions between 2 the same genes.
    filter_samegene <- function(fusion_df, fusion_col) {
        fusion_df |>
            mutate(fusion_norm = condense_ig_pair(pick({{ fusion_col }})[[1]])) |>
            separate( fusion_norm,
            into   = c("gene1", "gene2"),
            sep    = ":",
            remove = FALSE,
            extra  = "merge",
            fill   = "right" ) |>
            mutate(same_gene = gene1 == gene2) |>
            filter(!same_gene | is.na(same_gene)) |>   # keep rows with a missing partner; drop identical pairs
            select(-same_gene)
    }

    #reading the file with recurrent fusions from previous literature, filtering same gene fusions and normalising gene pairs for comparison
    nonmalignant_fusions <- read_csv("${params.nonmalignant_fusions_path}") |>
        mutate(fusion = paste0(up_gene, ':', dw_gene)) |>
        filter_samegene(fusion) 

    #reading JAFFAL output, filtering same gene fusions, and normalising gene pairs for comparison
    jaffal <- read_csv("${jaffal_csv}") |>
        mutate(`fusion genes` = str_replace_all(`fusion genes`, pattern = '::', replacement = ':')) |>
        filter(!(classification %in% c('PotentialTransSplicing', 'PotentialReadThrough'))) |>
        filter_samegene(`fusion genes`) |>
        mutate_at(vars(base1, base2), as.numeric) 
    
    # read every GF / SumGF line and pair them with a shared fusion_id - this will be used later in breakpoint-consensus filtering 
    longgf_in <- tibble(X1 = read_lines("${longgf_log}")) |>
      filter(str_detect(X1, "^(GF|SumGF)\t")) |>
      mutate(type      = str_extract(X1, "^[^\t]+"),
             genes = str_match(X1, "^[^\t]+\t([^ ]+)")[, 2],
             fusion_id = cumsum(type == "GF")) 

# reading CTAT-LR output, filtering same gene fusions, and normalising gene pairs for comparison
ctat <- read_tsv("${ctat_tsv}") |>
  mutate(`fusion genes` = str_replace_all(`#FusionName`, pattern = '--', replacement = ':')) |>
  filter_samegene(`fusion genes`) |>
  filter(!(fusion_norm %in% nonmalignant_fusions\$fusion_norm)) |>
  separate(LeftBreakpoint, into = c("chrom1", "base1", "strand1"), sep = ":") |>
  separate(RightBreakpoint, into = c("chrom2", "base2", "strand2"), sep = ":")

    # sanity check: each id has exactly one GF and one SumGF, with the same gene pair
    pair_check <- longgf_in |>
      group_by(fusion_id) |>
      summarise(n_gf  = sum(type == "GF"),
                n_sum = sum(type == "SumGF"),
                same  = n_distinct(genes) == 1)
    stopifnot(all(pair_check\$n_gf == 1), all(pair_check\$n_sum == 1), all(pair_check\$same))
    
    longgf_in <- longgf_in |> select(!c(genes, type))

    #formatting LongGF tibble using only rows starting with "SumGF", for consensus calling
    longgf <- longgf_in |>
      filter(grepl("SumGF", X1)) |>
      separate(col = X1, into = c("sumgf", "info"), sep = "\t") |>
      separate(col = info, into = c("fusion genes", "read_count", "chr1", "chr2"), sep = " ") |>
      separate(col = chr1, into = c("chrom1", "base1"), sep = ":") |>
      separate(col = chr2, into = c("chrom2", "base2"), sep = ":")

    #In LongGF SumGF tibble filtering same gene fusions, and normalising gene pairs for comparison
    longgf <- longgf |>
        filter_samegene(`fusion genes`) |>
        filter(!(fusion_norm %in% nonmalignant_fusions\$fusion_norm)) |>
      mutate(across(c(read_count, base1, base2), as.numeric))
    
    longgf_GF <- longgf_in |>
      filter(grepl("^GF", X1)) |>
      separate(col = X1, into = c("gf", "info"), sep = "\t") |>
      separate_wider_delim(cols = info, delim = " ", names_sep = "") |>
      rename("fusion genes" = info1) |>
      rename("read_count" = info2) |>
      filter_samegene(`fusion genes`) |>
  filter(!(fusion_norm %in% nonmalignant_fusions\$fusion_norm))
    
    ##############################################################
    #################  gene name-based filtering #################
    ##############################################################

    #generating list of gene partners in consensus between the 3 tools
    consensus_jaffal_longgf <- intersect(jaffal\$fusion_norm, longgf\$fusion_norm)
    consensus_jaffal_ctat <- intersect(jaffal\$fusion_norm, ctat\$fusion_norm)
    consensus_longgf_ctat <- intersect(longgf\$fusion_norm, ctat\$fusion_norm)
    consensus <- unique(c(consensus_jaffal_longgf, consensus_jaffal_ctat, consensus_longgf_ctat))

    #filtering out fusions recurrent in normal tissues genes and keeping only fusions are in consensus (based on gene fusion pairs names)
    jaffal_consensus <- jaffal |>
        filter(fusion_norm %in% consensus)  #selecting only fusion in consensus based on gene names

    #formatting LongGF tibble using rows starting with "GF" to use further as input to AGFusion annotator, and filtering based on gene name consensus
    longgf_consensus <- longgf_GF |>
      filter(fusion_norm %in% consensus) |> #filtering only fusion in consensus based on gene names
      select(!c(fusion_norm, gene1, gene2, fusion_id)) |> 
      unite(col = "X1", gf:`fusion genes`, sep = "\t") |>
      unite(col = "X1", everything(), sep = " ") 
    
    #formatting LongGF tibble using rows starting with "GF" to use further as input to AGFusion annotator, and filtering based on gene name consensus
    ctat_consensus <- ctat |>
        filter(fusion_norm %in% consensus)

    ##############################################################
    ##### read, HighConfidence and gene name-based filtering #####
    ##############################################################
    
    #JAFFAL: HighConfidence, enough spanning reads, not a known non-malignant fusion
    jaffal_filtered <- jaffal |>
        filter(classification == 'HighConfidence') |>
        filter(`spanning reads` >= ${params.min_read_count})

    #LongGF: enough supporting reads (needs read_count converted to numeric in `longgf`)
    longgf_read_filtered <- longgf |> filter(read_count >= ${params.min_read_count})
    #CTAT-LR: enough supporting reads 
    ctat_read_filtered <- ctat |> filter(num_LR >= ${params.min_read_count})

    #generating list of gene partners in consensus between the 3 tools, after quality-based filtering
    consensus_filtered_jaffal_longgf <- intersect(jaffal_filtered\$fusion_norm, longgf_read_filtered\$fusion_norm)
    consensus_filtered_jaffal_ctat <- intersect(jaffal_filtered\$fusion_norm, ctat_read_filtered\$fusion_norm)
    consensus_filtered_longgf_ctat <- intersect(longgf_read_filtered\$fusion_norm, ctat_read_filtered\$fusion_norm)
    consensus_filtered <- unique(c(consensus_filtered_jaffal_longgf, consensus_filtered_jaffal_ctat, consensus_filtered_longgf_ctat))

    jaffal_consensus_filtered <- jaffal_filtered |>
      filter(fusion_norm %in% consensus_filtered)

    ctat_consensus_filtered <- ctat_read_filtered |>
        filter(fusion_norm %in% consensus_filtered)

    ids_consensus_filtered <- longgf_read_filtered |>
    filter(fusion_norm %in% consensus_filtered) |>
      pull(fusion_id)
    
    longgf_consensus_filtered <- longgf_GF |>
      filter(fusion_id %in% ids_consensus_filtered) |>
      select(!c(fusion_norm, gene1, gene2, fusion_id)) |>
      unite(col = "X1", gf:`fusion genes`, sep = "\t") |>
      unite(col = "X1", everything(), sep = " ")

    ##############################################################
    ################# breakpoint-based filtering #################
    ##############################################################

    # put each tool's breakpoints in a common form: the two breakpoints are ordered
    # by (chromosome, position) so the same pair is comparable across tools
    canon_bp <- function(df) {
        df |>
            mutate(row_id = row_number(),
                across(c(base1, base2), as.numeric),
                swap  = chrom1 > chrom2 | (chrom1 == chrom2 & base1 > base2),
                chrA  = if_else(swap, chrom2, chrom1), posA = if_else(swap, base2, base1),
                chrB  = if_else(swap, chrom1, chrom2), posB = if_else(swap, base1, base2)) |>
            select(row_id, fusion_norm, chrA, posA, chrB, posB)
        }

    # all call pairs between two tools with the same gene pair, same chromosomes,
    # and both breakpoints within the tolerance
    match_pair <- function(a, b, tol) {
        inner_join(a, b, by = c("fusion_norm", "chrA", "chrB"),
                    suffix = c("_a", "_b"), relationship = "many-to-many") |>
            filter(abs(posA_a - posA_b) <= tol, abs(posB_a - posB_b) <= tol)
        }

    bp <- list(jaffal = canon_bp(jaffal_filtered),
                longgf = canon_bp(longgf_read_filtered),
                ctat   = canon_bp(ctat_read_filtered))

    # for every call, count how many OTHER tools have a matching breakpoint
    support <- map_dfr(combn(names(bp), 2, simplify = FALSE), function(p) {
        m <- match_pair(bp[[p[1]]], bp[[p[2]]], params.breakpoint_tolerance)
        bind_rows(tibble(tool = p[1], row_id = m\$row_id_a, partner = p[2]),
                    tibble(tool = p[2], row_id = m\$row_id_b, partner = p[1]))
    }) |>
        distinct() |>
        count(tool, row_id, name = "n_other_tools")   # 1 = supported by 2 of 3 tools, 2 = all 3

        # min_support = 1 gives 2-of-3 consensus; min_support = 2 gives all-3 consensus
    ids_for <- function(t, min_support = 1) {
        support |>
            filter(tool == t, n_other_tools >= min_support) |>
            pull(row_id) |>
            sort()
    }

    jaffal_breakpoint_consensus <- jaffal_filtered[ids_for("jaffal"), ]
    ctat_breakpoint_consensus   <- ctat_read_filtered[ids_for("ctat"), ]
    ids_breakpoint              <- longgf_read_filtered\$fusion_id[ids_for("longgf")]

    #formatting LongGF tibble using rows starting with "GF" to use further as input to AGFusion annotator, and filtering based on breakpoint consensus
    longgf_breakpoint_consensus <- longgf_GF |>
        filter(fusion_id %in% ids_breakpoint) |>
        select(!c(fusion_norm, gene1, gene2, fusion_id)) |>
        unite(col = "X1", gf:`fusion genes`, sep = "\t") |>
        unite(col = "X1", everything(), sep = " ", na.rm = TRUE)

    ##############################################################
    ####################### saving outputs #######################
    ##############################################################
        
    jaffal_breakpoint_consensus <- jaffal_breakpoint_consensus |> select(!c(fusion_norm, gene1, gene2))
    
    jaffal_consensus <- jaffal_consensus |> select(!c(fusion_norm, gene1, gene2))
    
    jaffal_consensus_filtered <- jaffal_consensus_filtered |> select(!c(fusion_norm, gene1, gene2))
    
    ctat_consensus <- ctat_consensus |> 
        select(!c(fusion_norm, gene1, gene2, `fusion genes`)) |> 
        unite(col = "LeftBreakpoint", chrom1:strand1, sep = "\t") |>
        unite(col = "RightBreakpoint", chrom2:strand2, sep = "\t")

    ctat_consensus_filtered <- ctat_consensus_filtered |> 
        select(!c(fusion_norm, gene1, gene2, `fusion genes`)) |> 
        unite(col = "LeftBreakpoint", chrom1:strand1, sep = "\t") |>
        unite(col = "RightBreakpoint", chrom2:strand2, sep = "\t")

    ctat_breakpoint_consensus <- ctat_breakpoint_consensus |> 
        select(!c(fusion_norm, gene1, gene2, `fusion genes`)) |> 
        unite(col = "LeftBreakpoint", chrom1:strand1, sep = "\t") |>
        unite(col = "RightBreakpoint", chrom2:strand2, sep = "\t")

    write_csv(jaffal_consensus, 'jaffal_consensus_NOnonmalignant.csv')

    write_csv(jaffal_consensus_filtered, 'jaffal_consensus_read_and_classification_filtered_NOnonmalignant.csv')

    write_csv(jaffal_breakpoint_consensus, 'jaffal_breakpoint_consensus_NOnonmalignant.csv')

    write_csv(longgf_consensus, 'longgf_consensus_NOnonmalignant.csv')

    write_csv(longgf_consensus_filtered, 'longgf_consensus_read_and_classification_filtered_NOnonmalignant.csv')
    
    write_csv(longgf_breakpoint_consensus, 'longgf_breakpoint_consensus_NOnonmalignant.csv')

    write_csv(ctat_consensus, 'ctat_consensus_NOnonmalignant.csv')

    write_csv(ctat_consensus_filtered, 'ctat_consensus_read_and_classification_filtered_NOnonmalignant.csv')

    write_csv(ctat_breakpoint_consensus, 'ctat_breakpoint_consensus_NOnonmalignant.csv')

    """
    stub:
    """
    touch jaffal_consensus_NOnonmalignant.csv
    touch jaffal_breakpoint_consensus_NOnonmalignant.csv
    touch jaffal_consensus_read_and_classification_filtered_NOnonmalignant.csv
    touch longgf_consensus_NOnonmalignant.csv
    touch longgf_consensus_read_and_classification_filtered_NOnonmalignant.csv
    touch longgf_breakpoint_consensus_NOnonmalignant.csv
    touch ctat_consensus_NOnonmalignant.csv
    touch ctat_consensus_read_and_classification_filtered_NOnonmalignant.csv
    touch ctat_breakpoint_consensus_NOnonmalignant.csv
    """
}
#!/usr/bin/env nextflow
process Consensus {
    cpus 1
    tag "${name}" 

    input:
    tuple val(name), path(jaffal_csv), path(longgf_log)

    output:
    tuple val(name), path("jaffal_consensus_NOnonmalignant.csv"), path("jaffal_breakpoint_consensus_NOnonmalignant.csv"), path("longgf_breakpoint_consensus_NOnonmalignant.csv")

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

    find_breakpoint_matches <- function(df1, df2, tolerance = 10, consensus_by = "fusion_norm") {
        df1 <- df1 |>
            rowwise() |>
            mutate(low_break = min(base1, base2), 
                high_break = max(base1, base2)) |>
            ungroup() |>
            rename_with(~ paste0(.x, "_1"))
        df2 <- df2 |>
            rowwise() |>
            mutate(low_break = min(base1, base2), 
                high_break = max(base1, base2)) |>
            ungroup() |>
            rename_with(~ paste0(.x, "_2"))
        by_cols <- setNames(
            paste0(consensus_by, "_2"),
            paste0(consensus_by, "_1"))
        df1 |>
            inner_join(
            df2,
            by = by_cols,
            relationship = "many-to-many") |> #produces all combinations between all fusion calls with different breakpoints in both tools
            filter(
            abs(low_break_1  - low_break_2)  <= tolerance,
            abs(high_break_1 - high_break_2) <= tolerance) |> #The tolerance filter keeps any combination where both breakpoints are close
            select(ends_with("_1")) |>
            rename_with(~ sub("_1\$", "", .x)) |>
            distinct() |>
            select(!c(high_break, low_break))
        }

    #reading the file with recurrent fusions from previous literature, filtering same gene fusions and normalising gene pairs for comparison
    nonmalignant_fusions <- read_csv("${params.nonmalignant_fusions_path}") |>
        mutate(fusion = paste0(up_gene, ':', dw_gene)) |>
        filter_samegene(fusion) 

    #reading jaffal output, filtering same gene fusions, and normalising gene pairs for comparison
    jaffal <- read_csv("${jaffal_csv}") |>
        mutate(`fusion genes` = str_replace_all(`fusion genes`, pattern = '::', replacement = ':')) |>
        filter(!(classification %in% c('PotentialTransSplicing', 'PotentialReadThrough'))) |>
        filter_samegene(`fusion genes`) |>
        mutate_at(vars(base1, base2), as.numeric) 

    #reading LongGF output
    longgf_in <- read_tsv("${longgf_log}", col_names = FALSE) 

    #formating LongGF tibble using only rows starting with "SumGF", for consensus calling
    longgf <- longgf_in |>
        filter(grepl("SumGF", X1)) |>
        separate(
            col = X1,
            into = c("sumgf", "info"),
            sep = "\t") |>
        separate(
            col = info,
            into = c("fusion genes", "read_count", "chr1", "chr2"),
            sep = " ") |>
        separate(
            col = chr1,
            into = c("chrom1", "base1"),
            sep = ":") |>
        separate(
            col = chr2,
            into = c("chrom2", "base2"),
            sep = ":")

    #In LongGF SumGF tibble filtering same gene fusions, and normalising gene pairs for comparison
    longgf <- longgf |>
        filter_samegene(`fusion genes`) |>
        mutate_at(vars(base1, base2), as.numeric) 

    #generating list of gene partners in consensus between the 2 tools
    consensus <- intersect(jaffal\$fusion_norm, longgf\$fusion_norm)

    #filtering out fusions recurrent in normal tissues genes and keeping only fusions are in consensus (based on gene fusion pairs names)
    jaffal_consensus <- jaffal |>
        filter(!(fusion_norm %in% nonmalignant_fusions\$fusion_norm)) |>
        filter(fusion_norm %in% consensus) |>
        select(!c(fusion_norm, gene1, gene2))
    
    #filtering genes that are in consensus based on gene fusion pairs names as well as breakpoint distance between tools, and filteting based on spanning reads
    jaffal_breakpoint_consensus <- jaffal |>
        filter(!(fusion_norm %in% nonmalignant_fusions\$fusion_norm)) |>
        filter(fusion_norm %in% consensus) |>
        find_breakpoint_matches(longgf, tolerance = ${params.breakpoint_tolerance}, consensus_by = "fusion_norm") |>
        filter(classification == 'HighConfidence') |>
        filter(`spanning reads` >= ${params.min_read_count})
    
    #formatting LongGF tibble using rows starting with "GF" to use further as input to AGFusion annotator, and filtering based on breakpoint consensus
    longgf_breakpoint_consensus <- longgf_in |>
        filter(grepl("^GF", X1)) |>
        separate(
            col = X1,
            into = c("gf", "info"),
            sep = "\t") |>
        separate_wider_delim(
            cols = info,
            delim = " ",
            names_sep = "") |>
        rename("fusion genes" = info1) |>
        rename("read_count" = info2) |>
        filter_samegene(`fusion genes`) |>
        filter(fusion_norm %in% jaffal_breakpoint_consensus\$fusion_norm) |>
        select(!c(fusion_norm, gene1, gene2)) |> 
        filter(read_count >= ${params.min_read_count}) |>
        unite(
            col = "X1",
            gf:`fusion genes`,
            sep = "\t") |>
        unite(
            col = "X1",
            everything(),
            sep = " ") 
    
    jaffal_breakpoint_consensus <- jaffal_breakpoint_consensus |> select(!c(fusion_norm, gene1, gene2))

    write_csv(jaffal_consensus, 'jaffal_consensus_NOnonmalignant.csv')

    write_csv(jaffal_breakpoint_consensus, 'jaffal_breakpoint_consensus_NOnonmalignant.csv')

    write_csv(longgf_breakpoint_consensus, 'longgf_breakpoint_consensus_NOnonmalignant.csv')

    """
    stub:
    """
    touch jaffal_consensus_NOnonmalignant.csv
    touch jaffal_breakpoint_consensus_NOnonmalignant.csv
    touch longgf_breakpoint_consensus_NOnonmalignant.csv
    """
}
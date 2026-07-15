#!/usr/bin/env nextflow
process Consensus {
    cpus 1
    tag "${name}" 

    input:
    tuple val(name), path(jaffal_csv), path(longgf_log)

    output:
    tuple val(name), path("consensus_NOnonmalignant.csv"), path("consensus_NOnonmalignant_breakpoints.csv")

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    #defining functions for the script to normalise gene pairs alphabetically (for comparison of fusions) and to match fusions within specified distance between breakpoints called by 2 tools
    normalize_pair <- function(x) {
        genes <- strsplit(x, ":")
        sapply(genes, function(g) paste(sort(g), collapse = ":"))
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
            paste0(consensus_by, "_1")
        )
        
        df1 |>
            inner_join(
            df2,
            by = by_cols,
            relationship = "many-to-many"
            ) |>
            filter(
            abs(low_break_1  - low_break_2)  <= tolerance,
            abs(high_break_1 - high_break_2) <= tolerance
            ) |>
            select(ends_with("_1")) |>
            rename_with(~ sub("_1\$", "", .x)) |>
            distinct() |>
            select(!c(high_break, low_break))
    }
        
    #reading the file with recurrent fusions from previous literature 
    nonmalignant_fusions <- read_csv("${params.nonmalignant_fusions_path}") |>
        mutate(fusion = paste0(up_gene, ':', dw_gene)) |>
        mutate(fusion_norm = normalize_pair(fusion))

    #reading and adjusting the format of jaffal and longgf results, filtering based on read count, and normalising gene pairs for comparison 
    jaffal <- read_csv("${jaffal_csv}") |>
        filter(!(classification == 'PotentialTransSplicing')) |>
        mutate(fusion_norm = normalize_pair(`fusion genes`)) |>
        mutate_at(vars(base1, base2), as.numeric) |>
        filter(`spanning reads` > ${params.min_read_count})

    longgf_in <- read_tsv("${longgf_log}", col_names = FALSE) 

    longgf <- longgf_in |>
        filter(grepl("SumGF", X1)) |>
        separate(
            col = X1,
            into = c("sumgf", "info"),
            sep = "\t"
        ) |>
        separate(
            col = info,
            into = c("fusion genes", "read_count", "chr1", "chr2"),
            sep = " "
        ) |>
        separate(
            col = chr1,
            into = c("chrom1", "base1"),
            sep = ":"
        )  |>
        separate(
            col = chr2,
            into = c("chrom2", "base2"),
            sep = ":"
        )

    longgf <- longgf |>
        mutate(fusion_norm = normalize_pair(`fusion genes`)) |>
        mutate_at(vars(base1, base2), as.numeric) |>
        filter(read_count > ${params.min_read_count})

    #generating list of gene partenrs in consensus between the 2 tools
    consensus <- intersect(jaffal\$fusion_norm, longgf\$fusion_norm)

    #filtering genes that are in consensus (based on gene fusion pairs) and that are not recurrent in normal tissues
    jaffal_consensus <- jaffal |>
        filter(!(fusion_norm %in% nonmalignant_fusions\$fusion_norm)) |>
        filter(fusion_norm %in% consensus) |>
        select(!c(fusion_norm))
    
    #filtering genes that are in consensus based on gene fusion pairs as well as breakpoint distance within 
    jaffal_breakpoint_consensus <- jaffal |>
        filter(!(fusion_norm %in% nonmalignant_fusions\$fusion_norm)) |>
        filter(fusion_norm %in% consensus) |>
        find_breakpoint_matches(longgf, tolerance = ${params.breakpoint_tolerance}, consensus_by = "fusion_norm") |>
        select(!c(fusion_norm))
    
    write_csv(jaffal_consensus, 'consensus_NOnonmalignant.csv')

    write_csv(jaffal_breakpoint_consensus, 'consensus_NOnonmalignant_breakpoints.csv')

    """
    stub:
    """
    touch consensus_NOnonmalignant.csv
    touch consensus_NOnonmalignant_breakpoints.csv
    """
}
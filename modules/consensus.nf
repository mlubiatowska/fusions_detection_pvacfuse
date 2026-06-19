#!/usr/bin/env nextflow
process Consensus {
    cpus params.cpus
    tag "${name}" 

    input:
    tuple val(name), path(jaffal_csv), path(longgf_log)

    output:
    tuple val(name), path("jaffal_longgfconsensus.csv")

    script:
    """
    #!/usr/bin/R

    library(tidyverse)

    jaffal <- read_csv("${jaffal_csv}") 

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

    normalize_pair <- function(x) {
    genes <- strsplit(x, ":")
    sapply(genes, function(g) paste(sort(g), collapse = ":"))
    }

    jaffal <- jaffal |>
    mutate(fusion_norm = normalize_pair(\`fusion genes\`))

    longgf <- longgf |>
    mutate(fusion_norm = normalize_pair(\`fusion genes\`))

    consensus <- intersect(jaffal\$fusion_norm, longgf\$fusion_norm)

    jaffal_out <- jaffal |>
    filter(fusion_norm %in% consensus) |>
    select(!fusion_norm)

    write_csv(jaffal_out, "${jaffal_longgfconsensus.csv}")
    
    """
    stub:
    """
    touch jaffal_longgfconsensus.csv
    """
}
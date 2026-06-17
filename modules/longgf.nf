#!/usr/bin/env nextflow

process Longgf {
    cpus params.cpus
    tag "${name}" 

    module params.samtools

    input:
    tuple val(name), path(unfiltered_rna_bam)

    output:
    tuple val(name), path("LongGF.${name}.log")

    script:
    """
    samtools sort -n -@ 2 \
        -o ${name}_Baseline_reads_aln_name_sorted.bam \
        ${unfiltered_rna_bam}

    LongGF \
    ${name}_Baseline_reads_aln_name_sorted.bam \
    ${params.gtf} \
    100 50 200 0 0 2 \
    > LongGF.${name}.log
    
    """
    stub:
    """
    touch LongGF.${name}.log
    """
}



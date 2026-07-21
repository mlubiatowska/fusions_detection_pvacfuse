#!/usr/bin/env nextflow
process LonggfPrep {
    cpus 1
    tag "${name}" 

    module params.samtools

    input:
    tuple val(name), path(unfiltered_rna_bam)

    output:
    tuple val(name), path("${name}_Baseline_reads_aln_name_sorted.bam")

    script:
    """
    samtools sort -n -@ 2 \
        -o ${name}_Baseline_reads_aln_name_sorted.bam \
        ${unfiltered_rna_bam}
    
    """
    stub:
    """
    touch ${name}_Baseline_reads_aln_name_sorted.bam
    """
}


process Longgf {
    tag "${name}" 

    input:
    tuple val(name), path(sorted_rna_bam)

    output:
    tuple val(name), path("LongGF.${name}.log")

    script:
    """

    LongGF \
    ${sorted_rna_bam} \
    ${params.gtf} \
    100 50 200 0 0 2 \
    > LongGF.${name}.log
    
    """
    stub:
    """
    touch LongGF.${name}.log
    """
}



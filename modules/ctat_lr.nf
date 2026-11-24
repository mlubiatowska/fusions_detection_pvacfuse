#!/usr/bin/env nextflow

process CTAT_LR {
    cpus 4
    tag "${name}" 

    input:
    tuple val(name), path(rna_fastq)

    output:
    tuple val(name), path("ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_predictions.tsv"), path("ctat_LR_fusion_outdir/ctat-LR-fusion.fusion_inspector_web.html")

    script:
    """
    ctat-LR-fusion \
        -T ${rna_fastq} \
        --genome_lib_dir  ${params.genome_lib} \
        --CPU 4 \
        --vis
    """
    stub: 
    """
    touch 
    """
}
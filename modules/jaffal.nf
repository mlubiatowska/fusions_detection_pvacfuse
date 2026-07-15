#!/usr/bin/env nextflow

process Jaffal {
    cpus 4
    tag "${name}" 

    module params.python

    input:
    tuple val(name), path(rna_fastq)

    output:
    tuple val(name), path("jaffa_results.csv")

    script:
    """
    bpipe run /JAFFA/JAFFAL.groovy ${rna_fastq}
    """
    stub: 
    """
    touch jaffa_results.csv
    """
}

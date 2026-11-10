#!/usr/bin/env nextflow

process Jaffal {
    tag "${name}" 

    module params.python

    input:
    tuple val(name), path(rna_fastq)

    output:
    tuple val(name), path("jaffa_results.csv")

    script:
    """
    bpipe run -p genome=hs1 -p annotation=catLiftOffGenesV1 /JAFFA/JAFFAL.groovy ${rna_fastq}
    """
    stub: 
    """
    touch jaffa_results.csv
    """
}

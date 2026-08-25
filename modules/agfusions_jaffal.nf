#!/usr/bin/env nextflow

process AGFusion_Jaffal {
    cpus 1
    tag "${name}"  

    input:
    tuple val (name), path(consensus_breakpoints)

    output:
    tuple val(name), path("agfusion_jaffal")

    script:
    """
    agfusion batch \
      -f ${consensus_breakpoints} \
      -a jaffa \
      -db ${params.db}  \
      -o agfusion_jaffal \
      --middlestar \
      --noncanonical

    """
    stub:
    """
    mkdir -p agfusion_jaffal
    """
}
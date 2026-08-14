#!/usr/bin/env nextflow

process AGFusion_LongGF {
    cpus 1
    tag "${name}"  

    input:
    tuple val (name), path(), path(longgf_consensus_breakpoints)

    output:
    tuple val(name), path("agfusion_longgf")

    script:
    """
    agfusion batch \
      -f ${longgf_consensus_breakpoints} \
      -a jaffa \
      -db ${params.db}  \
      -o agfusion_longgf \
      --middlestar \
      --noncanonical

    """
    stub:
    """
    mkdir -p agfusion
    """
}
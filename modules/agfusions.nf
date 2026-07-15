#!/usr/bin/env nextflow

process AGFusion {
    cpus 1
    tag "${name}"  

    input:
    tuple val (name), path(consensus_jaffal), path(consensus_breakpoints) //LongGF.${name}.log

    output:
    tuple val(name), path("agfusion")

    script:
    """
    agfusion batch \
      -f consensus_breakpoints \
      -a longgf \
      -db ${params.db}  \
      -o agfusion \
      --middlestar \
      --noncanonical

    """
    stub:
    """
    mkdir -p agfusion
    """
}
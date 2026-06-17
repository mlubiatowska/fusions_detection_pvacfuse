#!/usr/bin/env nextflow

process NanoFG_Somatic {
    cpus params.cpus
    tag "${name}" 

    module load SAMtools/1.11
    module load  Python/3.11.5-GCCcore-13.2.0

    input:
    tuple val(name), path(nanofg_bam)

    output:
    tuple val(name), path("")

    script:
    """"

    source ${NANOFG_VENV}/activate

    bash ${NANOFG}/NanoFG.sh \
        -b ${OUT}/${SAMPLE}_somatic_baseline_nochr.bam \
        -n ${SAMPLE} \
        -venv ${NANOFG_VENV}/activate \
        -o ${OUT}/BAM_out \
        -gtf ${GTF} \
        -t 6 \
        --without_last

    """
    stub:
    """
    
    """
}



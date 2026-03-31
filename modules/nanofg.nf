#!/usr/bin/env nextflow

process NanoFG {
    cpus params.cpus
    tag "${name}" 

    module load SAMtools/1.11
    module load  Python/3.11.5-GCCcore-13.2.0

    input:
    tuple val(name), path(dna_bam)

    output:
    tuple val(name), path("")

    script:
    """
    mkdir ${name}_fastq

    samtools fastq \
        --threads ${params.cpus} \
        -n -0 ./${name}_fastq/${name}_baseline.fastq \
        ${dna_bam}

    source ${NANOFG_VENV}/activate

    bash ${NANOFG}/NanoFG.sh \
        -f ./${name}_fastq/ \
        -n ${name} \
        -o . \
        -t ${params.cpus}
    """
    stub:
    """
    
    """
}



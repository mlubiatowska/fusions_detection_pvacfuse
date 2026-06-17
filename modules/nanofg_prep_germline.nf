#!/usr/bin/env nextflow

process NanoFG_BAM {
    cpus params.cpus
    tag "${name}" 

    module load SAMtools/1.11

    input:
    tuple val(name), path(dna_bam)

    output:
    tuple val(name), path("")

    script:
    """
    samtools view -H ${DNA_BAM} | sed -E 's/SN:chr([0-9XYM]+)/SN:\1/' > ${OUT}/header.sam

    samtools reheader ${OUT}/header.sam ${DNA_BAM} > ${OUT}/${SAMPLE}_somatic_baseline_nochr.bam

    samtools index ${OUT}/${SAMPLE}_somatic_baseline_nochr.bam
        
    """
    stub:
    """
    
    """
}



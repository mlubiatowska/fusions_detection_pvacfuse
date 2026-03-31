#!/usr/bin/env nextflow
nextflow.enable.dsl     = 2
/*
* Pipeline parameters
*/


include { Jaffal } from './modules/jaffal.nf'
include { Longgf } from './modules/longgf.nf'
include { NanoFG } from './modules/nanofg.nf'
include { Consensus } from './modules/consensus.nf'
include { AGFusion } from './modules/agfusions.nf'
include { PvacFuse } from './modules/pvacfuse.nf'


workflow {
    main:
    // Parse input CSV and create two separate channels   
    input_data = Channel
        .fromPath(params.input)
        .splitCsv(header: ['name', 'unfiltered_rna_bam', 'rna_fastq', 'dna_bam', 'hla_alleles'], sep: ",")
        .map{ row -> tuple(
            row.name, 
            file(row.unfiltered_rna_bam, checkIfExists: true),      // Convert to Path object
            file(row.rna_fastq, checkIfExists: true),      // Convert to Path object
            file(row.dna_fastq, checkIfExists: true),      // Convert to Path object
            file(row.hla_alleles, checkIfExists: true)      // Convert to Path object
        ) }
        //.map{ row -> tuple(row.name, row.normal_vcf, row.tumour_vcf, row.hg_bam, row.hla_alleles) }

    // Channel 1: name, RNA-seq BAM 
    bam_channel = input_data
        .map{ name, unfiltered_rna_bam, rna_fastq, dna_bam, hla_alleles -> tuple(name, unfiltered_rna_bam) }

    // Channel 2: name, RNA-seq FASTQ 
    fastq_channel = input_data
        .map{ name, unfiltered_rna_bam, rna_fastq, dna_bam, hla_alleles -> tuple(name, rna_fastq) }

    //Channel 3: name, DNA-seq BAM
    dna_channel = input_data
        .map{ name, unfiltered_rna_bam, rna_fastq, dna_bam, hla_alleles -> tuple(name, dna_bam) }

    //Channel 3: name, hla_alleles
    alleles_channel = input_data
        .map{ name, unfiltered_rna_bam, rna_fastq, dna_bam, hla_alleles -> tuple(name, hla_alleles)}


    //Workflow logic
    //Start by calling with multiple fusion gene calles 
    Jaffal(fastq_channel)

    Longgf(bam_channel)

    NanoFG(dna_channel)

    //find consensus between outputs, and filter out only consensus call from within Longgf (or Jaffal) output
    fusion_output_ch = Jaffal.out.join(Longgf.out)

    Consensus(fusion_output_ch)

    //Input into AGFusion
    AGFusion(Consensus.out)

    //Run pVACfuse on AGFusion output (optional)
    if (params.containsKey('include_pvacfuse') && params.include_pvacfuse) {
        pvacfuse_input_ch = AGFusion.out.join(alleles_channel)
        PvacFuse(pvacfuse_input_ch)
    }

    //tumour_channel.view()
    publish:
    jaffal                = Jaffal.out
    longgf                = Longgf.out
    nanofg                = NanoFG.out
    consensus             = Consensus.out
    agfusions             = AGFusion.out
    pvacfuse_neoag        = (params.containsKey('include_pvacfuse') && params.include_pvacfuse) ? PvacFuse.out : Channel.empty()
}

output {    

    jaffal              { path { name, jaffal           -> "${name}/jaffal" } }
    longgf              { path { name, longgf           -> "${name}/longgf" } }
    nanofg              { path { name, nanofg           -> "${name}/nanofg" } }
    consensus           { path { name, consensus        -> "${name}/consensus" } }
    agfusions           { path { name, agfusions        -> "${name}/agfusions" } }
    pvacfuse_neoag      { path { name, pvacfuse_neoag   -> "${name}/pvacfuse" } }
}

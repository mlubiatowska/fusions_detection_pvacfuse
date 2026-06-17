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
        .splitCsv(header: ['name', 'unfiltered_rna_bam', 'rna_fastq'], sep: ",")
        .map{ row -> tuple(
            row.name, 
            file(row.unfiltered_rna_bam, checkIfExists: true),      // Convert to Path object
            file(row.rna_fastq, checkIfExists: true)      // Convert to Path object
        ) }
        //.map{ row -> tuple(row.name, row.normal_vcf, row.tumour_vcf, row.hg_bam, row.hla_alleles) }

    // Channel 1: name, RNA-seq BAM 
    bam_channel = input_data
        .map{ name, unfiltered_rna_bam, rna_fastq -> tuple(name, unfiltered_rna_bam) }

    // Channel 2: name, RNA-seq FASTQ 
    fastq_channel = input_data
        .map{ name, unfiltered_rna_bam, rna_fastq -> tuple(name, rna_fastq) }

    //Workflow logic
    //Start by calling with multiple fusion gene calles 
    Jaffal(fastq_channel)

    Longgf(bam_channel)

    //find consensus between outputs, and filter out only consensus call from within Longgf (or Jaffal) output
    fusion_output_ch = Jaffal.out.join(Longgf.out)

    Consensus(fusion_output_ch)

    //Run pVACfuse on AGFusion output (optional)
    if (params.containsKey('include_pvacfuse') && params.include_pvacfuse) {

        //Input into AGFusion
        AGFusion(Consensus.out)

        //Input into pVACfuse
        pvacfuse_input_ch = AGFusion.out.join(alleles_channel)
        PvacFuse(pvacfuse_input_ch)
    }

    //tumour_channel.view()
    publish:
    jaffal                = Jaffal.out
    longgf                = Longgf.out
    consensus             = Consensus.out
    pvacfuse_neoag        = (params.containsKey('include_pvacfuse') && params.include_pvacfuse) ? PvacFuse.out : Channel.empty()
}

output {    

    jaffal              { path { name, jaffal           -> "${name}/jaffal" } }
    longgf              { path { name, longgf           -> "${name}/longgf" } }
    consensus           { path { name, consensus        -> "${name}/consensus" } }
    pvacfuse_neoag      { path { name, pvacfuse_neoag   -> "${name}/pvacfuse" } }
}

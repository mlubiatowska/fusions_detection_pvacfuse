#!/usr/bin/env nextflow
nextflow.enable.dsl     = 2
/*
* Pipeline parameters
*/


include { Jaffal } from './modules/jaffal.nf'
include { LonggfPrep } from './modules/longgf.nf'
include { Longgf } from './modules/longgf.nf'
include { Consensus } from './modules/consensus.nf'
include { AGFusion_Jaffal } from './modules/agfusions_jaffal.nf'
include { PvacFuse_Jaffal } from './modules/pvacfuse_jaffal.nf'
include { AGFusion_LongGF } from './modules/agfusions_longgf.nf'
include { PvacFuse_LongGF } from './modules/pvacfuse_longgf.nf'
//include { FilterAgfusion_Jaffal } from './modules/filter_agfusion_jaffal_proteinfa.nf'
include { FilterAgfusion_LongGF } from './modules/filter_agfusion_jaffal_proteinfa.nf'

workflow {
    main:
    // Parse input CSV and create two separate channels   
    input_data = Channel
        .fromPath(params.input)
        .splitCsv(header: ['name', 'unfiltered_rna_bam', 'rna_fastq', 'hla_alleles'], sep: ",")
        .map{ row -> tuple(
            row.name, 
            file(row.unfiltered_rna_bam, checkIfExists: true),      // Convert to Path object
            file(row.rna_fastq, checkIfExists: true),      // Convert to Path object
            file(row.hla_alleles, checkIfExists: true)      // Convert to Path object
        ) }
        //.map{ row -> tuple(row.name, row.normal_vcf, row.tumour_vcf, row.hg_bam, row.hla_alleles) }

    // Channel 1: name, RNA-seq BAM 
    bam_channel = input_data
        .map{ name, unfiltered_rna_bam, rna_fastq, hla_alleles -> tuple(name, unfiltered_rna_bam) }

    // Channel 2: name, RNA-seq FASTQ 
    fastq_channel = input_data
        .map{ name, unfiltered_rna_bam, rna_fastq, hla_alleles -> tuple(name, rna_fastq) }

    //Channel 3: name, hla_alleles
    alleles_channel = input_data
        .map{ name, unfiltered_rna_bam, rna_fastq, hla_alleles -> tuple(name, hla_alleles)}


    //Workflow logic
    //Start by calling with multiple fusion gene calles 
    Jaffal(fastq_channel)

    LonggfPrep(bam_channel)

    Longgf(LonggfPrep.out)

    //find consensus between outputs, and filter out only consensus call from within Longgf (or Jaffal) output
    fusion_output_ch = Jaffal.out.join(Longgf.out)

    Consensus(fusion_output_ch)

    //Run pVACfuse on AGFusion output (optional)
    if (params.containsKey('include_pvacfuse') && params.include_pvacfuse) {

        //Input into AGFusion and pvacfuse -for jaffal output 
        pvac_jaffal_ch = Consensus.out
            .map{ name, jaffal_consensus, jaffal_consensus_breakpoints, longgf_consensus_breakpoints -> tuple(name, jaffal_consensus_breakpoints)}
        AGFusion_Jaffal(pvac_jaffal_ch)
        //FilterAgfusion_Jaffal(AGFusion_Jaffal.out)
        jaffal_pvacfuse_input_ch = AGFusion_Jaffal.out.join(alleles_channel)
        PvacFuse_Jaffal(jaffal_pvacfuse_input_ch)

        //Input into AGFusion and pvacfuse -for longgf output
        pvac_longgf_ch = Consensus.out
            .map{ name, jaffal_consensus, jaffal_consensus_breakpoints, longgf_consensus_breakpoints -> tuple(name, longgf_consensus_breakpoints) }
        AGFusion_LongGF(pvac_longgf_ch)
        FilterAgfusion_LongGF(AGFusion_LongGF.out)
        longgf_pvacfuse_input_ch = FilterAgfusion_LongGF.out.join(alleles_channel)
        PvacFuse_LongGF(longgf_pvacfuse_input_ch)
        //testing if overlapping causes problems

        //filtered_ch = AGFusion_Jaffal.out.join(FilterAgfusion_LongGF.out)
        //pvacfuse_input_ch = filtered_ch.join(alleles_channel)
        //PvacFuse_Jaffal(pvacfuse_input_ch)
        
    }

    //tumour_channel.view()
    publish:
    jaffal                = Jaffal.out
    longgf                = Longgf.out
    consensus             = Consensus.out
    jaffal_agfusion       = AGFusion_Jaffal.out
    longgf_agfusion       = AGFusion_LongGF.out
    //filtered_jaffal_agfusion  = FilterAgfusion_Jaffal.out
    filtered_longgf_agfusion  = FilterAgfusion_LongGF.out
    jaffal_pvacfuse_neoag        = (params.containsKey('include_pvacfuse') && params.include_pvacfuse) ? PvacFuse_Jaffal.out : Channel.empty()
    longgf_pvacfuse_neoag        = (params.containsKey('include_pvacfuse') && params.include_pvacfuse) ? PvacFuse_LongGF.out : Channel.empty()
}

output {    

    jaffal              { path { name, jaffal                                   -> "${name}/jaffal" } }
    longgf              { path { name, longgf                                   -> "${name}/longgf" } }
    consensus           { path { name, jaffal_consensus, jaffal_consensus_breakpoints, longgf_consensus_breakpoints         -> "${name}/${params.consensus_outdir}" } }
    jaffal_agfusion     { path { name, agfusion                                 -> "${name}/agfusion" } }
    longgf_agfusion     { path { name, agfusion                                 -> "${name}/agfusion" } }
    //filtered_jaffal_agfusion    { path { name, agfusion_filtered, problematic_transcripts_report, missing_exons_report -> "${name}/agfusion" } }
    filtered_longgf_agfusion    { path { name, agfusion_filtered, problematic_transcripts_report, missing_exons_report -> "${name}/agfusion" } }
    jaffal_pvacfuse_neoag      { path { name, jaffal_pvacfuse_neoag                           -> "${name}/${params.pvacfuse_outdir}/jaffal" } }
    longgf_pvacfuse_neoag      { path { name, longgf_pvacfuse_neoag                           -> "${name}/${params.pvacfuse_outdir}/longgf" } } 
}
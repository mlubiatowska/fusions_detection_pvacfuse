# Fusion Gene consensus calling and neoantigen prediction
## Pipeline for fusion calling in myeloma samples
- In the future this will be followed by filtering and neoantigen prediciton from the consensus-called fusions.
- currently only working for fusion calling and consensus, using JAFFAL and LongGF
- requires local nextflow.config with correct pathways 

## Running the pipeline 
### Command
```
nextflow run https://github.com/mlubiatowska/fusions_detection_pvacfuse.git -r only_concordance \
    --input ${INPUT_FILE} \
    --outdir ${OUTPUT_DIR} \
	--include_pvacfuse
```
### Parameters:
```
 --input: csv file with data for all sampes (no header). Earch row should be 1 sample and must contain sample_name, unfiltered_rna_bam, rna_fastq, hla_alleles (this order)
 --outdir: path to output directory. Sub-derectories for each sample_name will be created there
 --include_pvacfuse: OPTIONAL argument to run AGFusion (fusion annotation) and pVACfuse (neoantigen prediction)
```

## Additional info about the pipeline 
### Tools 
Consensus called based on gene names and their order in fusion candidates. Fusion calling using:
- JAFFAL (long read version of JAFFA) v2.3,
- LongGF v0.1.2 with gencode.v29.chr_patch_hapl_scaff.annotation.gtf,
- NanoFG (no version - installed 30/03/2023, with dependencies including python Python/3.11, SAMtools/1.11, minimap2-2.30, LAST (from util-linux 2.32.1), NanoSV v1.2.4, PyVCF v0.6.8, pysam 0.23.3, wtdbg2.5, primer3 v 2.6.1). Reference used was /data/reference-data/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa and lastdb database was built based on it.
	- currently not included in the pipeline 
- AGFusion 
- pVACtools 

### INPUT: 
- RNAseq BAM and 
- RNAseq FASTQ files form , 
- HLA-LA output 
- * BAM files for WGS - currently not needed

### OUTPUT: 
- predicted fusions, 
- consensus fusions, 
- neoantigens predicted from the fusions.
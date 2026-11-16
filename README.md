# Fusion Gene consensus calling and neoantigen prediction
## Pipeline for fusion calling in myeloma samples

## Running the pipeline 
### Command
```
nextflow pull https://github.com/mlubiatowska/fusions_detection_pvacfuse.git -r main
nextflow run https://github.com/mlubiatowska/fusions_detection_pvacfuse.git -r main \
    --input sample_file_paths.csv \
    --outdir output \
    --min_read_count 3 \
    --breakpoint_tolerance 3 \
    --consensus_outdir consensus_ensembl115 \
    --nonmalignant_fusions_path fusions.csv \
    --include_pvacfuse \
    --pvacfuse_outdir pvacfuse \
    -c <path to custom config file>

```
### Parameters:
```
--input: csv file with data for all sampes (no header). Earch row should be 1 sample and must contain sample_name, unfiltered_rna_bam, rna_fastq, hla_alleles (this order) (required)
--outdir: path to output directory. Sub-derectories for each sample_name will be created there (required)
--min_read_count: minimum number of fusion gene supporting reads in JAFFAL and LongGF (testing with 5) (required)
--breakpoint_tolerance: maximum distance between breakpoint in JAFFAL and LongGF (testing with 5) (required)
--consensus_outdir: name of directories that contain consensus output (will be saved under directory with sample name) (required)
--nonmalignant_fusions_path: path to recurrent non-malignant fusions to be filtered out (downloaded from Babiceanu et al Nucleic Acids Res 2016) (required)
--include_pvacfuse: OPTIONAL argument to run AGFusion (fusion annotation) and pVACfuse (neoantigen prediction). Can be added after fusion prediction has been run using nextflow parameter 
--pvacfuse_outdir: OPTIONAL (required with option include_pvacfuse) name of directories that contain pvafuse output (will be saved under directory with sample name)
-c: path to nextflow.config
```

## Additional info about the pipeline 
### Tools 
Consensus called based on gene names and their order in fusion candidates. Fusion calling using:
- JAFFAL (long read version of JAFFA) v2.5 with GENCODE v49 database
- LongGF v0.1.2 with gencode.v49.chr_patch_hapl_scaff.annotation.gtf,
- AGFusion with Ensembl release 115 database (loaded with SciPy-bundle/2023.11-gfbf-2023b)
- pVACtools v6.0.5
- SAMtools/1.11 (for name-sorting the BAM files for LongGF input)
- R v4.5.2 with tidyverse v2.0.0

### INPUT: 
Input is a CSV file in which each row contains sample name and paths to required files in this order:
- Sample name,
- RNAseq BAM,
- RNAseq FASTQ,
- HLA-LA output

### OUTPUT: 
- JAFFAL and LongGF output files with fusion genes, 
- consensus fusions in JAFFAL and LongGF format (lines starting with ^GF for AGFUsion input):
    1. after filtering out the recurrent nonmalingant fusions ((jaffal/longgf)_consensys_NOnonmalignant.csv), JAFFAL calls clasified as 'PotentialTransSplicing' and 'PotentialReadThrough', fusions between 2 genes with the same names (includes IGH:IGH, IGL:IGL, IGK:IGK)
    2. above + only calls passing min_read_count filter and JAFFAL HighConfidence calls (exon-aligning calls, present in less than 20 GTEx tissues OR present in COSMIC, min 10kb distance, no mitochondrial genes)
    3. above + only calls passing breakpoint_tolerance filter
- AGFusion annotations and pVAFfuse neoantigens predicted from the fusions in JAFFAL and LongGF-format of filtered consensus calls
    - NOTE: JAFFAL and LongGF formats will not include exactly the same fusions, since the breakpoint called by each tool will differ. This is why AGFusion and pVACfuse is run on outputs in both formats, preserving all possible fusion neoantigen configuarations.
    - NOTE2: all outputs are saved, but only the 2nd one is passed further to AGFusion and pVACFuse 
        - this is because adding breakpoint agreement of 3bp filters out few genes, and the ones that are filtered out were shown to be true fusions. 

Pipeline for fusion calling in myeloma samples, followed by neoantigen prediciton from the consensus-called fusions.

Consensus called based on gene names and their order in fusion candidates. 
Fusion calling using:
- JAFFAL (long read version of JAFFA) v2.3, 
- LongGF v0.1.2 with gencode.v29.chr_patch_hapl_scaff.annotation.gtf, 
- NanoFG (no version - installed 30/03/2023, with dependencies including python Python/3.11, SAMtools/1.11, minimap2-2.30, LAST (from util-linux 2.32.1), NanoSV v1.2.4, PyVCF v0.6.8, pysam 0.23.3, wtdbg2.5, primer3 v 2.6.1). Reference used was /data/reference-data/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa and lastdb database was built based on it.  


INPUT: BAM and FASTQ files form RNAseq, and BAM files for WGS. OUTPUT: predicted fusions, consensus fusions, neoantigens predicred from the fusions.

# MitoMonitor
Mitochondrial genome assembly, annotation and taxonomic classification from short read metagenomic datasets

## Installation
mamba create -n MitoMonitor -c conda-forge -c bioconda snakemake mitos bwa samtools seqkit seqtk megahit minimap2 mamba r-bold r-readr r-parallel bioconductor-biostrings novoplasty

## MITOS database 

[download here](https://zenodo.org/record/2683856) and extract into root

## Reference mitochondrial genomes
Download target mitochondrial genomes from NCBI GenBank and extract to data/ref/

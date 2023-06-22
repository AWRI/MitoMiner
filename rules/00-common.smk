import pandas as pd


configfile: "config.yaml"


# read samples and fastq files
samples = pd.read_table("samples.tsv", dtype=str).set_index("sample", drop=False)
sampleList = samples.index.tolist()


# Wildcard constraints, required for getFq functions
wildcard_constraints:
    sample="|".join(samples.index)


# utility functions
def getFq1(wildcards):
    """Get R1 fastq reads for a sample."""
    fastqs = samples.loc[wildcards.sample].dropna()
    return fastqs.fq1

def getFq2(wildcards):
    """Get R2 fastq reads for a sample."""
    fastqs = samples.loc[wildcards.sample].dropna()
    return fastqs.fq2

rule mapReads:
    """
    Map to combined mitochondrial reference database
    """
    input:
        getFq1,
        getFq2
    output:
        mapDir + '{sample}.sorted.bam'
    params:
        config['ref']['fasta']
    threads:
        16
    log:
        errDir + '{sample}.mapReads.stderr'
    shell:
        ' bwa mem -t {threads} {params} {input} 2> {log} | samtools view -hb -F 4 | samtools sort -@ {threads} > {output}  ' 


rule getReads:
    """
    Map to combined mitochondrial reference database
    """
    input:
        fq1 = getFq1,
        fq2 = getFq2,
        bam = mapDir + '{sample}.sorted.bam'
    output:
        readNames = mapDir + '{sample}.readnames',
        R1 = mapDir + '{sample}_R1_mapped.fastq.gz',
        R2 = mapDir + '{sample}_R2_mapped.fastq.gz'
    params:
        config['ref']['fasta']
    threads:
        2
    log:
        errDir + '{sample}.mapReads.stderr'
    run:
        shell("samtools view {input.bam} | cut -f1 | sort | uniq > {output.readNames}")
        shell("seqtk subseq {input.fq1} {output.readNames} | bgzip > {output.R1}")
        shell("seqtk subseq {input.fq2} {output.readNames} | bgzip > {output.R2}")


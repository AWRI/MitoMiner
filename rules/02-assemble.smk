rule assemble:
    """
    assemble mitochondrial reads
    """
    input:
        R1 = mapDir + '{sample}_R1_mapped.fastq.gz',
        R2 = mapDir + '{sample}_R2_mapped.fastq.gz'
    output:
        file = assembleDir + '{sample}.contigs.fa'
    threads:
        4
    params:
        sample = '{sample}',
        out = assembleDir + '{sample}_megahit/'
    log:
        errDir + '{sample}.megahit.stderr'
    run:
        shell(" megahit -t {threads} -1 {input.R1} -2 {input.R2} -o {params.out} --out-prefix {params.sample} 2> {log}")
        shell("mv {params.out}/{params.sample}.contigs.fa {output}")


rule mapContigs:
    """
    map contigs to mito refs
    """
    input:
        assembly = assembleDir + '{sample}.contigs.fa'
    output:
        file = assembleDir + '{sample}.sorted.bam'
    threads:
        4
    params:
        config['ref']['fasta']
    log:
        errDir + '{sample}.minimap2.stderr'
    shell:
        'minimap2 -t {threads} -ax splice --cs {params} {input} 2> {log} | samtools view -hb -F 4 | samtools sort > {output}'

rule getContigs:
    """
    extract mitochondrial contigs 
    """
    input:
        assembly = assembleDir + '{sample}.contigs.fa',
        bam = assembleDir + '{sample}.sorted.bam'
    output:
        readNames = assembleDir + '{sample}.tigNames',
        tigs = outDir + '{sample}.fa'
    params:
        config['ref']['fasta']
    threads:
        1
    log:
        errDir + '{sample}.mapReads.stderr'
    run:
        shell("samtools view {input.bam} | cut -f1 | sort | uniq > {output.readNames}")
        shell("seqtk subseq {input.assembly} {output.readNames} > {output.tigs}")
        
        
        
rule annotateMito:
    """
    map contigs to mito refs
    """
    input:
        tigs = outDir + '{sample}.fa'
    output:
        file = outDir + '{sample}.mitoGenes.fasta',
        file2 = annotateDir + '{sample}.mitoGenes.fasta'
    threads:
        config['mitosCPU']
    params:
        sample = '{sample}',
        genCode = config['genCode']
    run:
        import queue, re
        from Bio import SeqIO
        from multiprocessing import Process, Lock, JoinableQueue
        
        fasta_sequences = SeqIO.parse(open(input.tigs),'fasta')
        outDir = str("data/processing/annotate/" + params.sample + "/" )
        if os.path.isdir(outDir):
                    shell(' rm -r {outDir} ')
        os.makedirs(outDir)

        # set up for multithreading and make queue
        mutex = Lock()
        q = JoinableQueue()
        
        def worker():
            """
            Function to itterate across fasta records and pass to mitos.
            """
            while True:
                fasta = q.get()
                if fasta is None:
                    break
                out = str("data/processing/annotate/" + params.sample + "/" + params.sample +  "_" + fasta.id + ".fasta")
                print(fasta.id)
                SeqIO.write(fasta, out, 'fasta')
                mitosOut = str("data/processing/annotate/" + params.sample + "/" + params.sample +  "_" + fasta.id + "_mitos")
                mitosLog = str("data/processing/annotate/" + params.sample + "/" + params.sample +  "_" + fasta.id + "_mitos.log")
                if os.path.isdir(mitosOut):
                    shell(' rm -r {mitosOut} ')
                    
                os.makedirs(mitosOut)
                shell(' runmitos.py -r mitos1-refdata/ -i {out} -o {mitosOut} -c {params.genCode} > {mitosLog}')
                
                fromFa = str("data/processing/annotate/" + params.sample + "/" + params.sample +  "_" + fasta.id + "_mitos/result.fas")
                toFa = str("data/processing/annotate/" + params.sample + "_" + fasta.id + "_results.fas")
                shell(' mv {fromFa} {toFa} ')
                q.task_done()
                
        def runMitos():
            """
            Function to queue jobs and run in multithread
            """
            
            if __name__ == 'snakemake.workflow':
                
                ## multithreading run
                # define threads and pass worker function
                jobs = []
                for n in range(threads):
                    j = Process(target=worker)
                    jobs.append(j)
                    j.start()
                
                # add to the queue and initiate processes
                # params added to only look up from local env
                for fasta in fasta_sequences:
                    q.put(fasta)    
                
                # wait for queue to empty    
                q.join()
                
                # purge the queue
                for i in range(threads):
                    q.put(None)
                
                # backup wait for threads
                for job in jobs:
                    job.join()
            else:
                return 
        
        runMitos()
        loc = str("data/processing/annotate/" + params.sample + "*fas")
        shell(' cat {loc} > {output.file}')
        shell(' cat {loc} > {output.file2}')
            
         
rule boldIdentification:
    """
    identify COX1 sequences with bold 
    """
    input:
        assembly = outDir + '{sample}.mitoGenes.fasta'
    output:
        file = outDir + '{sample}_boldIdentification.tsv'
    params:
        sample = '{sample}'
    threads:
        4
    log:
        errDir + '{sample}.boldID.stderr'
    shell:
        ' Rscript scripts/classifyBold.R {input.assembly} {params.sample} {threads} '
        

rule buildBlastDB:
    """
    make blast DB for reference 
    """
    input:
        config['ref']['fasta']
    output:
        config['ref']['fasta'] + '.nin'
    threads:
        1
    shell:
        ' makeblastdb -in {input} -dbtype nucl '
        
rule megaBLAST:
    """
    megaBLAST against references 
    """
    input:
        data = outDir + '{sample}.fa',
        ref = config['ref']['fasta'] + '.nin'
    output:
        outDir + '{sample}_blast.tsv'
    params:
        config['ref']['fasta']
    threads:
        4
    log:
        errDir + '{sample}.blast.stderr'
    shell:
        ' blastn -db {params} -query {input.data} -out {output} -outfmt "6 qseqid sseqid pident length evalue qcovs" -num_threads {threads}'


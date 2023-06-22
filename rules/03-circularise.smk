rule extractCox1:
    """
    Extract cox1 for each library 
    """
    input:
        data = outDir + '{sample}.mitoGenes.fasta'
    output:
        file = circularDir + '{sample}.cox1.fasta'
    threads:
        1
    shell:
        ' seqkit grep -nrip ".*cox1.*" {input.data} > {output.file} '
        
        
rule novoPlasty:
    """
    try and circularise each cox1 
    """
    input:
        data = circularDir + '{sample}.cox1.fasta',
        R1 = mapDir + '{sample}_R1_mapped.fastq.gz',
        R2 = mapDir + '{sample}_R2_mapped.fastq.gz'
    output:
        file = outDir + '{sample}.novoplasty.fasta'
    params:
        sample = '{sample}',
        NOVOplasty = config['NPexec']
    threads:
        1
    run:
        from Bio import SeqIO
        import re
        
        fasta_sequences = SeqIO.parse(open(input.data),'fasta')
        
        out = str("data/processing/circularise/" + params.sample + "/" ) 
        if os.path.isdir(out):
            shell('rm -r {out}')
        os.makedirs(out)
        
        for fa in fasta_sequences:
            print(len(fa))
            if not len(fa) > 150:
                print('breaking')
                continue
            else:    
                faN = re.sub(";", "", fa.id)
                outFAN = out + faN + '/'
                if os.path.isdir(outFAN):
                    shell('rm -r {outFAN}')
                os.makedirs(outFAN)
                fa.id = faN
                SeqIO.write(fa, outFAN + faN + '_cox1' + '.fasta', 'fasta')

            
            ##### make config file
                with open("data/processing/circularise/" + params.sample + "/" + faN + '/' 'NovoPlasty_config.txt', 'a') as config_file:
                    config_file.write('Project:\n')
                    config_file.write('-----------------------\n')
                    config_file.write('Project name          = ' + params.sample + faN + '\n')
                    config_file.write('Type                  = mito\n')
                    config_file.write('Genome Range          = 12000-25000\n')
                    config_file.write('K-mer                 = 33\n')
                    config_file.write('Max memory            = \n')
                    config_file.write('Extended log          = 0\n')
                    config_file.write('Save assembled reads  = yes\n')
                    config_file.write('Seed Input            = ' + out + faN + '/' + faN + '_cox1' + '.fasta' + '\n')
                    config_file.write('Extend seed directly  = no\n')
                    config_file.write('Reference sequence    = \n')
                    config_file.write('Variance detection    = \n')
                    config_file.write('Chloroplast sequence  = \n')
                    config_file.write('\n')
                    config_file.write('Dataset 1:\n')
                    config_file.write('-----------------------\n')
                    config_file.write('Read Length           = 151\n')
                    config_file.write('Insert size           = 300\n')
                    config_file.write('Platform              = illumina\n')
                    config_file.write('Single/Paired         = PE\n')
                    config_file.write('Combined reads        = \n')
                    config_file.write('Forward reads         = ' + input.R1 + '\n')
                    config_file.write('Reverse reads         = ' + input.R2 + '\n')
                    config_file.write('Store Hash            = \n')
                    config_file.write('\n')
                    config_file.write('Heteroplasmy:\n')
                    config_file.write('-----------------------\n')
                    config_file.write('MAF                   = \n')
                    config_file.write('HP exclude list       = \n')
                    config_file.write('PCR-free              = \n')
                    config_file.write('\n')
                    config_file.write('Optional:\n')
                    config_file.write('-----------------------\n')
                    config_file.write('Insert size auto      = yes\n')
                    config_file.write('Use Quality Scores    = no\n')
                    config_file.write('Output path           = ' + out + faN + '/' + '\n' )
                    
                shell('perl ' + params.NOVOplasty + ' -c ' "data/processing/circularise/" + params.sample + "/" + faN + '/' 'NovoPlasty_config.txt')
                

        shell(' echo "" | cat `find ' + out + ' -name "Circularized*"` | seqkit replace -p ".+" -r "CM_{{nr}}" > ' + output.file)

rule circBLAST:
    """
    megaBLAST circular against references 
    """
    input:
        data = outDir + '{sample}.novoplasty.fasta',
        ref = config['ref']['fasta'] + '.nin'
    output:
        outDir + '{sample}_circular_blast.tsv'
    params:
        config['ref']['fasta']
    threads:
        1
    log:
        errDir + '{sample}.blast.stderr'
    shell:
        ' blastn -db {params} -query {input.data} -out {output} -outfmt "6 qseqid sseqid pident length evalue qcovs" -num_threads {threads}'


rule out:
    """
    Map to combined mitochondrial reference database
    """
    input:
        oo = expand(outDir + '{sample}_circular_blast.tsv', sample=sampleList)
    output:
        outDir + 'counted_marker_kmers.tsv',
        outDir + 'test'
    threads:
        1
    log:
        errDir + 'out.stderr'
    shell:
        ' echo "cancel"'
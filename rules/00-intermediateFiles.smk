# directories
errDir = 'errorLogs/'
dataDir = 'data/'
processingDir = dataDir + 'processing/'
mapDir = processingDir + 'mapReads/'
assembleDir = processingDir + 'assembly/'
annotateDir = processingDir + 'annotate/'
circularDir = processingDir + 'circularise/'
outDir = 'out/'


finalOutput = [ outDir + 'counted_marker_kmers.tsv', outDir + 'test', expand(outDir + '{sample}_blast.tsv', sample=sampleList), expand(outDir + '{sample}_boldIdentification.tsv', sample=sampleList) ]


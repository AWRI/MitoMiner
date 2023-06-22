library(Biostrings)
library(bold)
library(parallel)
library(readr)


args <- commandArgs(trailingOnly=TRUE)

fa <- args[1]
sample <- args[2]
threads <- as.numeric(args[3])

dna <- readDNAStringSet(fa)

dna <- dna[grepl("cox1", names(dna))]

query <- mclapply(1:length(dna), mc.cores=10 , function(x){
  
  sDNA <- dna[x]
  
  df <- bold_identify(as.character(sDNA), db = "COX1")[[1]]
  
  if(length(df))
    df$contigName <- names(sDNA)
  df
})

query <- do.call(rbind, query)
query <- as.data.frame(query)

write_tsv(query, paste0("out/", sample, "_boldIdentification.tsv"))

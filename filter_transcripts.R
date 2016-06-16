# Start writing to an output file
logfile <- file(paste("logs/filter_transcripts.R.log",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "txt", sep = "."))
sink(logfile)
sink(logfile, type="message")
source("code/functions.R")
library(data.table)
library(Biostrings)

# create list of final transcript ID's actually used in the finial analysis so it can be used in previous steps to filter transcripts so as to not process more than required

observed.transcripts <- readRDS("data/observed.transcripts.rds")

print(paste(nrow(observed.transcripts),"observed unique transcripts"))

# Assign gene to transcripts
library(biomaRt)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", dataset="hsapiens_gene_ensembl")
updated.annotations <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id",'gene_biotype','transcript_biotype','hgnc_symbol','external_gene_name'), mart = ensembl, uniqueRows=TRUE)
updated.annotations <- data.table(updated.annotations)
#ensembl.map[,.N,by=transcript_biotype][order(-N)]
#ensembl.map[transcript_biotype=="protein_coding"][hgnc_symbol!=external_gene_name][sample(.N,99)]
setnames(updated.annotations, make.names(c("Ensembl Gene ID","Ensembl Transcript ID","Gene type","Transcript type","HGNC Symbol","Associated Gene Name")))

setkey(updated.annotations,Ensembl.Transcript.ID)

# setkey(observed_variants,transcript.id)
# updated.annotations[observed_variants][!is.na(Ensembl.Gene.ID) & Transcript.type=="protein_coding"]

# all rows in observed transcripts that have Gene.IDs in GRCh38 and are protein coding
observed.transcripts <- updated.annotations[observed.transcripts][!is.na(Ensembl.Gene.ID) & Transcript.type=="protein_coding"]
print(paste(nrow(observed.transcripts),"observed unique transcripts that still exist and are protein coding"))

# Get transcript lengths
fasta <- readDNAStringSet("data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz")

# Which sequences are not multiple of 3
multiple.of.3 <- width(fasta) %% 3L == 0L
# Which sequences are not ACTG only
clean.alphabet <- alphabetFrequency(fasta, baseOnly=TRUE)[,'other'] == 0

# Remove uncalculatable sequences (104,763 - 20,870 = 83,893 left)
fasta <- fasta[multiple.of.3 & clean.alphabet]

transcript.lengths <- data.table("transcript.id"=substr(names(fasta),0,15),"cds.length"=lengths(fasta))
setkey(transcript.lengths,transcript.id)

# Inner Join, all rows with ok transcript length & alphabet, with transcript that still exist and were observed
ok.transcripts <- observed.transcripts[transcript.lengths, nomatch=0]
print(paste(nrow(ok.transcripts),"transcripts left after removing not multiple of 3 or N in coding sequence"))
setkey(ok.transcripts,Ensembl.Gene.ID)

# For each gene, find max cds length
max.length.by.gene <- ok.transcripts[,.(max.cds = max(cds.length)),by=Ensembl.Gene.ID]
print(paste(nrow(max.length.by.gene),"unique genes mapped to ok transcripts"))

# For all rows in max.length.by.gene, add cds values where cds_length=max(cds_length) by gene, for each gene (multiples as sometimes multiple transcript have max cds length)
setkey(max.length.by.gene,Ensembl.Gene.ID,max.cds)
setkey(ok.transcripts,Ensembl.Gene.ID,cds.length)
transcript.per.gene <- ok.transcripts[max.length.by.gene,]
# 23,685

paste(sum(duplicated(transcript.per.gene)),"duplicate entries (same gene and cds length) after getting dNdS from canonical (longest) transcript")

# set key as gene id
setkey(transcript.per.gene,Ensembl.Gene.ID)

# Chose random transcript if two transcripts match to same gene and same length
unique <- unique(transcript.per.gene)
# 18,829

print("Genes with duplicate gene IDs")
unique[Associated.Gene.Name %in% unique[duplicated(unique)]$Associated.Gene.Name]

# Still some duplicate gene names as two ensembl gene IDs have the same name
setkey(unique,Associated.Gene.Name)
unique <- unique(unique)
print(paste(nrow(unique),"unique gene IDs mapped to transcripts"))

archive.file("data/final.transcript.list.rds")
saveRDS(unique,"data/final.transcript.list.rds",compress=FALSE)

sessionInfo()

sink(type="message")
sink()
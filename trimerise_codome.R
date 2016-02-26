# Start writing to an output file
logfile <- file(paste("logs/trimerise_codome.R.log",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "txt", sep = "."))
sink(logfile)
sink(logfile, type="message")

source("code/functions.R")

library(data.table)
print("Creating codome data")
print(Sys.time())
library(Biostrings)
print("Reading in reference genome")
fasta <- readDNAStringSet("data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz")

# Which sequences are not multiple of 3
multiple.of.3 <- width(fasta) %% 3L == 0L
# Which sequences are not ACTG only
clean.alphabet <- alphabetFrequency(fasta, baseOnly=TRUE)[,'other'] == 0

# Remove uncalculatable sequences (104,763 - 20,870 = 83,893 left)
fasta <- fasta[multiple.of.3 & clean.alphabet]

unique <- readRDS("data/final.transcript.list.rds")

# extract each trimer in codome with transcript id and original codon
load.transcripts <- function(x){
metadata <- substr(x,0,15) # transcriptID
# as.list(strsplit(x, "[: ]")[[1]])
if (metadata %in% unique){ # check if will be used in final analysis or not
transcript <- stringset[[x]]
length <- length(transcript)
views <- Views(transcript,seq(1,length-2,1),width=3)

# All trimers in transcript
transcript.table <- data.table("base_motif"=as.character(views),"transcript"=as.character(metadata))

# List of original codons
transcript.codons <- as.character(Views(transcript,seq(1,length(transcript)-2,3),width=3))

# position of middle nucleotie in original codon
codon.position <- rep(c("1","2","3"),length(transcript.codons))
# remove first and last
codon.position <- codon.position[c(-1,-length(codon.position))]
transcript.table$codon.position <- codon.position

# which was parent codon of middle nucleotide
original.codons <- rep(transcript.codons,each=3)
# Remove first and last codon
original.codons <- original.codons[c(-1,-length(original.codons))]
transcript.table$original.codons <- original.codons
return(transcript.table)
}else{}
}


stringset <- fasta
test <- lapply(names(stringset),load.transcripts)
# 11.5sec for 1000, *83 16 min without removing unused
# 480.202 - 8 mins (for 23k)

codome <- rbindlist(test)

# makes subsetting faster
setkey(codome,base_motif)

archive.file("data/codome.rds")
saveRDS(codome,"data/codome.rds",compress=FALSE)
print("Finished creating codome data")
print(Sys.time())


# standardise base motif (middle nt always pyrimidine)
standardise <- function(motif){
if (substr(motif,2,2) %in% c("G","A")){
	standardised.base_motif <- toString(reverseComplement(DNAString(motif)))
}else{
	standardised.base_motif <- motif
}
return(standardised.base_motif)
}

# count each trimer
motif.counts <- codome[,.N,by=base_motif]

# create conversion table 64 -> 32 standardised
standard_motif <- sapply(unique(codome$base_motif),standardise)
motif.names <- data.table("base_motif"=names(standard_motif),standard_motif)

# add standard names to counts
motif.counts <- motif.counts[motif.names]

standard.motif.counts <- motif.counts[,.("coding.trimer.counts"=sum(N)),by=standard_motif]

setnames(standard.motif.counts,"standard_motif","base_motif")

archive.file("data/coding.trimer.counts.rds")
saveRDS(standard.motif.counts,"data/coding.trimer.counts.rds",compress=FALSE)

sessionInfo()

sink(type="message")
sink()
source("code/functions.R")
# Start writing to an output file
logfile <- file(paste("logs/calculate_expected_variants.R.log",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "txt", sep = "."))
sink(logfile)
sink(logfile, type="message")
print(Sys.time())

library("Biostrings")
library(data.table)

motif.probabilities <- readRDS("data/motif.probabilities.rds")

genetic.code <- c("TTT"="F", "TTC"="F", "TTA"="L", "TTG"="L",
       "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S",
       "TAT"="Y", "TAC"="Y", "TAA"="STOP", "TAG"="STOP",
       "TGT"="C", "TGC"="C", "TGA"="STOP", "TGG"="W",
       "CTT"="L", "CTC"="L", "CTA"="L", "CTG"="L",
       "CCT"="P", "CCC"="P", "CCA"="P", "CCG"="P",
       "CAT"="H", "CAC"="H", "CAA"="Q", "CAG"="Q",
       "CGT"="R", "CGC"="R", "CGA"="R", "CGG"="R",
       "ATT"="I", "ATC"="I", "ATA"="I", "ATG"="M",
       "ACT"="T", "ACC"="T", "ACA"="T", "ACG"="T",
       "AAT"="N", "AAC"="N", "AAA"="K", "AAG"="K",
       "AGT"="S", "AGC"="S", "AGA"="R", "AGG"="R",
       "GTT"="V", "GTC"="V", "GTA"="V", "GTG"="V",
       "GCT"="A", "GCC"="A", "GCA"="A", "GCG"="A",
       "GAT"="D", "GAC"="D", "GAA"="E", "GAG"="E",
       "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G")

is.nonsynon <- function(codon1,codon2) {genetic.code[codon1] != genetic.code[codon2]}


# All possible 5mers
letters = c("A","T","G","C")
five.mers <- do.call(paste, c(expand.grid(letters, letters, letters, letters, letters), list(sep='')))
stopifnot(length(five.mers) == 1024)


calculate.fivemer.probability <- function(){
fivemer.probability.calc = data.table()
for (fivemer in (five.mers)){

fivemer.probability = 0
fivemer.probability.dt <- data.table()

for (codon.position in (2:4)){

	purine = FALSE
	
	# Get original trinucleotide in context for counting
	if (substr(fivemer, codon.position, codon.position) %in% c("G","A")){
		purine = TRUE
		original.trinucleotide = toString(reverseComplement(subseq(DNAString(fivemer), start=codon.position-1, end=codon.position+1)))
	} else {original.trinucleotide = substr(fivemer,codon.position-1,codon.position+1)}

codon.probability = 0

	#
	for (new.nt in letters[letters != substr(fivemer, codon.position, codon.position)]){		
		# Substitute new letter in codon.position
		m.fivemer <- fivemer
		substr(m.fivemer, codon.position, codon.position) <- new.nt


		# check if mutation is nonsynon or synon
		nonsynon = is.nonsynon(substr(m.fivemer,2,4),substr(fivemer,2,4))
		synon = 1 - nonsynon
		
		# Get mutation in context for counting
		if (purine == TRUE){
			new.trinucleotide = toString(reverseComplement(subseq(DNAString(m.fivemer), start=codon.position-1, end=codon.position+1)))
		} else{new.trinucleotide = substr(m.fivemer,codon.position-1,codon.position+1)}


		# Double check mutated base is referenced as pyrimidine
		stopifnot(substr(original.trinucleotide,2,2) %in% c("C","T"))
		
		# Format mutation as NN N.N
		mutation.in.context = paste(substr(original.trinucleotide,2,2),substr(new.trinucleotide,2,2)," ",substr(original.trinucleotide,1,1),".",substr(original.trinucleotide,3,3),sep='')
		
		codon.probability = 0
		# Outer product(list of mutation probabilities (one per project) * (synon,nosynon)) -> two lists (nonsynon prob & synon prob for each project)
		codon.probability <- motif.probabilities[mutation==mutation.in.context,mutation.probability] %o% c("nonsynon.prob"=nonsynon,"synon.prob"=synon)
		
		# For each site and letter add up probabilities
		fivemer.probability <- fivemer.probability + codon.probability
		} # each letter

		} # each codon.position
		
		# Convert to data table to join
		fivemer.probability.dt <- data.table(fivemer.probability)
		# add project codes
		fivemer.probability.dt$project_code <- motif.probabilities[mutation==mutation.in.context,project_code]
		# add fivemer
		fivemer.probability.dt$fivemer <- fivemer

	  
		# Add next fivemer set
	  fivemer.probability.calc <-  rbindlist(list(fivemer.probability.calc,fivemer.probability.dt))

	print(fivemer)
} # each fivemer
return(fivemer.probability.calc)
} # close function

fivemer.probabilities <- calculate.fivemer.probability()
# for some reason TAA is added to end of first two column names
names(fivemer.probabilities) <- c("nonsynon.prob","synon.prob","project_code","fivemer")


print("Reading in reference genome")
fasta <- readDNAStringSet("data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz")

# Sum over all projects
fivemer.probabilities.sum <- fivemer.probabilities[,.(nonsynon.probability=sum(nonsynon.prob),synon.probability=sum(synon.prob)),by=fivemer]
setkey(fivemer.probabilities.sum,fivemer)
unique <- readRDS("data/final.transcript.list.rds")
setkey(unique,Ensembl.Transcript.ID)

# check if transcript will be used in final analysis or not (reduce loops from 83893 to 18803 i.e. 1/5th)
fasta.names <- data.table(names(fasta))
fasta.names$transcript.id <- substr(fasta.names$V1,0,15)
setkey(fasta.names,transcript.id)
names <- fasta.names[unique]
# remove transcripts not in fasta (ESNTR etc)
print("transcripts not in fasta reference genome")
names[is.na(V1)]
names <- names[!is.na(V1)]


count.nonsynon <- function(x) {
transcript <- fasta[[x]]
# add N to start and beginning
n <- DNAString("NN")
subseq(n, start=2, end=1) <- transcript

# overlapping 5mers with codons in centre
views <- Views(n, start=seq(1,length(n)-3,3), width=5)

#SLOW
view.chr <- data.table(fivemer=as.character(views))
setkey(view.chr,fivemer)

# Sum always returns two NA values due to first and last fivemer having N at start/beginning
result <- fivemer.probabilities.sum[view.chr][,.("nonsynon.probability"=sum(nonsynon.probability,na.rm=TRUE),"synon.probability"=sum(synon.probability,na.rm=TRUE))]
result[,attributes:= x]
result[,cds.length:= length(transcript)]
return(result)

}

print(Sys.time())
print("Summing values per fivemer")
nonsynon.count <- lapply(names$V1,count.nonsynon)
# 6.5 mins

print(Sys.time())
print("Tidying and saving data")
nonsynon.count.dt <- rbindlist(nonsynon.count)

# Split attributes
nonsynon.count.dt[, c("transcript.id","cds.type","chromosome","gene.id","gene.biotype","transcript.biotype") := tstrsplit(attributes, " ", fixed=TRUE)]
# remove old attributes
nonsynon.count.dt$attributes <- NULL

# Remove variable names from entries
cnames <- names(nonsynon.count.dt)[4:9]
nonsynon.count.dt[, cnames := lapply(nonsynon.count.dt[,cnames,with=FALSE],function(x){gsub("[a-z]+_?[a-z]+:", "",x)}),with=FALSE]

# Split variant location field
nonsynon.count.dt[, c("reference.genome","chromosome","chromosome.start","chromosome.stop","strand") := tstrsplit(chromosome, ":", fixed=TRUE)]

print(paste(nrow(nonsynon.count.dt),"rows in final table"))

# remove useless (refgenome,strand) and untrustworthy (cds.type,gene.id,gene/transcript biotype) annotations 
nonsynon.count.dt$cds.type <- NULL
nonsynon.count.dt$gene.id <- NULL
nonsynon.count.dt$gene.biotype <- NULL
nonsynon.count.dt$transcript.biotype <- NULL
nonsynon.count.dt$strand <- NULL
nonsynon.count.dt$reference.genome <- NULL

# Add gene.name, mappability, cds length, gene.id
setkey(unique,Ensembl.Transcript.ID)
setkey(nonsynon.count.dt,transcript.id)

nonsynon.count.dt <- nonsynon.count.dt[unique[,.(Ensembl.Transcript.ID,"Ensembl.Gene.ID"=gene.id,mappability,Associated.Gene.Name)]]

# Save
archive.file("data/expected_variants_per_transcript.tsv")
write.table(nonsynon.count.dt, "data/expected_variants_per_transcript.tsv", sep="\t",row.names=FALSE,quote=FALSE)


# Quality Control Checks:
expected.variants <- readRDS("data/expected_variants_per_transcript.rds")
old.expected <- fread("data/archive/N_per_transcript_clean.tsv")
summary(expected.variants)
str(expected.variants)
expected.variants[nonsynon.probability>5,]

library(ggplot2)

# Set catagoricals as factors
expected.variants[, c("reference.genome","chromosome","strand","gene.biotype","transcript.biotype","cds.type") := lapply(expected.variants[,c("reference.genome","chromosome","strand","gene.biotype","transcript.biotype","cds.type"),with=FALSE],as.factor),with=FALSE]

# Correlation between old and new nonsynon values
setnames(old.expected,"transcript","transcript.id")
setkey(old.expected,transcript.id)
setnames(expected.variants,"Ensembl.Transcript.ID","transcript.id")
setkey(expected.variants,transcript.id)
merged <- merge(expected.variants,old.expected,all=TRUE)
merged$old.new.ratio <- merged$nonsynon.probability / merged$nonsynonymous_sites

pdf(paste("results/QC.nonsynon-vs-cds.length",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "pdf", sep = "."),onefile = TRUE)
plot(expected.variants$cds.length,expected.variants$nonsynon.probability,xlim=c(0,7000),ylim=
c(0,1.6))
new <- lm(expected.variants$nonsynon.probability~expected.variants$cds.length)
abline(new, col="red")
text(x=4000,y=0.2,label=paste("Intercept:",signif(new[[1]][1],3),"\nGradient:",signif(new[[1]][2],3)))

plot(old.expected$cds_length,old.expected$nonsynonymous_sites,xlim=c(0,7000),ylim=c(0,6000))
old <- lm(old.expected$nonsynonymous_sites~old.expected$cds_length)
abline(old, col="red")
text(x=4000,y=500,label=paste("Intercept:",signif(old[[1]][1],3),"\nGradient:",signif(old[[1]][2],3)))


plot(merged$nonsynonymous_sites,merged$nonsynon.probability)
#,col= ifelse(merged$old.new.ratio >= 0.00000161, "red","black")
plot(merged$nonsynonymous_sites,merged$nonsynon.probability,xlim=c(0,15000),ylim=c(0,4),col=rgb(0,0,0,alpha=0.2))
plot(merged$cds_length-merged$nonsynonymous_sites,merged$synon.probability,xlim=c(0,4000),ylim=c(0,1.2),col=rgb(0,0,0,alpha=0.2))

dev.off()

#merged[order(-old.new.ratio)][is.na(old.new.ratio)==FALSE][nonsynonymous_sites>12000]

sessionInfo()

sink(type="message")
sink()
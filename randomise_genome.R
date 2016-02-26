library(data.table)
library(Biostrings) # to reverse complement
library(argparser)
p <- arg_parser("Create randomised genome cohort")
p <- add_argument(p, "jobNumber", help = "Job number ")
cArgs <- parse_args(p)

print(paste("loading codome:",Sys.time()))
codome <- readRDS("data/codome.rds")
# 1.7GB, ~25 seconds to load, pre-keyed
print(paste("finished loading codome:",Sys.time()))

# Load motif probabilities
motif.probabilities <- readRDS("data/motif.probabilities.rds")

# Pretend all mutations came from one exome
motif.probabilities$single.exome.mutation.probaility <- motif.probabilities$mutation_count / motif.probabilities$coding.trimer.counts

# probability of base_motif mutating (to anything), over all patients in group
base.mutation.probability <- motif.probabilities[,.("base.mutation.probability"=sum(single.exome.mutation.probaility)),by=.(project_code,base_motif)]

# Relative probability of each mutation given base motif
intermotif.prob <- motif.probabilities[,.(mutation,"intermotif.prob"=single.exome.mutation.probaility/sum(single.exome.mutation.probaility)),by=.(project_code,base_motif)]
# Some are NaN (0/0)

# All possible trimers
nucleotides <- c("A","C","G","T") # alphabetical order!
trimers <- do.call(paste, c(expand.grid(nucleotides, nucleotides, nucleotides), list(sep='')))


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

projects = unique(motif.probabilities$project_code)


randomise.mutations <- function(){
	vcf <- data.frame()

for (project in projects){ # x55
	print(project)

for (motif in trimers){ # x64

#change to pyrimidine
if (substr(motif,2,2) %in% c("G","A")){
	standardised.base_motif <- toString(reverseComplement(DNAString(motif)))
}else{
	standardised.base_motif <- motif
}

# lookup base mutaiton prob
base.motif.probability <- base.mutation.probability[project_code==project & base_motif==standardised.base_motif,base.mutation.probability]

if (base.motif.probability == 0){next} # skip if no mutations to be chosen/made

# don't need the whole table - much faster (x30)
subset <- codome[base_motif==motif]

motif.count <- nrow(subset)
# toal number of mutations for this base_motif
length <- round(motif.count*base.motif.probability)

if (length == 0){next} # skip if no mutations to be chosen/made

# select new positions for mutations
mutation.positions <- sample(motif.count,size=length,replace=TRUE)

# select only mutated bases
subset <- subset[mutation.positions]

standardised.old.base <- substr(standardised.base_motif, 2, 2)
standardised.alt.bases <- nucleotides[-which(nucleotides==standardised.old.base)]

# intermotif probability
mutation.weightings <- intermotif.prob[mutation %in% paste0(standardised.old.base,standardised.alt.bases," ",substr(standardised.base_motif,1,1),".",substr(standardised.base_motif,3,3)) & project_code==project]$intermotif.prob

old.base <- substr(motif, 2, 2)
alt.bases <- nucleotides[-which(nucleotides==old.base)]


# relpace mutated bases with new nt
subset[,new.base :=sample(alt.bases,length,replace=TRUE,prob=mutation.weightings)]

# make a copy of codon to modify
subset$new.codon <- subset$original.codons

# Replace old base with new
substr(subset$new.codon, subset[,codon.position], subset[,codon.position]) <- subset$new.base

# annotate with S or N
subset$is.nonsynon <- is.nonsynon(subset$original.codons, subset$new.codon)

subset$project_code <- project

dflist = list(vcf,subset)
vcf <- rbindlist(dflist)
}
}

return(vcf)
}

vcf <- randomise.mutations()
# 206 seconds 3.4 mins

print("Summary of result:")
summary(vcf)

print("Structure of result:")
str(vcf)

# total mutations observed from ICGC data
print(paste(motif.probabilities[,sum(mutation_count)],"mutations in original ICGC data"))

# compare mutation number by project
print("mutations by project:")
data.table(vcf[,.("simulated"=.N),by=project_code],motif.probabilities[,.("observed"=sum(mutation_count)),by=project_code])

print("motifs mutated:")
vcf[,.("simulated"=.N),by=base_motif]

# exome length used
print(paste(sum(as.numeric(motif.probabilities$coding.trimer.counts))/(3*55),"nt exome length in motif probabilities"))

print(paste(nrow(codome),"nt exome length (codome)"))

# Table of synon and nonsynon counts by transcript
counts <- vcf[,.("N"=sum(is.nonsynon),"S"=.N-sum(is.nonsynon)),by=transcript]

transcripts <- data.table("transcript"=unique(codome$transcript))

setkey(counts,transcript)
counts <- counts[transcripts][order(transcript)]

counts$job <- cArgs$jobNumber

write.table(counts,paste0("data/randomised_genomes/",cArgs$jobNumber,".tsv"),sep="\t",row.names=FALSE,quote=FALSE)

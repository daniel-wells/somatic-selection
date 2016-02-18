# Start writing to an output file
logfile <- file(paste("logs/calculate_dNdS.R.log",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "txt", sep = "."))
sink(logfile)
sink(logfile, type="message")

source("code/functions.R")
library(data.table)
observed_variants <- fread("data/observed_variants_by_transcript.filtered.tsv",sep="|", header=FALSE)

setnames(observed_variants, c("gene.name", "gene.id", "strand", "transcript.name", "transcript.id", "protein.id", "variant.class", "nt.mutation", "aa.mutation", "projects","donors.affected","mutation","project.count","tested.donors"))
setkey(observed_variants, transcript.id)

# Check number of each variant.class
observed_variants[,.N,by=variant.class][order(-N)]

synon.count <- observed_variants[variant.class=="synonymous_variant",list(S=sum(donors.affected)), by=list(transcript.id)]

nonsynon.count <- observed_variants[variant.class %in% c("missense_variant","frameshift_variant","disruptive_inframe_deletion","disruptive_inframe_insertion","inframe_deletion","inframe_insertion","start_lost","stop_lost","stop_gained"), list(N=sum(donors.affected)), by=list(transcript.id)]

setkey(synon.count,transcript.id)
setkey(nonsynon.count,transcript.id)

# Full outer join
counts <- merge(synon.count,nonsynon.count, all=TRUE)
setkey(counts, transcript.id)
print(paste(nrow(counts[is.na(S),]),"transcripts with synon count of 0"))
print(paste(nrow(counts[is.na(N),]),"transcripts with nonsynon count of 0"))

# No rows with both N and S NA - i.e. all NA's are actually Zeros - gene must have been sequenced if one of the two catagories are >0
# Non sequenced genes are just not in this list as it's generated from variant observed list
print(paste(nrow(counts[is.na(N) & is.na(S),]),"transcripts with both nonsynon and synon count of 0"))

# Therefore set NA's to 0
counts$N[is.na(counts$N)] <- 0
counts$S[is.na(counts$S)] <- 0

expected_variants <- fread("data/expected_variants_per_transcript.tsv", header=TRUE)
setkey(expected_variants, transcript.id)

print(paste(nrow(expected_variants),"expected variants per transcript loaded"))

updated.annotations <- fread(input = 'zcat < data/raw/mart_export.txt.gz')
setnames(updated.annotations, make.names(names(updated.annotations)))
setkey(updated.annotations, "Ensembl.Transcript.ID")

# All rows in expected_variants + any rows from updated.annotations which match
expected_variants <- updated.annotations[expected_variants,]

# Remove non protein coding transcripts AND transcripts which are not present in GRCh38
expected_variants <- expected_variants[Transcript.type=="protein_coding",]

print(paste(nrow(expected_variants),"expected variants per transcript after filtering"))

# All rows in counts with rows from expected_variants which match
expected_variants <- expected_variants[counts,nomatch=0]
# some transcripts with observed variations do not have calculated nonsynon sites due to N or not multiple of 3, generating NA in e.g. chromosome column as these transcripts were not passed through to N per transcript. Inner Join
print(paste(nrow(expected_variants),"expected variants per transcript with measured actual variants"))

# Calculate dNdS
expected_variants$dS <- expected_variants$S / expected_variants$synon.probability
expected_variants$dN <- expected_variants$N / expected_variants$nonsynon.probability
expected_variants$dNdS <- expected_variants$dN / expected_variants$dS

print(paste(nrow(expected_variants[N==0]),"transcripts have N==0 generating a dNdS of 0"))
print(paste(nrow(expected_variants[S==0]),"transcripts have S==0 generating a dNdS of Inf"))
print(paste(nrow(expected_variants[S==0 & N==0]),"transcripts have S and N ==0"))
print(paste(nrow(expected_variants[is.na(dS)]),"transcripts have dS as NaN generating a dNdS of NaN:"))
expected_variants[synon.probability==0 | nonsynon.probability==0]

print(paste(nrow(expected_variants),"transcripts to be written to file"))


# Remove dNdS NA, NaN and Infinite rows
expected_variants <- expected_variants[is.na(dNdS)==FALSE & is.finite(dNdS)==TRUE,]
sessionInfo()

write.table(expected_variants,paste("data/dNdS_by_transcript",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "tsv", sep = ".") , sep="\t",row.names=FALSE,quote=FALSE)sink(type="message")
sink(type="message")
sink()
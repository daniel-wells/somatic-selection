# Start writing to an output file
logfile <- file(paste("logs/calculate_dNdS.R.log",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "txt", sep = "."))
sink(logfile)
sink(logfile, type="message")

source("code/functions.R")
library(data.table)

observed_variants <- readRDS("data/coding.mutations.rds")

# "gene.id" -> "gene_affected"
# "transcript.id" -> "transcript_affected"
# "variant.class" -> "consequence_type"
setnames(observed_variants,c("icgc_mutation_id","icgc_donor_id","project_code","chromosome","chromosome_start","chromosome_end","chromosome_strand","mutation_type","reference_genome_allele","mutated_from_allele","mutated_to_allele","variant.class","aa_mutation","cds_mutation","gene.id","transcript.id","sequencing_strategy"))

# columns removed in load_mutations.R creating duplicate rows
setkey(observed_variants,icgc_mutation_id,icgc_donor_id,transcript.id)
observed_variants <- unique(observed_variants)

# which project/primary site has the most mutations?
project.info <- fread("data/raw/ICGC_projects.tsv")
setnames(project.info, "Project Code", "project_code")
setnames(project.info, "Primary Site", "primary_site")
project.groupings <- project.info[,.(project_code,primary_site)]

project.groupings[order(primary_site)]

setkey(project.groupings,project_code)
setkey(observed_variants,project_code)

observed_variants <- project.groupings[observed_variants]

observed_variants[,.N,by=project_code][order(-N)]
observed_variants[,.N,by=primary_site][order(-N)]

setkey(observed_variants,icgc_mutation_id,icgc_donor_id,transcript.id)
# Check number of each variant.class
observed_variants[,.N,by=variant.class][order(-N)]

synon.count <- observed_variants[variant.class=="synonymous_variant",.("S"=.N), by=list(transcript.id)]

nonsynon.count <- observed_variants[variant.class %in% c("missense_variant","frameshift_variant","disruptive_inframe_deletion","disruptive_inframe_insertion","inframe_deletion","inframe_insertion","start_lost","stop_lost","stop_gained"), .("N"=.N), by=list(transcript.id)]

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

# Load expected variants per transcript
expected_variants <- fread("data/expected_variants_per_transcript.tsv", header=TRUE)
setkey(expected_variants, transcript.id)

print(paste(nrow(expected_variants),"expected variants per transcript loaded"))

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

# Remove dNdS NA, NaN
expected_variants <- expected_variants[is.na(dNdS)==FALSE,]
print(paste(nrow(expected_variants),"transcripts to be written to file"))

archive.file("data/dNdS_by_transcript.tsv")
write.table(expected_variants, "data/dNdS_by_transcript.tsv", sep="\t", row.names=FALSE, quote=FALSE)

sessionInfo()

sink(type="message")
sink()
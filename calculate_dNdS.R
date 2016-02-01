library(data.table)
observed_variants <- fread("data/observed_variants_by_transcript.filtered.tsv",sep="|", header=FALSE)

setnames(observed_variants, c("gene.name", "gene.id", "strand", "transcript.name", "transcript.id", "protein.id", "variant.class", "nt.mutation", "aa.mutation", "projects","donors.affected","mutation","project.count","tested.donors"))
setkey(observed_variants, transcript.id)

# Check number of each variant.class
observed_variants[,.N,by=variant.class][order(-N)]

# Remove exon variants
observed_variants <- observed_variants[variant.class!="exon_variant",]

# Calculate and add randomised versions of annotations

randomise <- function(){
  type <- c("synonymous_variant","missense_variant")
  x <- as.data.frame(sample(type, nrow(observed_variants), replace=TRUE))
  return(x)
}

# For each row randomly assign synon or missense
# NBNBNB Need to back calculate from number of (non)synon sites - it's not counts which should happen at equal freq but rate if selection = 0!
y <- as.data.table(replicate(10, randomise()))
replication.names <- paste("replication",1:ncol(y),sep=".")
setnames(y,replication.names)

# Append randomised annotations to main table
observed_variants <- cbind(observed_variants,y)

# Set up table with list of transcripts
random.counts <- as.data.table(unique(observed_variants$transcript.id))
setnames(crandom.counts,"transcript.id")

# Count total synon per transcript for each randomisation set
for (colj in replication.names){
	print(colj)
	X <- observed_variants[get(colj)=="synonymous_variant",list(S=sum(donors.affected)), by=list(transcript.id)]
	setnames(X,c("transcript.id",colj))
	setkey(random.counts,transcript.id)
	setkey(X,transcript.id)
	random.counts <- merge(colvar,X, all=TRUE)
}

# In future need to get patient non aggregated so each row is an observed mutation, then just sum col==synon, then will get 0's and can more easily join lists as equal length and preserve order


synon.count <- observed_variants[variant.class=="synonymous_variant",list(S=sum(donors.affected)), by=list(transcript.id)]

nonsynon.count <- observed_variants[variant.class=="missense_variant" | variant.class=="frameshift_variant" | variant.class=="disruptive_inframe_deletion" | variant.class=="disruptive_inframe_insertion" | variant.class=="inframe_deletion" | variant.class=="inframe_insertion" | variant.class=="start_lost" | variant.class=="stop_lost", list(N=sum(donors.affected)), by=list(transcript.id)]

setkey(synon.count,transcript.id)
setkey(nonsynon.count,transcript.id)

# Full outer join
counts <- merge(synon.count,nonsynon.count, all=TRUE)
setkey(counts, transcript.id)
# counts[is.na(S),]
# 3,270, now 3366
# counts[is.na(N),]
# 514, now 413

transcript.gene <- unique(observed_variants[variant.class!="exon_variant",gene.name, by=transcript.id])
setkey(transcript.gene, transcript.id)

# All rows in transcript.gene + rows that match from counts - should have equal length 81,216
named.counts <- counts[transcript.gene,]
setkey(named.counts, transcript.id)


expected_variants <- fread("data/nonsynonymous_sites_per_transcript_clean.tsv", header=TRUE)
setkey(expected_variants, transcript)

updated.annotations <- fread(input = 'zcat < data/raw/mart_export.txt.gz')
setnames(updated.annotations, make.names(names(updated.annotations)))
setkey(updated.annotations, "Ensembl.Transcript.ID")

# All rows in expected_variants + any rows from updated.annotations which match
expected_variants <- updated.annotations[expected_variants,]
# length = 72,127

# Remove non protein coding transcripts
expected_variants <- expected_variants[Transcript.type=="protein_coding" | is.na(Transcript.type),]
# 71,147, 980 removed

# All rows in named.counts with rows from expected_variants which match
expected_variants <- expected_variants[named.counts,nomatch=0]
# some transcripts with observed variations do not have calculated nonsynon sites due to N or not multiple of 3, generating NA in e.g. chromosome column as these transcripts were not passed through to N per transcript. Inner Join
# 63,747

# Calculate synonymous sites
expected_variants$synonymous_sites <- expected_variants$cds_length - expected_variants$nonsynonymous_sites

# Calculate dNdS
expected_variants$dNdS <- (expected_variants$N/expected_variants$nonsynonymous_sites)/(expected_variants$S/expected_variants$synonymous_sites)

# for 208 transcripts expected_variants[is.na(N), generating a dNdS of 0
# for 1931 transcripts S==NA, generating a dNdS of Inf
# No rows have both N and S ==0

# Remove dNdS NA rows
expected_variants <- expected_variants[is.na(dNdS)==FALSE,]
#61608 cds

write.table(expected_variants, "data/dNdS_by_transcript.tsv", sep="\t",row.names=FALSE,quote=FALSE)
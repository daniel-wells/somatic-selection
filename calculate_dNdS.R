library(data.table)
variations <- fread("data/observed_variants_by_transcript.filtered.tsv",sep="|", header=FALSE)

setnames(variations, c("gene.name", "gene.id", "strand", "transcript.name", "transcript.id", "protein.id", "variant.class", "nt.mutation", "aa.mutation", "projects","donors.affected","mutation","project.count","tested.donors"))
setkey(variations, transcript.id)

# Check number of each variant.class
variations[,.N,by=variant.class][order(-N)]

synon.count <- variations[variant.class=="synonymous_variant",list(S=sum(donors.affected)), by=list(transcript.id)]

nonsynon.count <- variations[variant.class=="missense_variant" | variant.class=="frameshift_variant" | variant.class=="disruptive_inframe_deletion" | variant.class=="disruptive_inframe_insertion" | variant.class=="inframe_deletion" | variant.class=="inframe_insertion" | variant.class=="start_lost" | variant.class=="stop_lost", list(N=sum(donors.affected)), by=list(transcript.id)]

# NorS.variations <- variations[variant.class=="missense_variant" | variant.class=="synonymous_variant",]

#ENST00000000233
# ENST00000288602 #BRAF

#synon.count <- NorS.variations[,list(S=sum(variant.class=="synonymous_variant")), by=list(transcript.id)]
#nonsynon.count <- NorS.variations[, list(N=sum(variant.class=="missense_variant")), by=list(transcript.id)]
# 81,216 unique transcripts

#nonsynon.count <- variations[variant.class=="missense_variant", list(N=.N), by=list(transcript.id)]
#synon.count <- variations[variant.class=="synonymous_variant", list(S=.N), by=list(transcript.id)]
setkey(synon.count,transcript.id)
setkey(nonsynon.count,transcript.id)

# All rows in nonsynon.count with rows in synon.count which match (as 0 counts included previously length stays same)
#counts <- synon.count[nonsynon.count,]

# Full outer join
counts <- merge(synon.count,nonsynon.count, all=TRUE)
setkey(counts, transcript.id)
# counts[is.na(S),]
# 3,270, now 3366
# counts[is.na(N),]
# 514, now 413

variations[,gene.name, by=transcript.id]

#transcript.gene <- unique(variations[variant.class=="missense_variant" | variant.class=="synonymous_variant",gene.name, by=transcript.id])

transcript.gene <- unique(variations[variant.class!="exon_variant",gene.name, by=transcript.id])
setkey(transcript.gene, transcript.id)

# All rows in transcript.gene + rows that match from counts - should have equal length 81,216
named.counts <- counts[transcript.gene,]
setkey(named.counts, transcript.id)


cds <- fread("data/N_per_transcript_clean.tsv", header=TRUE)
setkey(cds, transcript)

updated.annotations <- fread(input = 'zcat < data/raw/mart_export.txt.gz')
setnames(updated.annotations, make.names(names(updated.annotations)))
setkey(updated.annotations, "Ensembl.Transcript.ID")

# All rows in cds + any rows from updated.annotations which match
cds <- updated.annotations[cds,]
# length = 72,127

# Remove non protein coding transcripts
cds <- cds[Transcript.type=="protein_coding" | is.na(Transcript.type),]
# 71,147, 980 removed

# All rows in named.counts with rows from cds which match
cds <- cds[named.counts,nomatch=0]
# some transcripts with observed variations do not have calculated nonsynon sites due to N or not multiple of 3, generating NA in e.g. chromosome column as these transcripts were not passed through to N per transcript. Inner Join
# 63,747

# Calculate synonymous sites
cds$synonymous_sites <- cds$cds_length - cds$nonsynonymous_sites

# Calculate dNdS
cds$dNdS <- (cds$N/cds$nonsynonymous_sites)/(cds$S/cds$synonymous_sites)

# for 208 transcripts cds[is.na(N), generating a dNdS of 0
# for 1931 transcripts S==NA, generating a dNdS of Inf
# No rows have both N and S ==0

# Remove dNdS NA rows
cds <- cds[is.na(dNdS)==FALSE,]
#61608 cds

write.table(cds, "data/dNdS-2016.01.29.tsv", sep="\t",row.names=FALSE,quote=FALSE)
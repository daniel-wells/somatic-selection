# Start writing to an output file
logfile <- file(paste("logs/calculate_dNdS.R.log",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "txt", sep = "."))
sink(logfile)
sink(logfile, type="message")

source("code/functions.R")
library(data.table)

observed_variants <- readRDS("data/coding.mutations.filtered.rds")
setkey(observed_variants,icgc_mutation_id,icgc_donor_id,Ensembl.Transcript.ID)

print("observed_variants summary:")
print(summary(observed_variants))

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

observed_variants[,Ensembl.Transcript.ID:=as.character(Ensembl.Transcript.ID)]
setkey(observed_variants,icgc_mutation_id,icgc_donor_id,Ensembl.Transcript.ID)

# Check number of each variant.class
observed_variants[,.N,by=variant.class][order(-N)]

synon.count <- observed_variants[variant.class=="synonymous_variant",.("S"=.N), by=list(Ensembl.Transcript.ID,project_code)]

nonsynon.count <- observed_variants[variant.class %in% c("missense_variant","frameshift_variant","disruptive_inframe_deletion","disruptive_inframe_insertion","inframe_deletion","inframe_insertion","start_lost","stop_lost","stop_gained"), .("N"=.N), by=list(Ensembl.Transcript.ID,project_code)]

setkey(synon.count,Ensembl.Transcript.ID,project_code)
setkey(nonsynon.count,Ensembl.Transcript.ID,project_code)

# Full outer join
counts <- merge(synon.count,nonsynon.count, all=TRUE)

print(paste(nrow(counts[is.na(S),]),"project-transcripts with synon count of 0"))
print(paste(nrow(counts[is.na(N),]),"project-transcripts with nonsynon count of 0"))

# No rows with both N and S NA - i.e. all NA's are actually Zeros - gene must have been sequenced if one of the two catagories are >0
# Non sequenced genes are just not in this list as it's generated from variant observed list
print(paste(nrow(counts[is.na(N) & is.na(S),]),"transcripts with both nonsynon and synon count of 0"))

stopifnot(nrow(counts[is.na(N) & is.na(S),])==0)

# Therefore set NA's to 0
counts$N[is.na(counts$N)] <- 0
counts$S[is.na(counts$S)] <- 0

# Add project grouping information
setkey(counts, project_code)
counts <- project.groupings[counts]

print("Counts by project summary:")
print(summary(counts))

# counts by primary site
counts.bysite <- counts[,.(S=sum(S),N=sum(N)),by=.(primary_site,Ensembl.Transcript.ID)]

print("Counts by site summary:")
print(summary(counts.bysite))
print(paste(nrow(counts.bysite[S==0]),"site-transcripts with synon count of 0"))
print(paste(nrow(counts.bysite[N==0]),"site-transcripts with nonsynon count of 0"))

# counts overall
counts.pancancer <- counts[,.(S=sum(S),N=sum(N)),by=Ensembl.Transcript.ID]

print("Counts pancancer summary:")
print(summary(counts.pancancer))
print(paste(nrow(counts.pancancer[S==0]),"site-transcripts with synon count of 0"))
print(paste(nrow(counts.pancancer[N==0]),"site-transcripts with nonsynon count of 0"))

# Load expected variants per transcript
expected_variants <- fread("data/expected_variants_per_transcript.tsv", header=TRUE)
print(paste(nrow(expected_variants),"expected variants per transcript loaded"))

# Add project grouping information
setkey(expected_variants, project_code)
expected_variants <- project.groupings[expected_variants]
setkey(expected_variants,Ensembl.Transcript.ID,primary_site)

print("expected_variants summary:")
print(summary(expected_variants))

# aggregate by primary_site
expected_variants.bysite <- unique(expected_variants[,.(nonsynon.probability=sum(nonsynon.probability),synon.probability=sum(synon.probability),cds.length,chromosome,chromosome.start,chromosome.stop,Ensembl.Gene.ID,Associated.Gene.Name),by=.(Ensembl.Transcript.ID,primary_site)])

# aggregate over all cancers
expected_variants.pancancer <- unique(expected_variants[,.(nonsynon.probability=sum(nonsynon.probability),synon.probability=sum(synon.probability),cds.length,chromosome,chromosome.start,chromosome.stop,Ensembl.Gene.ID,Associated.Gene.Name),by=Ensembl.Transcript.ID])

# Join expected with actual counts
setkey(expected_variants,Ensembl.Transcript.ID,project_code)
setkey(counts,Ensembl.Transcript.ID,project_code)
# All rows in counts with rows from expected_variants which match
expected_variants.byproject <- expected_variants[counts,nomatch=0]
# some transcripts with observed variations do not have calculated nonsynon sites due to N or not multiple of 3, generating NA in e.g. chromosome column as these transcripts were not passed through to N per transcript. Inner Join
print(paste(nrow(expected_variants.byproject),"expected variants per project-transcript with measured actual variants"))

setkey(expected_variants.bysite,Ensembl.Transcript.ID,primary_site)
setkey(counts.bysite,Ensembl.Transcript.ID,primary_site)
expected_variants.bysite <- expected_variants.bysite[counts.bysite,nomatch=0]

setkey(expected_variants.pancancer,Ensembl.Transcript.ID)
setkey(counts.pancancer,Ensembl.Transcript.ID)
expected_variants.pancancer <- expected_variants.pancancer[counts.pancancer,nomatch=0]

# Calculate dNdS
expected_variants.pancancer[,expected.ratio:=nonsynon.probability/synon.probability]
expected_variants.pancancer[,observed.ratio:=N/S]
expected_variants.pancancer[,odds.ratio:=observed.ratio/expected.ratio]

expected_variants.bysite[,expected.ratio:=nonsynon.probability/synon.probability]
expected_variants.bysite[,observed.ratio:=N/S]
expected_variants.bysite[,odds.ratio:=observed.ratio/expected.ratio]

expected_variants.byproject[,expected.ratio:=nonsynon.probability/synon.probability]
expected_variants.byproject[,observed.ratio:=N/S]
expected_variants.byproject[,odds.ratio:=observed.ratio/expected.ratio]


# log ratio
expected_variants.pancancer[,log2.odds.ratio:=log(odds.ratio,2)]
expected_variants.bysite[,log2.odds.ratio:=log(odds.ratio,2)]
expected_variants.byproject[,log2.odds.ratio:=log(odds.ratio,2)]
# Some with N=0 leads to infinite log.dNdS values

print(paste(nrow(expected_variants.bysite[N==0]),"site-transcripts have N==0 generating a dNdS of 0"))
print(paste(nrow(expected_variants.bysite[S==0]),"site-transcripts have S==0 generating a dNdS of Inf"))
print(paste(nrow(expected_variants.bysite[S==0 & N==0]),"site-transcripts have S and N ==0"))
print(paste(nrow(expected_variants.bysite[is.na(expected.ratio)]),"site-transcripts have dS as NaN generating a dNdS of NaN:"))

print(paste(nrow(expected_variants.byproject[N==0]),"project-transcripts have N==0 generating a dNdS of 0"))
print(paste(nrow(expected_variants.byproject[S==0]),"project-transcripts have S==0 generating a dNdS of Inf"))
print(paste(nrow(expected_variants.byproject[S==0 & N==0]),"project-transcripts have S and N ==0"))
print(paste(nrow(expected_variants.byproject[is.na(expected.ratio)]),"project-transcripts have dS as NaN generating a dNdS of NaN:"))

print(paste(nrow(expected_variants.pancancer[N==0]),"pancancer-transcripts have N==0 generating a dNdS of 0"))
print(paste(nrow(expected_variants.pancancer[S==0]),"pancancer-transcripts have S==0 generating a dNdS of Inf"))
print(paste(nrow(expected_variants.pancancer[S==0 & N==0]),"pancancer-transcripts have S and N ==0"))
print(paste(nrow(expected_variants.pancancer[is.na(expected.ratio)]),"pancancer-transcripts have dS as NaN generating a dNdS of NaN:"))
expected_variants[synon.probability==0 | nonsynon.probability==0]


##### calculate p-values #####
expected_variants.pancancer[,"total.mut":=S+N]
expected_variants.pancancer[,"prob.s":=synon.probability/(nonsynon.probability+synon.probability)]

expected_variants.bysite[,"total.mut":=S+N]
expected_variants.bysite[,"prob.s":=synon.probability/(nonsynon.probability+synon.probability)]

expected_variants.byproject[,"total.mut":=S+N]
expected_variants.byproject[,"prob.s":=synon.probability/(nonsynon.probability+synon.probability)]

print("Calculating P-values")
test <- function(x, p, n){binom.test(x, n, p, alternative="two.sided")$p.value}
expected_variants.pancancer[,p.value:=mapply(test, S, prob.s, total.mut)]
expected_variants.bysite[,p.value:=mapply(test, S, prob.s, total.mut)]
expected_variants.byproject[,p.value:=mapply(test, S, prob.s, total.mut)]

# Adjust by BH method
expected_variants.pancancer[,q.value:=p.adjust(p.value, method = "fdr")]
expected_variants.bysite[,q.value:=p.adjust(p.value, method = "fdr"),by=primary_site]
expected_variants.byproject[,q.value:=p.adjust(p.value, method = "fdr"),by=project_code]

print("By project summary:")
print(summary(expected_variants.byproject))

print("By site summary:")
print(summary(expected_variants.bysite))

print("Pancancer summary:")
print(summary(expected_variants.pancancer))

# Check odds ratio is not 0/0
stopifnot(nrow(expected_variants.byproject[is.na(odds.ratio)==TRUE,])==0)
stopifnot(nrow(expected_variants.bysite[is.na(odds.ratio)==TRUE,])==0)
stopifnot(nrow(expected_variants.pancancer[is.na(odds.ratio)==TRUE,])==0)


archive.file("data/dNdS_byproject.tsv")
write.table(expected_variants.byproject[order(p.value)], "data/dNdS_byproject.tsv", sep="\t", row.names=FALSE, quote=FALSE)

archive.file("data/dNdS_bysite.tsv")
write.table(expected_variants.bysite[order(p.value)], "data/dNdS_bysite.tsv", sep="\t", row.names=FALSE, quote=FALSE)

archive.file("data/dNdS_pancancer.tsv")
write.table(expected_variants.pancancer[order(p.value)], "data/dNdS_pancancer.tsv", sep="\t", row.names=FALSE, quote=FALSE)



################# Plot distribution of synon.variants ##################

# remove non-numeric to get nt position of mut
observed_variants[,cds_position:=as.numeric(gsub("[^0-9]","",cds_mutation))]

plot.mutdist <- function(transcript.list){
ggplot(observed_variants[Ensembl.Transcript.ID %in% transcript.list & variant.class=="synonymous_variant",.(cds_position,Associated.Gene.Name)], aes(cds_position)) + 
	geom_histogram() + 
	facet_wrap(~Associated.Gene.Name,scales="free") + 
	labs(title="Distribution of Synonymous variants")
	# todo add colour = $is.synonymous
}

plot.mutdist.N <- function(transcript.list){
ggplot(observed_variants[Ensembl.Transcript.ID %in% transcript.list & variant.class!="synonymous_variant",.(cds_position,Associated.Gene.Name)], aes(cds_position)) + 
	geom_histogram() + 
	facet_wrap(~Associated.Gene.Name,scales="free") + 
	labs(title="Distribution of Nonsynonymous variants")
	# todo add colour = $is.synonymous
}

transcript.list <- observed_variants[order(p.value)][odds.ratio<1][1:36]$Ensembl.Transcript.ID

transcript.list.top <- observed_variants[order(p.value)][odds.ratio>1][1:36]$Ensembl.Transcript.ID


pdf("results/QC_mutation_distribution.pdf",width=16, height=9, onefile = TRUE)
plot.mutdist(transcript.list.top)
plot.mutdist(transcript.list)
plot.mutdist.N(transcript.list.top)
plot.mutdist.N(transcript.list)
dev.off()


sessionInfo()

sink(type="message")
sink()
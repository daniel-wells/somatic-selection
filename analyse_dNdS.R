source("code/functions.R")
# Start writing to an output file
logfile <- file(paste("logs/analyse_dNdS.R.log",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "txt", sep = "."))
sink(logfile)
sink(logfile, type="message")
# module load apps/R/3.2.2/gcc-4.4.7+lapack-3.5.0+blas-20110419

# Using http mirror london
# install.packages('ggplot2')
# install.packages("ggrepel")
# install.packages('data.table')
# install.packages('metRology')

library(data.table)
cds <- fread("data/dNdS_by_transcript.tsv", header=TRUE)
print(paste(nrow(cds),"transcripts dNdS values loaded."))

setnames(cds,"gene.id","gene")
setnames(cds,"Associated.Gene.Name","gene.name")
setnames(cds,"cds.length","cds_length")

# include genes with observed S==0, calculate new dNdS by setting S==0.5
cds[is.infinite(dNdS),S.2:=0.5]
cds[S.2==0.5,dNdS:=dN/(S.2/synon.probability)]

dNdS.by.gene <- cds

# set key as all columns
setkey(dNdS.by.gene,gene)

# Calculate mean dNdS per gene to remove duplicates
dNdS.by.gene <- unique(dNdS.by.gene[is.finite(dNdS),.(dNdS = mean(dNdS,na.rm=TRUE),gene.name,chromosome,chromosome.start,strand,Ensembl.Transcript.ID,cds_length,uS=mean(S),uN=mean(N),ucS=mean(synon.probability),ucN=mean(nonsynon.probability)),by=gene])

# add ranking and order
dNdS.by.gene[,ranking:=rank(dNdS,ties.method="first")]
dNdS.by.gene <- dNdS.by.gene[order(-ranking)]
print(paste(nrow(dNdS.by.gene),"genes"))

# log dNdS as it's a ratio
dNdS.by.gene$Log.dNdS <- log10(dNdS.by.gene$dNdS)
# Some with N=0 leads to infinite log.dNdS values

# Remove genes with poor average Alignability of 100mers (GEM from ENCODE/CRG(Guigo))
blacklist <- fread('zcat < data/raw/exons.hg19.mappability100.bed.gz')
setnames(blacklist,c("chromosome","start","end","exon","V5","strand","gene.id","gene.name","mappability"))

blacklist <- blacklist[,.("mappability"=mean(mappability)),by=gene.id]

# hist(blacklist$mappability,breaks=1000,ylim=c(0,500))

# for each row in blacklist, add data from dNdS.by.gene
setkey(dNdS.by.gene,gene)
dNdS.by.gene <- dNdS.by.gene[blacklist][!is.na(dNdS)]
dNdS.by.gene <- dNdS.by.gene[mappability>0.75]
paste(nrow(blacklist[mappability<0.75]),"genes removed due to mappability <0.75")
paste(nrow(dNdS.by.gene),"genes remaining")

######################
###### Annotate ######
######################

##### annotate COSMIC cancer genes #####
cancer.genes <- fread("data/raw/cancer_gene_census.csv", header=TRUE)
setnames(cancer.genes, make.names(names(cancer.genes)))


# Correct an filter COSMIC genes via HGNC
hgnc <- fread("data/raw/HGNC.tsv", header=TRUE)
setnames(hgnc,make.names(names(hgnc)))
setkey(hgnc,Entrez.Gene.ID)
setkey(cancer.genes,Entrez.GeneId)
cancer.genes <- hgnc[cancer.genes]

print("Genes in Cancer Census with no Ensembl.ID")
hgnc[cancer.genes][is.na(Approved.Symbol),.(Approved.Symbol,Gene.Symbol,Locus.Group,Entrez.Gene.ID,Ensembl.ID)]

# Remove non protein coding genes from cancer list
cancer.genes <- cancer.genes[Locus.Group=="protein-coding gene"]

# Remove non somaticly mutated genes from cancer list
cancer.genes <- cancer.genes[Somatic=="yes"]

dNdS.by.gene$cancer.gene <- dNdS.by.gene$gene %in% cancer.genes$Ensembl.ID
# NB some duplicates, 518 unique

print(paste(nrow(dNdS.by.gene[cancer.gene==TRUE]),"cancer genes with dNdS value"))

print("Cancer genes not in dNdS.by.gene database (S or N = 0, or mappability bad)")
cancer.genes[!(cancer.genes$Ensembl.ID %in% dNdS.by.gene$gene),.(Gene.Symbol,Entrez.Gene.ID,Ensembl.ID)]

##### annotate vogelstein cancer genes ######
cancer.genes.vogelstein <- fread("data/raw/vogelstein_driver_genes.tdv", header=TRUE)

# Correct Names
cancer.genes.vogelstein[gene_name=="FAM123B",gene_name:="AMER1"]
cancer.genes.vogelstein[gene_name=="MLL2",gene_name:="KMT2D"]
cancer.genes.vogelstein[gene_name=="MLL3",gene_name:="KMT2C"]

dNdS.by.gene$cancer.gene.vogelstein <- dNdS.by.gene$gene.name %in% cancer.genes.vogelstein$gene_name

print(paste(nrow(dNdS.by.gene[cancer.gene.vogelstein==TRUE]),"vogelstein cancer genes with dNdS value"))

# Which genes in vogelstein do not map to dNdS.by.gene database
print("Genes in vogelstein list not matched to dNdS value:")
cancer.genes.vogelstein[!(cancer.genes.vogelstein$gene_name %in% dNdS.by.gene$gene.name),gene_name]

# For all rows in dNdS.by.gene, add classification from cancer.genes.vogelstein if possible
setkey(cancer.genes.vogelstein,gene_name)
setkey(dNdS.by.gene,gene.name)
dNdS.by.gene$classification <- cancer.genes.vogelstein[dNdS.by.gene,classification]


##### add expression data ######
RNAseq <- fread("data/RNAseq.by.gene.tsv", header=TRUE)
setnames(RNAseq,c("gene.name","mean.expression","mean.of.stdev.of.expression","expression.ranking","expression.percent.rank","log10.expression"))
setkey(RNAseq,gene.name)

# Which dNdS.by.gene gene names not in RNAseq database
print(paste(length(dNdS.by.gene[!(dNdS.by.gene$gene.name %in% RNAseq$gene.name),.(gene.name)]$gene.name),"genes not in RNAseq database"))

# For all rows in dNdS.by.gene, add expression from RNAseq if avaliable
dNdS.by.gene <- RNAseq[dNdS.by.gene,]


# to detect bias calculate dN, dS, and ranking given a cutoff
dNdS.by.gene$dS <- dNdS.by.gene$uS / dNdS.by.gene$ucS
dNdS.by.gene$dN <- dNdS.by.gene$uN / dNdS.by.gene$ucN

# Re-rank after subsetting uS
dNdS.by.gene[uS>3,ranking.1:=rank(dNdS,ties.method="first")]
dNdS.by.gene[uS>10,ranking.2:=rank(dNdS,ties.method="first")]
dNdS.by.gene[uS>15,ranking.3:=rank(dNdS,ties.method="first")]
dNdS.by.gene[uS>3 & uN>3,ranking.4:=rank(dNdS,ties.method="first")]


#########################
###### Plot Graphs ######
#########################

library(ggplot2)
library(ggrepel)
archive.file("results/results.pdf")
pdf("results/results.pdf", width=16, height=9, onefile = TRUE)

# S Histogram (Background mutation rate power)
ggplot(dNdS.by.gene, aes(uS)) + geom_histogram(binwidth = 1) + scale_x_continuous(limits = c(0, 100)) + labs(x="mean synonymous mutations per gene",title="Distribution of S per gene")

# N Histogram (Background mutation rate power)
ggplot(dNdS.by.gene, aes(uN)) + geom_histogram(binwidth = 1) + scale_x_continuous(limits = c(0, 100)) + labs(x="mean nonsynonymous mutations per gene",title="Distribution of N per gene")

# Graphs of dS by ranking
ggplot(dNdS.by.gene, aes(ranking, dS)) + geom_point(alpha=0.3) + geom_smooth() + ylim(0,750) + labs(title="dS by ranking for all genes") + annotate("text", x = 1000, y=0, label = paste(nrow(dNdS.by.gene[is.na(ranking)==FALSE]),"genes"))

ggplot(dNdS.by.gene, aes(ranking.1, dS)) + geom_point(alpha=0.3) + geom_smooth() + ylim(0,750) + labs(title="dS by ranking for uS>3") + annotate("text", x = 1000, y=0, label = paste(nrow(dNdS.by.gene[is.na(ranking.1)==FALSE]),"genes"))

ggplot(dNdS.by.gene, aes(ranking.4, dS)) + geom_point(alpha=0.3) + geom_smooth() + ylim(0,750) + labs(title="dS by ranking for uS>3 & uN>3") + annotate("text", x = 1000, y=0, label = paste(nrow(dNdS.by.gene[is.na(ranking.4)==FALSE]),"genes"))

ggplot(dNdS.by.gene, aes(ranking.2, dS)) + geom_point(alpha=0.3) + geom_smooth() + ylim(0,750) + labs(title="dS by ranking for uS>10") + annotate("text", x = 1000, y=0, label = paste(nrow(dNdS.by.gene[is.na(ranking.2)==FALSE]),"genes"))

ggplot(dNdS.by.gene, aes(ranking.3, dS)) + geom_point(alpha=0.3) + geom_smooth() + ylim(0,750) + labs(title="dS by ranking for uS>15") + annotate("text", x = 1000, y=0, label = paste(nrow(dNdS.by.gene[is.na(ranking.3)==FALSE]),"genes"))


# Funnel plot (uS by dNdS)
# with uS as 'study size' (size approximating power which is background mutation rate, more than just cds_length)
up.conf <- function(x) {(5/(x+3))+0.1}
up.conf.l <- function(x) {(5/(x+3))+0.12}

low.conf <- function(x) {(-7/(x+2))-0.3}
low.conf.l <- function(x) {(-7/(x+2))-0.3}

ggplot(dNdS.by.gene[is.finite(Log.dNdS)], aes(uS,Log.dNdS)) + geom_point(aes(colour = cancer.gene),alpha=0.3) + xlim(0,100) + ylim(-1.6,+2.1) + scale_colour_manual(name="In COSMIC cancer gene census?",  values =c("black", "red")) + theme(legend.position="bottom") + stat_function(fun = up.conf)+ stat_function(fun = low.conf) + geom_label_repel(data = dNdS.by.gene[is.finite(Log.dNdS)][Log.dNdS<low.conf.l(uS) | Log.dNdS>up.conf.l(uS)], aes(label = gene.name), size = 2, box.padding = unit(0.5, "lines"), point.padding = unit(0.1, "lines"), force=1,segment.size=0.2,segment.color="blue")

print(paste("Total Vogelstein cancer genes:",sum(dNdS.by.gene$cancer.gene.vogelstein)))
print(paste("Vogelstein cancer genes in top reigon:",sum(dNdS.by.gene[Log.dNdS>up.conf(uS)]$cancer.gene.vogelstein)))
print(paste("All genes in top reigon:",nrow(dNdS.by.gene[Log.dNdS>up.conf(uS)])))
print(paste("Total cancer genes:",sum(dNdS.by.gene$cancer.gene)))
print(paste("Cancer genes in top reigon:",sum(dNdS.by.gene[Log.dNdS>up.conf(uS)]$cancer.gene)))


# Order chromosomes
dNdS.by.gene$chromosome <- factor(dNdS.by.gene$chromosome, levels = c(1:22,"X","Y","MT"))

# Graphs of dNdS by chromosome position
ggplot(dNdS.by.gene, aes(chromosome.start, Log.dNdS)) + geom_point(alpha=0.3,size=0.05) + facet_wrap(~chromosome,scales="free_x")

# Overall Dist
ggplot(dNdS.by.gene[uS>3], aes(dNdS)) + geom_histogram(binwidth = 0.01) + geom_vline(xintercept = 1,color = "red") + theme_grey(base_size = 30) + labs(x="mean dN/dS per gene",title="Overall dN/dS Distribution") + scale_x_log10()

# Fit data to students-t and normal dist
library(MASS)
paramst <- as.list(MASS::fitdistr(dNdS.by.gene[is.finite(Log.dNdS),Log.dNdS], "t")$estimate)
paramsn <- as.list(MASS::fitdistr(dNdS.by.gene[is.finite(Log.dNdS),Log.dNdS], "normal")$estimate)

# QQ plot t dist
ggplot(dNdS.by.gene, aes(sample = Log.dNdS)) + stat_qq(distribution = qt, dparams = paramst["df"])+ labs(title="QQ plot, students-t distribution")

# QQ plot normal dist
ggplot(dNdS.by.gene, aes(sample = Log.dNdS)) + stat_qq(distribution = qnorm, dparams = paramsn)+ labs(title="QQ plot, normal distribution")

library(metRology)
# Fitting t dist over actual, args from paramst
ggplot(dNdS.by.gene, aes(Log.dNdS)) + geom_histogram(aes(y=..density..),binwidth=0.01) + geom_vline(xintercept = 0,color = "blue") + theme_grey(base_size = 30) + labs(x="Log10(mean dN/dS per gene)",title="Overall Distribution")  + stat_function(geom = "line", fun = dt.scaled, arg = list(df = paramst$df, mean = paramst$m, sd = paramst$s), colour = "red") + annotate("text", x = 1, y=2.5, label = paste("df:",signif(paramst$df,3),"\n Mean: ",signif(paramst$m,3),"\n sd:",signif(paramst$s,3)))

print(paste("T-distribution parameters (of Log.dNdS): df:",signif(paramst$df,3),"\n Mean: ",signif(paramst$m,3)))

# Cancer vs normal density dist
ggplot(dNdS.by.gene[uS>3], aes(x=dNdS,fill=cancer.gene.vogelstein)) + geom_density(alpha=0.3) + geom_vline(xintercept = 1,color = "red") + theme_grey(base_size = 30) + labs(x="mean dN/dS per gene",title="Distribution of known cancer genes") + scale_x_log10() + theme(legend.position="bottom")

# Cancer gene dist TSG vs Onco vs normal
ggplot(dNdS.by.gene[uS>3], aes(x=dNdS,fill=classification)) + geom_density(alpha=0.3) + geom_vline(xintercept = 1,color = "red") + theme_grey(base_size = 30) + labs(x="mean dN/dS per gene",title="Distribution of known cancer genes") + scale_x_log10() + theme(legend.position="bottom")

# Ranked points graph (sideways S)
ggplot(dNdS.by.gene[uS>3], aes(x=ranking,y=dNdS)) + geom_point(aes(colour = cancer.gene),alpha=0.3,shape=21,size=0.5)  + geom_hline(yintercept = 1,color = "red") + scale_y_log10() + xlim(-700,20000)+ theme_grey(base_size = 10) +  geom_text_repel(data = subset(dNdS.by.gene[uS>3], dNdS>4.25 | dNdS<0.12), aes(label = gene.name), size = 2, box.padding = unit(0.5, "lines"), point.padding = unit(0.1, "lines"), force=1,segment.size=0.2) + scale_colour_manual(name="In COSMIC cancer gene census?",  values =c("black", "red")) + labs(y="mean dN/dS per gene",title="Overall Distribution") + theme(legend.position="bottom")  

dev.off()

bottom <- dNdS.by.gene[is.finite(Log.dNdS) & Log.dNdS<low.conf(uS),.(cancer.gene,uS,uN,dNdS,gene.name,ranking,expression.percent.rank)][order(ranking)]
# Bottom Hits



top <- dNdS.by.gene[Log.dNdS>up.conf(uS),.(cancer.gene,cancer.gene.vogelstein,uS,uN,dNdS,gene.name,ranking,expression.percent.rank)][order(-ranking)]
# Top Hits

# Save whole table
archive.file("results/dNdS.tsv")
write.table(dNdS.by.gene, "results/dNdS.tsv", sep="\t",row.names=FALSE,quote=FALSE)

# Export top and bottom 25 hits to tsv file
archive.file("results/top_dNdS.tsv")
write.table(top, "results/top_dNdS.tsv", sep="\t",row.names=FALSE,quote=FALSE)
archive.file("results/bottom_dNdS.tsv")
write.table(bottom, "results/bottom_dNdS.tsv", sep="\t",row.names=FALSE,quote=FALSE)

sessionInfo()

sink(type="message")
sink()
# Start writing to an output file
logfile <- file(paste("results/analyse_dNdS.log",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "txt", sep = "."))
sink(logfile)
sink(logfile, type="message")
# module load apps/R/3.2.2/gcc-4.4.7+lapack-3.5.0+blas-20110419

# Using http mirror london
# install.packages('ggplot2')
# install.packages("ggrepel")
# install.packages('data.table')
# install.packages('metRology')

library(data.table)
cds <- fread("data/dNdS_by_transcript.2016-02-10.12-27-58.tsv", header=TRUE)
print(paste(nrow(cds),"transcripts dNdS values loaded.",))

setnames(cds,"gene.id","gene")
setnames(cds,"cds.length","cds_length")

# For each gene, find max cds length
max.cds.by.gene <- unique(cds[is.finite(dNdS),.(max.cds = max(cds_length,na.rm=TRUE)),by=gene])

# For all rows in gene.maxlength, add cds values where cds_length=max(cds_length) by gene, for each gene (multiples as sometimes multiple transcript have max cds length)
setkey(max.cds.by.gene,gene,max.cds)
setkey(cds,gene,cds_length)
dNdS.by.gene <- cds[max.cds.by.gene,]

paste(nrow(dNdS.by.gene[duplicated(dNdS.by.gene)==TRUE,]),"duplicate entries (same gene and cds length) after getting dNdS from canonical (longest) transcript")

print(paste(nrow(dNdS.by.gene[,.(stddev = sd(dNdS)),by=gene][stddev>0.1]),"genes with high standard deviation between transcripts and dNdS value"))
# treat as seperate genes???

# set key as all columns
setkey(dNdS.by.gene)

# Calculate mean dNdS per gene to remove duplicates
dNdS.by.gene <- unique(dNdS.by.gene[is.finite(dNdS),.(dNdS = mean(dNdS,na.rm=TRUE),gene.name,chromosome,chromosome.start,strand,cds_length,uS=mean(S),uN=mean(N),ucS=mean(synon.probability),ucN=mean(nonsynon.probability)),by=gene])

# add ranking and order
dNdS.by.gene[,ranking:=rank(dNdS,ties.method="first")]
dNdS.by.gene <- dNdS.by.gene[order(-ranking)]
print(paste(nrow(dNdS.by.gene),"genes"))

# log dNdS as it's a ratio
dNdS.by.gene$Log.dNdS <- log10(dNdS.by.gene$dNdS)
# Some with N=0 leads to infinite log.dNdS values

######################
###### Annotate ######
######################

cancer.genes <- fread("data/raw/cancer_gene_census.csv", header=TRUE)
setnames(cancer.genes, make.names(names(cancer.genes)))
dNdS.by.gene$cancer.gene <- dNdS.by.gene$gene.name %in% cancer.genes$Gene.Symbol

# Which cancer gene names not in dNdS.by.gene database
cancer.genes[!(cancer.genes$Gene.Symbol %in% dNdS.by.gene$gene.name),.(Gene.Symbol,Entrez.GeneId)]

# 551 / 572
# Some because not in ensembl database at all, some because either S or N = 0

cancer.genes.strict <- fread("data/raw/cosmic_cancer_curated.tsv", header=FALSE)
dNdS.by.gene$cancer.gene.strict <- dNdS.by.gene$gene.name %in% cancer.genes.strict$V1

# All cancer gene stict are in
cancer.genes.strict[!(cancer.genes.strict$V1 %in% dNdS.by.gene$gene.name),V1]
sum(dNdS.by.gene$cancer.gene.strict)
# 178 / 178

cancer.genes.vogelstein <- fread("data/raw/vogelstein_driver_genes.tdv", header=TRUE)
dNdS.by.gene$cancer.gene.vogelstein <- dNdS.by.gene$gene.name %in% cancer.genes.vogelstein$gene_name

# Which genes in vogelstein do not map to dNdS.by.gene database
print("Genes in vogelstein list not matched to dNdS value:")
cancer.genes.vogelstein[!(cancer.genes.vogelstein$gene_name %in% dNdS.by.gene$gene.name),gene_name]
# 122 / 125

# For all rows in dNdS.by.gene, add classification from cancer.genes.vogelstein if possible
setkey(cancer.genes.vogelstein,gene_name)
setkey(dNdS.by.gene,gene.name)
dNdS.by.gene$classification <- cancer.genes.vogelstein[dNdS.by.gene,classification]


# Add expression data
RNAseq <- fread("results/RNAseq.by.gene.2016-01-26.14-26-42.tsv", header=TRUE)
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

pdf(paste("results/results",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "pdf", sep = "."), width=16, height=9, onefile = TRUE)

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
# with uS as 'study size'
ggplot(dNdS.by.gene, aes(uS,Log.dNdS)) + geom_point(alpha=0.3) + xlim(0,200)

# with cds_length as 'study size'
ggplot(dNdS.by.gene, aes(cds_length,Log.dNdS)) + geom_point(alpha=0.3) + xlim(0,15000)

# Identify outlier genes
dNdS.by.gene[uS>100 & Log.dNdS<(-0.5)]
# MMLT3, TBP, DSPP
dNdS.by.gene[uS>25 & Log.dNdS<(-1.5)]
# THEM6

# Order chromosomes
dNdS.by.gene$chromosome <- factor(dNdS.by.gene$chromosome, levels = c(1:22,"X","Y","MT"))

# Graphs of dNdS by chromosome position
ggplot(dNdS.by.gene, aes(chromosome.start, Log.dNdS)) + geom_point(alpha=0.3,size=0.05) + facet_wrap(~chromosome2,scales="free_x")

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


# Bottom 25
bottom <- dNdS.by.gene[ranking<100 & uS>3,.(cancer.gene,uS,uN,dNdS,gene.name,ranking,expression.percent.rank)][order(ranking)]

# Get all USP17L genes with uS>3
USP17L <- dNdS.by.gene[grep("USP17",dNdS.by.gene$gene.name),.(gene.name,chromosome,uS,uN,ranking,dNdS,Log.dNdS,expression.percent.rank)][uS>3][order(ranking)]

dNdS.by.gene$USP17L <- dNdS.by.gene$gene.name %in% USP17L$gene.name

# USP17L vs normal density dist + theme_grey(base_size = 30)
ggplot(dNdS.by.gene, aes(x=dNdS,fill=USP17L)) + geom_density(alpha=0.3) + geom_vline(xintercept = 1,color = "red")  + labs(x="mean dN/dS per gene",title="Distribution of USP17L genes vs normal") + scale_x_log10() + theme(legend.position="bottom")

dev.off()

# Top 75
top <- dNdS.by.gene[ranking>max(ranking)-400 & uS>3,.(cancer.gene,uS,uN,dNdS,gene.name,ranking,expression.percent.rank)][order(-ranking)]

# Save whole table
write.table(dNdS.by.gene, paste("results/dNdS",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "tsv", sep = "."), sep="\t",row.names=FALSE,quote=FALSE)

# Export top and bottom 25 hits to tsv file
write.table(top, paste("results/top_dNdS",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "tsv", sep = "."), sep="\t",row.names=FALSE,quote=FALSE)
write.table(bottom, paste("results/bottom_dNdS",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "tsv", sep = "."), sep="\t",row.names=FALSE,quote=FALSE)
write.table(USP17L, paste("results/USP17L",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "tsv", sep = "."), sep="\t",row.names=FALSE,quote=FALSE)

sessionInfo()

sink(type="message")
sink()
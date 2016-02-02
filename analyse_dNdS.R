# module load apps/R/3.2.2/gcc-4.4.7+lapack-3.5.0+blas-20110419

# Using http mirror london
# install.packages('ggplot2')
# install.packages("ggrepel")
# install.packages('data.table')
# install.packages('metRology')

library(data.table)
cds <- fread("data/dNdS-2016.01.29.tsv", header=TRUE)
# 61,553

# For each gene, find max cds length
max.cds.by.gene <- unique(cds[is.finite(dNdS),.(max.cds = max(cds_length,na.rm=TRUE)),by=gene])

# For all rows in gene.maxlength, add cds values where cds_length=max(cds_length) by gene, for each gene (multiples as sometimes multiple transcript have max cds length)
setkey(max.cds.by.gene,gene,max.cds)
setkey(cds,gene,cds_length)
dNdS.by.gene <- cds[max.cds.by.gene,]

# dNdS.by.gene[duplicated(dNdS.by.gene)==TRUE,]
# 5274 duplicates+ (by gene,cds_length)

# nrow(dNdS.by.gene[,.(stddev = sd(dNdS)),by=gene][stddev>0.1])
# 22 genes with high stddev per gene, treat as seperate genes???

# set key as all columns
setkey(dNdS.by.gene)

# Calculate mean dNdS per gene to remove duplicates
# dNdS.by.gene[uS<3.0000001,]
# 1618 genes with uS<3
# 17,448 genes left

# Still 56 genes with dNdS = 0 due to no N found
dNdS.by.gene[,ranking:=rank(dNdS,ties.method="first")]
dNdS.by.gene <- dNdS.by.gene[order(-ranking)]
#19,066 genes

cancer.genes <- fread("data/raw/cancer_gene_census.csv", header=TRUE)
setnames(cancer.genes, make.names(names(cancer.genes)))
dNdS.by.gene$cancer.gene <- dNdS.by.gene$gene.name %in% cancer.genes$Gene.Symbol

# Which cancer gene names not in dNdS.by.gene database
cancer.genes[!(cancer.genes$Gene.Symbol %in% dNdS.by.gene$gene.name),.(Gene.Symbol,Entrez.GeneId)]

#  Gene.Symbol Entrez.GeneId
#  1:        BCL5           603
#  2:     C12orf9         93669
#  3:    C15orf65        145788
#  4: CDKN2A(p14)          1029
#  5:      DUX4L1         22947
#  6:    HMGN2P46        283651
#  7:         IGH          3492
#  8:         IGK         50802
#  9:         IGL          3535
# 10:      MALAT1        378938
# 11:  RNF217-AS1          7955
# 12:       RPL22          6146
# 13:     RUNDC2A         84127
# 14:        SSX2          6757
# 15:        TCL6         27004
# 16:         TRA          6955
# 17:         TRB          6957
# 18:         TRD          6964
# 19:      ZNF198          7750
# 20:      ZNF278         23598

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
cancer.genes.vogelstein[!(cancer.genes.vogelstein$gene_name %in% dNdS.by.gene$gene.name),gene_name]
# [1] "FAM123B" "MLL2"    "MLL3"
# 122 / 125

# For all rows in dNdS.by.gene, add classification from cancer.genes.vogelstein if possible
setkey(cancer.genes.vogelstein,gene_name)
setkey(dNdS.by.gene,gene.name)
dNdS.by.gene$classification <- cancer.genes.vogelstein[dNdS.by.gene,classification]

dNdS.by.gene$Log.dNdS <- log10(dNdS.by.gene$dNdS)


# Add expression data

RNAseq <- fread("results/RNAseq.by.gene.2016-01-26.14-26-42.tsv", header=TRUE)
setnames(RNAseq,c("gene.name","mean.expression","mean.of.stdev.of.expression","expression.ranking","expression.percent.rank","log10.expression"))
setkey(RNAseq,gene.name)

# Which dNdS.by.gene gene names not in RNAseq database
# dNdS.by.gene[!(dNdS.by.gene$gene.name %in% RNAseq$gene.name),.(gene.name)]$gene.name
# 1,702 !

# For all rows in dNdS.by.gene, add expression from RNAseq if avaliable
dNdS.by.gene <- RNAseq[dNdS.by.gene,]

# Calculate sliding window mean of S over ranking
slideFunct <- function(data, window, step){
  total <- length(data)
  window = total / 200
  spots <- seq(from = 1, to = (total - window + 1), by = step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- mean(data[spots[i]:(spots[i] + window - 1)])
  }
  return(result)
}

sliding.mean.uS <- slideFunct(dNdS.by.gene$uS,100,100)
sliding.mean.uN <- slideFunct(dNdS.by.gene$uN,100,100)

sliding.mean.cds <- slideFunct(dNdS.by.gene$cds_length,100,100)

dNdS.by.gene$dN <- dNdS.by.gene$uN / dNdS.by.gene$cds

sliding.mean.dS <- slideFunct(dNdS.by.gene$dS,100,100)
sliding.mean.dS20 <- slideFunct(dNdS.by.gene[uS>20]$dS,37,37)
sliding.mean.dN <- slideFunct(dNdS.by.gene$dN,100,100)

pdf(width=16, height=9, onefile = TRUE)
plot(sliding.mean.uS,xlim=c(0,190),ylim=c(0,85),type="o")
lines(sliding.mean.uN,type="o",col="blue")
plot(sliding.mean.dS,type="o")
lines(sliding.mean.dS20,type="o",col="red")
plot(sliding.mean.dN,type="o")
plot(sliding.mean.cds,type="o")
dev.off()

# Plot Graphs!
library(ggplot2)
library(ggrepel)
pdf(paste("results/results",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "pdf", sep = "."), width=16, height=9, onefile = TRUE)

# S Histogram (Background mutation rate power)
ggplot(dNdS.by.gene, aes(uS)) + geom_histogram(binwidth = 1) + scale_x_continuous(limits = c(0, 100)) + labs(x="mean synonymous mutations per gene",title="Distribution of S per gene")

ggplot(dNdS.by.gene, aes(uN)) + geom_histogram(binwidth = 1) + scale_x_continuous(limits = c(0, 100)) + labs(x="mean nonsynonymous mutations per gene",title="Distribution of N per gene")

# Overall Dist
ggplot(dNdS.by.gene[uS>3], aes(dNdS)) + geom_histogram(binwidth = 0.01) + geom_vline(xintercept = 1,color = "red") + theme_grey(base_size = 30) + labs(x="mean dN/dS per gene",title="Overall dN/dS Distribution") + scale_x_log10()

# Fit data to students-t and normal dist
library(MASS)
paramst <- as.list(MASS::fitdistr(dNdS.by.gene$Log.dNdS, "t")$estimate)
paramsn <- as.list(MASS::fitdistr(dNdS.by.gene$Log.dNdS, "normal")$estimate)

# QQ plot t dist
ggplot(dNdS.by.gene, aes(sample = Log.dNdS)) + stat_qq(distribution = qt, dparams = paramst["df"])+ labs(title="QQ plot, students-t distribution")

# QQ plot normal dist
ggplot(dNdS.by.gene, aes(sample = Log.dNdS)) + stat_qq(distribution = qnorm, dparams = paramsn)+ labs(title="QQ plot, normal distribution")

library(metRology)
# Fitting t dist over actual, args from paramst
ggplot(dNdS.by.gene, aes(Log.dNdS)) + geom_histogram(aes(y=..density..),binwidth=0.01) + geom_vline(xintercept = 0,color = "blue") + theme_grey(base_size = 30) + labs(x="mean dN/dS per gene",title="Overall Distribution")  + stat_function(geom = "line", fun = dt.scaled, arg = list(df = paramst$df, mean = paramst$m, sd = paramst$s), colour = "red")
#4.838563, mean = -0.1599829, sd = 0.1631476

# Cancer vs normal density dist
ggplot(dNdS.by.gene[uS>3], aes(x=dNdS,fill=cancer.gene.vogelstein)) + geom_density(alpha=0.3) + geom_vline(xintercept = 1,color = "red") + theme_grey(base_size = 30) + labs(x="mean dN/dS per gene",title="Distribution of known cancer genes") + scale_x_log10() + theme(legend.position="bottom")

# Cancer gene dist TSG vs Onco vs normal
ggplot(dNdS.by.gene[uS>3], aes(x=dNdS,fill=classification)) + geom_density(alpha=0.3) + geom_vline(xintercept = 1,color = "red") + theme_grey(base_size = 30) + labs(x="mean dN/dS per gene",title="Distribution of known cancer genes") + scale_x_log10() + theme(legend.position="bottom")

# Ranked points graph (sideways S)
ggplot(dNdS.by.gene[uS>3], aes(x=ranking,y=dNdS)) + geom_point(aes(colour = cancer.gene),alpha=0.3,shape=21,size=0.5)  + geom_hline(yintercept = 1,color = "red") + scale_y_log10() + xlim(-700,20000)+ theme_grey(base_size = 10) +  geom_text_repel(data = subset(dNdS.by.gene[uS>3], dNdS>3.25 | dNdS<0.085), aes(label = gene.name), size = 2, box.padding = unit(0.5, "lines"), point.padding = unit(0.1, "lines"), force=1,segment.size=0.2) + scale_colour_manual(name="In COSMIC cancer gene census?",  values =c("black", "red")) + labs(y="mean dN/dS per gene",title="Overall Distribution") + theme(legend.position="bottom")  


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

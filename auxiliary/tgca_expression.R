# module load apps/R/3.2.2/gcc-4.4.7+lapack-3.5.0+blas-20110419
library(data.table)

RNAseq.files <- list.files(path="/users/bschuster/sharedscratch/Splicing/GeneExpRaw",pattern="*_hiseq_rnaseqv2_gene_exp-raw.tdv.gz",full.names=TRUE)
# Not including "ga_rnaseq" or "hiseq_totalrnaseqv2"
# 31 in total

# alternative? concatenate all files first (but doesn't output cancer type) pan.mutations <- rbindlist(lapply(paste('zcat < ',RNAseq.files), fread))

expression.by.gene.total = NULL

for (file in RNAseq.files ){
cancer.type <- sub("_hiseq_rnaseqv2_gene_exp-raw.tdv.gz","",sub("/users/bschuster/sharedscratch/Splicing/GeneExpRaw/", "", file))

RNAseq <- fread(input = paste('zcat < ',file))
print(paste("Read to data.table: ",cancer.type))

setnames(RNAseq, c("sample.id", "sample.type", "gene", "expression.score", "expression.level", "USCS.ID"))

setkey(RNAseq,gene)

# Aggregate by gene, mean expression
expression.by.gene <- RNAseq[sample.type=="TP",.(mean=mean(expression.level),stdev=sd(expression.level)),by=gene]
#20501 genes

expression.by.gene$cancer.type = cancer.type

# Append new data to totals table
expression.by.gene.total <- rbindlist(list(expression.by.gene.total,expression.by.gene),use.names=TRUE,fill=TRUE)
print(paste("Updated totals with: ",cancer.type))
}

# What types are in the second column?
# prad[,.N,by=sample.type]
#    sample.type       N
# 1:          TP 9963972
# 2:          NT 1066104
# 3:          TM   20502

# TP - Tumor Primary Solid
# NT - Normal Tissue Solid
# TM - Metastatic
# https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm?codeTable=sample%20type

write.table(expression.by.gene.total, paste("/users/dwells/sharedscratch/results/tgca_RNAseq.total",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "tsv", sep = "."), sep="\t",row.names=FALSE,quote=FALSE)

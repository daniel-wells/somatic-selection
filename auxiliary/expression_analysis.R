library(data.table)
RNAseq <- fread(input = '/users/dwells/sharedscratch/results/tgca_RNAseq.total.2016-01-26.13-39-57.tsv')

expression.by.gene <- RNAseq[,.(mean=mean(mean),mean.of.stdev=sd(stdev)),by=gene]

# Rank by expression, then examine top hits
expression.by.gene[,ranking:=rank(mean,ties.method="first")]
expression.by.gene$percent.rank <- expression.by.gene$ranking/20501
top <- expression.by.gene[ranking>max(ranking)-100,.(gene,mean,mean.of.stdev,ranking,percent.rank=(ranking/20501))][order(-ranking)]

#bottom hits
expression.by.gene[mean==0,gene]
#220 genes

# BRCA
# MGA - matrix Gla protein
# SCGB2A2 - mammaglobin-A, breast tissue marker
# TPT1 - tumor protein, translationally-controlled 1
# ACTB - Beta Actin
# GAPDH -

# PRAD
# KLK3 - PSA precursor?
# MSMB - Beta-microseminoprotein, prostate epithelia tissue expression
# NPY - neuropeptide Y receptors
# B2M - Beta-2-microglobulin

write.table(expression.by.gene, paste("/users/dwells/sharedscratch/results/RNAseq.by.gene",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "tsv", sep = "."), sep="\t",row.names=FALSE,quote=FALSE)
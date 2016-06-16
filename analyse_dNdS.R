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

library(data.table)
dNdS.by.gene.bysite <- fread("data/dNdS_bysite.tsv", header=TRUE)
dNdS.by.gene <- fread("data/dNdS_pancancer.tsv", header=TRUE)

# set key
setkey(dNdS.by.gene,Ensembl.Gene.ID)
setkey(dNdS.by.gene.bysite,Ensembl.Gene.ID)


# add ranking and order
dNdS.by.gene[,ranking:=rank(odds.ratio,ties.method="first")]

dNdS.by.gene[,or.rank:=rank(odds.ratio,ties.method="first")]
dNdS.by.gene[,p.rank:=rank(p.value,ties.method="first")]

dNdS.by.gene <- dNdS.by.gene[order(-ranking)]
print(paste(nrow(dNdS.by.gene),"genes with dNdS loaded"))

######################
###### Annotate ######
######################

dNdS.by.gene$is.significant <- dNdS.by.gene$q.value<0.1
dNdS.by.gene$is.interesting <- dNdS.by.gene$q.value>0.1 & dNdS.by.gene$p.value<0.01


##### annotate COSMIC cancer genes #####
cancer.genes <- fread("data/raw/cancer_gene_census.csv", header=TRUE)
setnames(cancer.genes, make.names(names(cancer.genes)))


# Correct an filter COSMIC genes via HGNC
hgnc <- fread("data/raw/HGNC.tsv", header=TRUE)
setnames(hgnc,make.names(names(hgnc)))
setnames(hgnc,c("Entrez.Gene.ID","Entrez.Gene.ID.supplied.by.NCBI."),c("Entrez.Gene.ID.manual","Entrez.Gene.ID"))
setnames(hgnc,c("Ensembl.Gene.ID","Ensembl.ID.supplied.by.Ensembl."),c("Ensembl.Gene.ID.manual","Ensembl.Gene.ID"))
setkey(hgnc,Entrez.Gene.ID)
setkey(cancer.genes,Entrez.GeneId)
cancer.genes <- hgnc[cancer.genes]

print("Genes in Cancer Census with no Approved.Symbol")
cancer.genes[is.na(Approved.Symbol),.(Approved.Symbol,Gene.Symbol,Locus.Group,Entrez.Gene.ID,Ensembl.Gene.ID)]

# Remove non protein coding genes from cancer list
cancer.genes <- cancer.genes[Locus.Group=="protein-coding gene"]

# Remove non somaticly mutated genes from cancer list
cancer.genes <- cancer.genes[Somatic=="yes"]

dNdS.by.gene$cancer.gene <- dNdS.by.gene$Ensembl.Gene.ID %in% cancer.genes$Ensembl.Gene.ID
# NB some duplicates, 518 unique

print(paste(nrow(dNdS.by.gene[cancer.gene==TRUE]),"cancer genes with dNdS value"))

print("Cancer genes not in dNdS.by.gene database (S or N = 0, or mappability bad)")
cancer.genes[!(cancer.genes$Ensembl.Gene.ID %in% dNdS.by.gene$Ensembl.Gene.ID),.(Gene.Symbol,Entrez.Gene.ID,Ensembl.Gene.ID)]

##### annotate vogelstein cancer genes ######
cancer.genes.vogelstein <- fread("data/raw/vogelstein_driver_genes.tdv", header=TRUE)

# Correct Names
cancer.genes.vogelstein[gene_name=="FAM123B",gene_name:="AMER1"]
cancer.genes.vogelstein[gene_name=="MLL2",gene_name:="KMT2D"]
cancer.genes.vogelstein[gene_name=="MLL3",gene_name:="KMT2C"]

dNdS.by.gene$cancer.gene.vogelstein <- dNdS.by.gene$Associated.Gene.Name %in% cancer.genes.vogelstein$gene_name

print(paste(nrow(dNdS.by.gene[cancer.gene.vogelstein==TRUE]),"vogelstein cancer genes with dNdS value"))

# Which genes in vogelstein do not map to dNdS.by.gene database
print("Genes in vogelstein list not matched to dNdS value:")
cancer.genes.vogelstein[!(cancer.genes.vogelstein$gene_name %in% dNdS.by.gene$Associated.Gene.Name),gene_name]

# For all rows in dNdS.by.gene, add classification from cancer.genes.vogelstein if possible
setkey(cancer.genes.vogelstein,gene_name)
setkey(dNdS.by.gene,Associated.Gene.Name)
dNdS.by.gene$classification <- cancer.genes.vogelstein[dNdS.by.gene,classification]


##### add expression data ######
if(file.exists("data/RNAseq.by.gene.tsv")){
RNAseq <- fread("data/RNAseq.by.gene.tsv", header=TRUE)
setnames(RNAseq,c("Associated.Gene.Name","mean.expression","mean.of.stdev.of.expression","expression.ranking","expression.percent.rank","log10.expression"))
setkey(RNAseq,Associated.Gene.Name)

# Which dNdS.by.gene gene names not in RNAseq database
print(paste(length(dNdS.by.gene[!(dNdS.by.gene$Associated.Gene.Name %in% RNAseq$Associated.Gene.Name),.(Associated.Gene.Name)]$Associated.Gene.Name),"genes not in RNAseq database"))

# For all rows in dNdS.by.gene, add expression from RNAseq if avaliable
dNdS.by.gene <- RNAseq[dNdS.by.gene,]
}

# Summary Statistics
# how many 
print(paste(nrow(dNdS.by.gene),"entries (genes)"))
print(paste(nrow(dNdS.by.gene[total.mut>37]),"genes have 80% power (>37 mutations)"))
#11,407 / 17,929

# mean effect size
print(paste(dNdS.by.gene[q.value<0.1 & odds.ratio>1,.(mean(log2.odds.ratio))],"mean log2 odds ratio"))

# how many sig genes
print(paste(nrow(dNdS.by.gene[q.value<0.1 & odds.ratio>1]),"genes with q<0.1 positive selection"))
print(paste(nrow(dNdS.by.gene[q.value<0.1 & odds.ratio<1]),"genes with q<0.1 negative selection"))

# how many known cancer genes
print(paste(sum(dNdS.by.gene[q.value<0.1 & odds.ratio>1,cancer.gene.vogelstein]),"vogelstein cancer genes with q<0.1"))

#########################
###### Plot Graphs ######
#########################

library(ggplot2)
library(ggrepel)

#p-value correction plot
archive.file("data/p_correction.pdf")
pdf("results/p_correction.pdf",width=16, height=9, onefile = TRUE)
ggplot(dNdS.by.gene, aes(p.value,q.value)) + 
	geom_point(size=0.5) +
	theme(text = element_text(size=17)) +
	labs(x="P-value",y="q value") + 
	geom_abline(slope = 1,col="red")
dev.off()

# Funnel plot
pdf("results/funnel.pdf", width=16, height=9, onefile = TRUE)
ggplot(dNdS.by.gene, aes(total.mut,odds.ratio)) + 
	geom_point(aes(colour = is.significant),alpha=0.5,size=0.5) + 
	theme(legend.position="bottom",text=element_text(size = 17)) + 
	scale_colour_manual(name="significant? (q<0.1)",  values =c("black", "red")) + 
	scale_y_log10() + 
	scale_x_log10() + 
	labs(x="Total number of mutations per gene",y="Odds Ratio") + 
	geom_label_repel(data = dNdS.by.gene[order(p.value)][1:20], aes(label = Associated.Gene.Name), size = 3, box.padding = unit(0.5, "lines"), point.padding = unit(0.1, "lines"), force=1,segment.size=0.2,segment.color="blue")
dev.off()

# Annotated q value Volcano Plot
ggplot(dNdS.by.gene, aes(log2.odds.ratio,abs(log(q.value,10)))) + 
	geom_point(aes(colour = cancer.gene),alpha=0.3) + 
	scale_colour_manual(name="In COSMIC cancer gene census?",
		values =c("black", "red")) + 
	theme(legend.position="bottom",
		text=element_text(size = 17)) + 
	scale_y_log10() + 
	labs(x=bquote(log[2]*'(n:s/N:S) / Fold Change'),
		y=bquote('Negative-'*log[10]~'q-values')) + 
	geom_segment(color = "red",linetype=2,aes(x=-5,xend=5,y=1,yend=1)) + 
	geom_label_repel(data = dNdS.by.gene[order(p.value)][1:20], 
		aes(label = Associated.Gene.Name), 
		size = 3, 
		box.padding = unit(0.5, "lines"), 
		point.padding = unit(0.1, "lines"), 
		force=1,
		segment.size=0.2,
		segment.color="blue")
dev.off()

# Annotated P value Volcano Plot
pdf("results/volcano.pdf", width=16, height=9, onefile = TRUE)
ggplot(dNdS.by.gene, aes(log2.odds.ratio,abs(log(p.value,10)))) + 
	geom_point(aes(colour = cancer.gene),alpha=0.3,size=0.5) + 
	scale_colour_manual(name="In COSMIC cancer gene census?",
		values =c("black", "red")) + 
	theme(legend.position="bottom",
		text=element_text(size = 17)) + 
	scale_y_log10() + 
	labs(x=bquote(log[2]*'(Odds Ratio)'),
		y=bquote('Negative'*log[10]~'(P-value)')) + 
	geom_label_repel(data = dNdS.by.gene[order(p.value)][1:20],
		aes(label = Associated.Gene.Name),
		size = 3,
		box.padding = unit(0.5, "lines"),
		point.padding = unit(0.1, "lines"),
		force=1,
		segment.size=0.2,
		segment.color="blue")
dev.off()


# Annotated P value Volcano Plot by primary site
pdf("results/volcano_bysite.pdf", width=16, height=9, onefile = TRUE)
plot.volcano <- function(data){
ggplot(data, aes(log2.odds.ratio,abs(log(p.value,10)))) + 
	geom_point(alpha=0.3,size=0.5) + 
	theme(legend.position="bottom") + 
	labs(x=bquote(log[2]*'(n:s/N:S) / Fold Change'),y=bquote('Negative-'*log[10]~'P-values'),title=paste("Volcano Plot",unique(data$primary_site))) +
	scale_y_log10(limits=c(NA,65)) + 
	geom_label_repel(data = data[order(p.value)][1:20], aes(label = Associated.Gene.Name), size = 2, box.padding = unit(0.5, "lines"), point.padding = unit(0.1, "lines"), force=1,segment.size=0.2,segment.color="blue")
}

# plot volcano per primary site
for (site in unique(dNdS.by.gene.bysite$primary_site)){
	print(plot.volcano(dNdS.by.gene.bysite[primary_site==site]))
}
dev.off()


# COMBINE
most.signif <- dNdS.by.gene[dNdS.by.gene[, .I[which.min(p.value)], by=Ensembl.Gene.ID]$V1]
# but only positive...?

# all genes using most significant from individual types
pdf(width=16, height=9, onefile = TRUE)
ggplot(most.signif, aes(log2.odds.ratio,abs(log(p.value,10)))) +
	geom_point(aes(colour = cancer.gene),alpha=0.3,size=0.5) +
	scale_colour_manual(name="In COSMIC cancer gene census?",  values =c("black", "red")) +
	theme(legend.position="bottom") +
	scale_y_log10(limits=c(0.05,NA)) +
	labs(x=bquote(log[2]*'(n:s/N:S) / Fold Change'),y=bquote('Negative-'*log[10]~'P-values'),title="Volcano Plot") +
	geom_segment(color = "red",linetype=2,aes(x=0.85,xend=5,y=2.6,yend=2.6)) +
	geom_label_repel(data = most.signif[p.value<(0.0015) & log(odds.ratio)>0.6], aes(label = Associated.Gene.Name), size = 2, box.padding = unit(0.5, "lines"), point.padding = unit(0.1, "lines"), force=1,segment.size=0.2,segment.color="blue")
dev.off()

# join individual and pancancer
signif.pancancer <- dNdS.by.gene
setkey(signif.pancancer,Associated.Gene.Name)
setkey(most.signif,Associated.Gene.Name)
merge.signif <- merge(signif.pancancer,most.signif,all=TRUE)
stopifnot(nrow(merge.signif[cancer.gene.x!=cancer.gene.y])==0)

# pancancer significance vs subtype significance
pdf(width=16, height=9, onefile = TRUE)
ggplot(merge.signif,aes(abs(log(p.value.y,10)),abs(log(p.value.x,10)))) + 	geom_point(aes(colour = cancer.gene.y)) + 
	scale_y_log10() + scale_x_log10() +
	scale_colour_manual(name="In COSMIC cancer gene census?",  values =c("black", "red")) +
	theme(legend.position="bottom") + 
	geom_label_repel(data = merge.signif[abs(log(p.value.y,10))>(9) | abs(log(p.value.x,10))>(9)], aes(label = Associated.Gene.Name), size = 2, box.padding = unit(0.5, "lines"), point.padding = unit(0.1, "lines"), force=1,segment.size=0.2,segment.color="blue")
dev.off()

#####################
# TODO: CONFIDENCE INTERVALS?

print(paste("All genes in top reigon:",nrow(dNdS.by.gene[log(q.value)<(-6) & log(odds.ratio)>0.6])))
print(paste("Total Vogelstein cancer genes:",sum(dNdS.by.gene$cancer.gene.vogelstein)))
print(paste("Vogelstein cancer genes in top reigon:",sum(dNdS.by.gene[log(q.value)<(-6) & log(odds.ratio)>0.6]$cancer.gene.vogelstein)))
print(paste("Total cancer genes:",sum(dNdS.by.gene$cancer.gene)))
print(paste("Cancer genes in top reigon:",sum(dNdS.by.gene[log(q.value)<(-6) & log(odds.ratio)>0.6]$cancer.gene)))


# S Histogram (Background mutation rate power)
archive.file("results/S_distribution.pdf")
pdf("results/S_distribution.pdf", width=15, height=8, onefile = TRUE)
ggplot(dNdS.by.gene, aes(S)) + 
	geom_histogram(binwidth = 1) + 
	scale_x_continuous(limits = c(0, 100)) + 
	labs(x="mean synonymous mutations per gene",title="Distribution of S per gene")
dev.off()

# N Histogram
archive.file("results/N_distribution.pdf")
pdf("results/N_distribution.pdf", width=15, height=8, onefile = TRUE)
ggplot(dNdS.by.gene, aes(N)) + 
	geom_histogram(binwidth = 1) + 
	scale_x_continuous(limits = c(0, 100)) + 
	labs(x="mean nonsynonymous mutations per gene",title="Distribution of N per gene")
dev.off()

# Total mutation Histogram
archive.file("results/mutation_distribution.pdf")
pdf("results/mutation_distribution.pdf", width=15, height=8, onefile = TRUE)
ggplot(dNdS.by.gene, aes(total.mut)) + 
	geom_histogram(binwidth = 1) + 
	scale_x_continuous(limits = c(0, 200)) + 
	labs(x="Number of mutations per gene",y="Number of genes") +
	theme(legend.position="bottom",text = element_text(size=17))
dev.off()

# Total mutation hist for largest primary site
ggplot(dNdS.by.gene.bysite[primary_site=="Skin"], aes(total.mut)) + 
	geom_histogram(binwidth = 1) + 
	scale_x_continuous(limits = c(-1, 100)) + 
	labs(x="Number of mutations per gene",y="Number of genes") +
	theme(legend.position="bottom",text = element_text(size=17))
dev.off()

# P histogram
archive.file("results/pvalue_distribution.pdf")
pdf("results/pvalue_distribution.pdf", width=15, height=8, onefile = TRUE)
ggplot(dNdS.by.gene, aes(p.value)) +
	geom_histogram(aes(y =..density..),binwidth=0.005) +
	labs(x="P-value",y="Density") +
	geom_hline(yintercept = 1,linetype="dotted",color="red") + 
	theme(legend.position="bottom",text = element_text(size=17))
dev.off()

# P histogrm for largest primary site (Skin)
ggplot(dNdS.by.gene.bysite[primary_site=="Skin"], aes(p.value)) +
	geom_histogram(aes(y =..density..),binwidth=0.005) +
	labs(x="P-value",y="Density") +
	ylim(0,10) +
	theme(legend.position="bottom",text = element_text(size=17))
dev.off()

# P-value QQ plot
# http://GettingGeneticsDone.blogspot.com/
# See http://gettinggeneticsdone.blogspot.com/p/copyright.html

# Define the function
ggd.qqplot = function(pvector, main=NULL, ...) {
    o = -log10(sort(pvector,decreasing=F))
    e = -log10( 1:length(o)/length(o) )
    plot(e,o,pch=19,cex=1, main=main, ...,
        xlab=expression(Expected~~-log[10](italic(p))),
        ylab=expression(Observed~~-log[10](italic(p))),
        xlim=c(0,max(e)), ylim=c(0,4)) # max(o)
    lines(e,e,col="red")
}

# Using the ggd.qqplot() function
archive.file("results/pvalule_qqplot.pdf")
pdf("results/pvalule_qqplot.pdf", width=15, height=8, onefile = TRUE)
ggd.qqplot(dNdS.by.gene$p.value)
# END http://GettingGeneticsDone.blogspot.com/
dev.off()


# Order chromosomes
dNdS.by.gene$chromosome <- factor(dNdS.by.gene$chromosome, levels = c(1:22,"X","Y","MT"))

# Graphs of dNdS by chromosome position
archive.file("results/QC_bychromosome.pdf")
pdf("results/QC_bychromosome.pdf", width=15, height=8, onefile = TRUE)
ggplot(dNdS.by.gene, aes(chromosome.start, log2.odds.ratio)) + 
	geom_point(alpha=0.3,size=0.05) + 
	facet_wrap(~chromosome,scales="free_x")

# Graphs of P value by chromosome position
ggplot(dNdS.by.gene, aes(chromosome.start, abs(log(q.value,10)))) + 
	geom_point(alpha=0.3,size=0.05) + 
	facet_wrap(~chromosome,scales="free_x") + 
	scale_y_log10()
dev.off()

# Overall Dist
archive.file("results/OR_distribution.pdf")
pdf("results/OR_distribution.pdf", width=15, height=8, onefile = TRUE)
ggplot(dNdS.by.gene, aes(odds.ratio)) + 
	geom_histogram(binwidth = 0.01) + 
	geom_vline(xintercept = 1,color = "red") + 
	theme_grey(base_size = 17) + 
	labs(x=bquote(Log[2]~'Odds Ratio'),y="Number of genes") + 
	scale_x_log10()
dev.off()

# Cancer vs normal density dist
archive.file("results/cancer_distribution.pdf")
pdf("results/cancer_distribution.pdf", width=15, height=8, onefile = TRUE)
ggplot(dNdS.by.gene[S>3], aes(x=odds.ratio,fill=cancer.gene.vogelstein)) + 
	geom_density(alpha=0.3) + 
	geom_vline(xintercept = 1,color = "red") + 
	theme_grey(base_size = 30) + 
	labs(x="mean dN/dS per gene",title="Distribution of known cancer genes") + 
	scale_x_log10() + 
	theme(legend.position="bottom")

# Cancer gene dist TSG vs Onco vs normal
ggplot(dNdS.by.gene[S>3], aes(x=odds.ratio,fill=classification)) + 
	geom_density(alpha=0.3) + 
	geom_vline(xintercept = 1,color = "red") + 
	theme_grey(base_size = 30) + 
	labs(x="mean dN/dS per gene",title="Distribution of known cancer genes") + 
	scale_x_log10() + 
	theme(legend.position="bottom")
dev.off()

# Bottom 100 Hits (by p-value)
bottom <- dNdS.by.gene[order(p.value)][odds.ratio<1][1:100]
#,.(Associated.Gene.Name,S,N,odds.ratio,log2.odds.ratio,p.value,q.value,cancer.gene,Ensembl.Transcript.ID,expression.percent.rank,cancer.gene.vogelstein)][1:100]


# Top 100 Hits (by p-value)
top <- dNdS.by.gene[order(p.value)][odds.ratio>1][1:100]
#,.(Associated.Gene.Name,S,N,odds.ratio,log2.odds.ratio,p.value,q.value,cancer.gene,Ensembl.Transcript.ID,expression.percent.rank,cancer.gene.vogelstein)][1:100]

# Save whole table
archive.file("results/dNdS.tsv")
write.table(dNdS.by.gene, "results/dNdS.tsv", sep="\t",row.names=FALSE,quote=FALSE)

# Export top and bottom 100 hits to tsv file
archive.file("results/top_dNdS.tsv")
write.table(top, "results/top_dNdS.tsv", sep="\t",row.names=FALSE,quote=FALSE)
archive.file("results/bottom_dNdS.tsv")
write.table(bottom, "results/bottom_dNdS.tsv", sep="\t",row.names=FALSE,quote=FALSE)

sessionInfo()

sink(type="message")
sink()
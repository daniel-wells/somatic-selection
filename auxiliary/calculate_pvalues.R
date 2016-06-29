# Start writing to an output file
logfile <- file(paste("logs/calculate_dNdS.R.log",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "txt", sep = "."))
sink(logfile)
sink(logfile, type="message")

source("code/functions.R")

library(data.table)
random.dnds <- fread("data/random.aggregated.tsv")

# 187363055 
# Remove 9963 header rows
random.dnds <- random.dnds[transcript!="transcript"]
# 187353092
# = 9964 * 18803

# Convert types back to int
random.dnds[,job:=as.integer(job)]
random.dnds[,N:=as.integer(N)]
random.dnds[,S:=as.integer(S)]

random.dnds[is.na(N)]
# 4768
random.dnds[is.na(S)]
# 4768
random.dnds[is.na(S)|is.na(N)]
# 4768

# replace NA with 0
random.dnds[is.na(N),N:=0]
random.dnds[is.na(S),S:=0]

# Cache results
saveRDS(random.dnds,"data/random.dnds.rds")
# random.dnds <- readRDS("data/random.dnds.rds")
# library(data.table)

setkey(random.dnds,transcript)

# Calculate mean of N and S per transcript
dnds.stats <- random.dnds[,.("mean.N"=mean(N,na.rm = TRUE),"mean.S"=mean(S,na.rm = TRUE)),by=transcript]
dnds.stats[,random.ratio:=mean.N/mean.S]


dNdS.by.gene <- fread("data/dNdS_pancancer.tsv", header=TRUE)

setkey(dNdS.by.gene,Ensembl.Transcript.ID)
setkey(dnds.stats,transcript)


# for every row in random.dnds add S and N from dNdS.by.gene
dnds.stats2 <- dNdS.by.gene[dnds.stats,.(transcript,S,N,nonsynon.probability,synon.probability,mean.N,mean.S,random.ratio,Associated.Gene.Name)][!is.na(S)]
dnds.stats2[,observed.ratio:=N/S]
dnds.stats2[,calculated.ratio:=nonsynon.probability/synon.probability]


# Correlation between calculated and random odds
archive.file("results/random_vs_calculated.pdf")
pdf("results/random_vs_calculated.pdf", width=15, height=8, onefile = TRUE)
ggplot(dnds.stats2, aes(random.ratio,calculated.ratio)) + 	geom_point(alpha=0.1,size=0.5) +
	geom_abline(slope = 1,col="red",size=0.2) +
	theme_grey(base_size = 17) + 
	scale_x_log10(breaks=c(2:6)) +
	scale_y_log10(breaks=c(2:6)) +
	labs(x='Expected odds by shuffling',y='Expected odds by calculation')
dev.off()


# compare distributions of S and N seperately
archive.file("results/QC_randomised.pdf")
pdf("results/QC_randomised.pdf", width=15, height=8, onefile = TRUE)

ggplot(dnds.stats2, aes(mean.S)) + geom_histogram(binwidth=1) + xlim(0,100)
ggplot(dnds.stats2, aes(synon.probability)) + geom_histogram(binwidth=0.005)+ xlim(0,0.5)
ggplot(dnds.stats2, aes(S)) + geom_histogram(binwidth=1) + xlim(0,100)
ggplot(dnds.stats2, aes(mean.N)) + geom_histogram(binwidth=1) + xlim(0,150)
ggplot(dnds.stats2, aes(nonsynon.probability)) + geom_histogram(binwidth=0.005)+ xlim(0,0.65)
ggplot(dnds.stats2, aes(N)) + geom_histogram(binwidth=1) + xlim(0,150)
dev.off()


new <- lm(dnds.stats2$mean.N~dnds.stats2$nonsynon.probability)

ggplot(dnds.stats2, aes(mean.N,nonsynon.probability)) + geom_point(alpha=0.1,size=0.5) + xlim(0,500) + ylim(0,2.5) + stat_smooth(method = "lm", col = "red",size=0.1)
ggplot(dnds.stats2, aes(mean.S,synon.probability)) + geom_point(alpha=0.1,size=0.5) + xlim(0,200) + ylim(0,1) + stat_smooth(method = "lm", col = "red",size=0.1)
dev.off()

ggplot(dnds.stats2, aes(mean.N,mean.S)) + geom_point(alpha=0.1,size=0.5)+ xlim(0,800) + ylim(0,250)
ggplot(dnds.stats2, aes(nonsynon.probability,synon.probability)) + geom_point(alpha=0.1,size=0.5)+ xlim(0,3) + ylim(0,1)
dev.off()

dnds.stats2[,"total.mut":=S+N]
dnds.stats2[,"prob.s.e":=synon.probability/(nonsynon.probability+synon.probability)]
dnds.stats2[,"prob.s.r":=mean.S/(mean.N+mean.S)]
dnds.stats2[,"or.e":=(N/S)/(nonsynon.probability/synon.probability)]
dnds.stats2[,"or.r":=(N/S)/(mean.N/mean.S)]

print("Calculating P-values")
test <- function(x, p, n){binom.test(x, n, p, alternative="two.sided")$p.value}
dnds.stats2[,p.value.e:=mapply(test, S, prob.s.e, total.mut)]
dnds.stats2[,p.value.r:=mapply(test, S, prob.s.r, total.mut)]

cbind(dnds.stats2[order(p.value.r)][or.r<1,.(gene.name,p.value.r,or.r)][1:35],dnds.stats2[order(p.value.e)][or.r<1,.(gene.name,p.value.e,or.e)][1:35])

# p-value dist
ggplot(dnds.stats2, aes(p.value.r)) +
	geom_histogram(binwidth=0.005)

# volcano plot 
ggplot(dnds.stats2, aes(log(or.r,2),abs(log(p.value.r,10)))) +
	geom_point(alpha=0.3,size=0.5) +
	scale_y_log10()

# volcano plot
ggplot(dnds.stats2, aes(log(or.e,2),abs(log(p.value.e,10)))) +
	geom_point(alpha=0.3,size=0.5) +
	scale_y_log10()
dev.off()



ggplot(dnds.stats2, aes(abs(log(p.value.e,10)),abs(log(p.value.r,10)))) + geom_point(alpha=0.1,size=0.5) + geom_abline(slope = 1,col="red") +
	scale_y_log10() +
	scale_x_log10()
dev.off()

ggplot(dnds.stats2, aes(prob.s.r,prob.s.e)) + geom_point(alpha=0.1,size=0.5) + geom_abline(slope = 1,col="red")
dev.off()

# calculate dnds - SHOULD CALCULATE BASED ON RANDOM SAMPLE?
dnds.stats$random.dnds <- dnds.stats[,(mean.N/nonsynon.probability)/(mean.S/synon.probability)]
dnds.stats$observed.dnds <- dnds.stats[,(N/nonsynon.probability)/(S/synon.probability)]

# (n/N)/(s/S) = (n/s)/(N/S), so precalculate N/S
dnds.stats$NS <- dnds.stats[,(nonsynon.probability/synon.probability)]

setkey(dnds.stats,transcript)


hist(log(random.dnds2[transcript=="ENST00000336138",dnds]),breaks=300)


# annotate randomised data, with NS
random.dnds2 <- dnds.stats[random.dnds,.(transcript,i.S,i.N,NS,job)]

# Calculate dnds for each randomisation row
random.dnds2$dnds <- random.dnds2[,(i.N/i.S)/NS]
# plot overall dnds distribution
hist(log(random.dnds2$dnds),breaks=1000)
# plot dnds for specific gene
hist(log(random.dnds2[transcript=="ENST00000336138",dnds]),breaks=300)
hist(random.dnds2[transcript=="ENST00000336138",dnds],breaks=1000,xlim=c(0,1.5))


setkey(random.dnds2,transcript)

# function to calculate one tail (lower) p-value
calculate.pvalue <- function(transcript.id){
	calculated.dnds <- dnds.stats[transcript==transcript.id,observed.dnds]
	tsx <- random.dnds2[transcript==transcript.id]
return(nrow(tsx[dnds <= calculated.dnds]) / nrow(random.dnds2[transcript==transcript.id]))
}

# calculate p-values for all 18,000 transcripts
p.values <- sapply(dnds.stats$transcript,calculate.pvalue)
# 74 seconds
# 133 seconds??!

# plot histogram of p-values
hist(p.values,breaks=1000)

# replace 0 with minimum p-value (conservative)
p.values2 <- replace(p.values, p.values==0, 1/9964)

# adjust p-values for mutliple testing by Benjamini & Hochberg method
p.values.adj <- p.adjust(p.values2, method = "fdr")

hist(p.values.adj,breaks=1000,ylim=c(0,300))

p.adj <- data.table("p"=p.values.adj,"transcript"=names(p.values.adj))
p.adj[p==0]
# 148
p.adj[p<0.05]
# 314

hist(log(dnds.stats$random.expected),breaks=1000)
plot(log(dnds.stats$random.dnds),log(dnds.stats$observed.dnds),xlim=c(-2,2),ylim=c(-2,2))

plot(log(dnds.stats$synon.probability/dnds.stats$nonsynon.probability),log(dnds.stats$mean.S/dnds.stats$mean.N))

plot(dnds.stats$synon.probability/dnds.stats$nonsynon.probability,dnds.stats$mean.S/dnds.stats$mean.N)

library(ggplot2)

pdf(width=16, height=9, onefile = TRUE)
# randomised
hist(random.dnds[1:10000000][,S/N],breaks=10000,xlim=c(0,1.5))

# observed
hist(dnds.stats[,S/N],breaks=7000,xlim=c(0,1.5))
dev.off()

pdf(width=16, height=9, onefile = TRUE)
hist(dnds.stats[,S/N],breaks=2000,xlim=c(0,1.5))
# calculated
hist(dnds.stats[,synon.probability/nonsynon.probability],breaks=100,xlim=c(0,1.5))
dev.off()


pdf(width=16, height=9, onefile = TRUE)
ggplot(dnds.stats, aes(S.dev)) + geom_histogram(binwidth = 0.1) + scale_x_continuous(limits = c(-20, 10)) + labs(x="deviation from expected by random in units of standard deviation",title="Distribution of deviation of expected synon")

cor(dnds.stats$mean.S,dnds.stats$S)

ggplot(dnds.stats, aes(mean.S,S)) + geom_point(alpha=0.1) + scale_x_continuous(limits = c(0, 250)) + scale_y_continuous(limits = c(0, 250)) + labs(x="random S mean",title="random vs actual S") + geom_abline(slope = 1,col="red")

ggplot(dnds.stats, aes(N.dev)) + geom_histogram(binwidth = 0.1) + scale_x_continuous(limits = c(-20, 10)) + labs(x="deviation from expected by random in units of standard deviation",title="Distribution of deviation of expected nonsynon")

ggplot(dnds.stats, aes(mean.N,N)) + geom_point(alpha=0.1) + scale_x_continuous(limits = c(0, 250)) + scale_y_continuous(limits = c(0, 250)) + labs(x="random N mean",title="random vs actual N") + geom_abline(slope = 1,col="red")
dev.off()


archive.file("data/p.values.tsv")
write.table(p.adj, "data/p.values.tsv", sep="\t", row.names=FALSE, quote=FALSE)

sessionInfo()

sink(type="message")
sink()
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
random.dnds <- readRDS("data/random.dnds.rds")
# library(data.table)

setkey(random.dnds,transcript)

# Calculate mean and sd of N and S per transcript
dnds.stats <- random.dnds[,.("mean.N"=mean(N,na.rm = TRUE),"sd.N"=sd(N, na.rm = FALSE),"mean.S"=mean(S,na.rm = TRUE),"sd.S"=sd(S, na.rm = FALSE)),by=transcript]

cds <- fread("data/dNdS_by_transcript.tsv", header=TRUE)

setkey(cds,Ensembl.Transcript.ID)
setkey(dnds.stats,transcript)

# for every row in random.dnds add S and N from cds
dnds.stats <- cds[dnds.stats,.(transcript,S,N,nonsynon.probability,synon.probability,mean.N,mean.S,sd.N,sd.S)]

# calculate dnds - SHOULD CALCULATE BASED ON RANDOM SAMPLE?
dnds.stats$random.dnds <- dnds.stats[,(mean.N/nonsynon.probability)/(mean.S/synon.probability)]
dnds.stats$observed.dnds <- dnds.stats[,(N/nonsynon.probability)/(S/synon.probability)]

# (n/N)/(s/S) = (n/s)/(N/S), so precalculate N/S
dnds.stats$NS <- dnds.stats[,(nonsynon.probability/synon.probability)]

setkey(dnds.stats,transcript)


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
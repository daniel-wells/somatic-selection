source("code/functions.R")
# Start writing to an output file
logfile <- file(paste("logs/power_analysis.log",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "txt", sep = "."))
sink(logfile)
sink(logfile, type="message")

calculate.power <- function(total.mut,LORs,graph){
total <- NULL
for (LOR in LORs){
# average prob of s
prob.s <- 0.28

# given number of synonymous observed, and odds ratio, calculate new n/s, and total mutations,
# given (n/s) / (N/S) = odds raito, n is only unknown
ratio.n.s = 2^(LOR) * (1-0.28)/0.28
# total mut = n + s
total.mut = obs.s * ratio.n.s + obs.s
new.prob.s = obs.s / total.mut

# Calculate no of mutations X which has P-value P(X|H0)<0.0005, two tailed alpha=0.001 as coded critical value
crit.S.l <- qbinom(0.0005,size=(total.mut),prob=prob.s)-1
crit.S.u <- qbinom(0.9995,size=(total.mut),prob=prob.s)
actual.p.l <- pbinom(crit.S.l,total.mut,prob.s)
actual.p.u <- 1-pbinom(crit.S.u,total.mut,prob.s)
total.p <- actual.p.l + actual.p.u

# Calculate new S for H1 wich has effect size of |log(OR,2)|=0.5
# (n/s) / (N/S) = 2 | 0.5
# S/(S+N)=0.28 , S+N=total.mut
S = prob.s * total.mut
N = total.mut - S

# new ratios of n/s (= r1 & r1)
new.ratio = 2^(LOR) * N/S

# find new s / prob.s
# sub n=total.mut-s & rearange

s = total.mut/(new.ratio+1)
new.prob.s = s/total.mut # = 1/(r1+1)

if (graph==TRUE){
# number of possible S  
all.S <- seq(0,total.mut,by=1)

# distribution of number of S, given total number of mutations in gene 50
y <- dbinom(all.S,total.mut,prob.s)
plot(all.S,y,ylim=c(0,0.25),type='o')
abline(v=crit.S.l,col="red")
abline(v=crit.S.u,col="red")
lines(0:total.mut,dbinom(all.S,total.mut,new.prob.s),col="blue",type='o')
legend('topright',lty=1,c("Null Hypothesis","Data"),col=c("black","blue"))
}

# probability of rejecting H0 if H1 is true
if(LOR>0){
power <- pbinom(crit.S.l,total.mut,new.prob.s) + 1-pbinom(crit.S.u,total.mut,new.prob.s)
}else{
power <- 1-pbinom(crit.S.u,total.mut,new.prob.s) + pbinom(crit.S.l,total.mut,new.prob.s)
}

# prob less than or equal to crit.S1
new <- cbind(LOR,new.ratio,crit.S.l,crit.S.u,actual.p.l,actual.p.u,total.p,total.mut,power,new.prob.s)
total <- rbind(new,total)
}
return(total)
}

# plot example graph
archive.file("results/binomial_distribution.pdf")
pdf("results/binomial_distribution.pdf", width=11.7, height=8.3, onefile = TRUE)
calculate.power(50,-1,TRUE) # 2fold
dev.off()

# dont plot graph, but calculate for ranges of total.mut and effect sizes
power <- calculate.power(0:200,c(-2,-1,1,2),FALSE)

library(data.table)
data.t <- data.table(power)
data.t[,fold.change:=signif(2^LOR,3)]
data.t[,log2.odd.ratio:=as.factor(LOR)]

archive.file("data/power_analysis.tsv")
write.table(data.t, "data/power_analysis.tsv", sep="\t", row.names=FALSE, quote=FALSE)

library(ggplot2)
archive.file("results/power_curve.pdf")
pdf("results/power_curve.pdf", width=15, height=8, onefile = TRUE)
ggplot(data.t,aes(total.mut,group=fold.change,color=log2.odd.ratio)) +
  geom_line(aes(y=power)) +
  ylim(0,1) +
  theme(legend.position="bottom",text = element_text(size=17)) + 
  scale_colour_manual(name=bquote(Log[2]~'Odds Ratio'),  values =c("navy", "lightskyblue","peru","orangered")) +
  labs(x="Total number of mutations in a gene",y="Statistical Power")
dev.off()


sessionInfo()
sink(type="message")
sink()
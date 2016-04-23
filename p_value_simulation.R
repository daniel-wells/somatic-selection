library(data.table)

# Simulate p-value distribution
p.distribution <- function(n,prob.s){
# sample total mutations as chisq dist
total.mut <- round(rchisq(n, 1.9, ncp = 5)*15)
# sample S as binom dist
data <- data.table("S"=rbinom(n,total.mut,prob.s),"total.mut"=total.mut)
data[,N:=total.mut-S]
# remove those with no mutations at all
data <- data[total.mut!=0]
# calculate p-values
test <- function(x, p, n){binom.test(x, n, p, alternative="two.sided")$p.value}
data[,p.value:=mapply(test, S, prob.s, total.mut)]
return(data)
}

data <- p.distribution(18000,0.28)

# p-value histogram
hist(data$p.value,breaks=200)

# which are the highest p-vaues
data[p.value>0.95 & p.value!=1]

# histogram of total mutation number
hist(round(rchisq(2000, 1.9, ncp = 5)*15),breaks=200)

# for skin
hist(round(rexp(20000, rate = 0.05)),breaks=200)
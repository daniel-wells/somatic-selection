
##### ICGC Project List and Summary Statistics

#install.packages("jsonlite")
library(jsonlite)
library(data.table)
# from https://dcc.icgc.org/projects/details, Orange button in top right, Export Table as TSV
icgc <- fromJSON("https://dcc.icgc.org/api/v1/projects?filters={}&from=1&size=100", flatten=TRUE)
icgc <- data.table(icgc$hits)[,.("Project Code"=id,"Project Name"=name,"Primary Site"=primarySite,"Donors (DCC)"=totalLiveDonorCount,"Donors (All)"=totalDonorCount,"SSM"=ssmTestedDonorCount,	"CNSM"=cnsmTestedDonorCount,	"STSM"=stsmTestedDonorCount,	"SGV"=sgvTestedDonorCount,	"METH-A"=methArrayTestedDonorCount,	"METH-S"=methSeqTestedDonorCount,	"EXP-A"=expArrayTestedDonorCount,	"EXP-S"=expSeqTestedDonorCount,	"PEXP"=pexpTestedDonorCount,	"miRNA-S"=mirnaSeqTestedDonorCount,	"JCN"=jcnTestedDonorCount)]
write.table(icgc, file = "data/raw/ICGC_projects.tsv", sep = "\t",row.names = FALSE,quote = FALSE)

##### Vogelstein Cancer Driver Gene List

#install.packages("gdata")
library(gdata)
# from supplementary of http://doi.org/10.1126/science.1235122
vogelstein <- read.xls("http://science.sciencemag.org/content/sci/suppl/2013/03/27/339.6127.1546.DC1/1235122TablesS1-4.xlsx", sheet=6,pattern="Gene Symbol")
setnames(vogelstein,c("Gene.Symbol","Classification.","Ocogene.score..","Tumor.Suppressor.Gene.score..","Process","X..Mutated.Tumor.Samples.."),c("gene_name","classification","oncogene_score","tsg_score","process","tumour_sample_count"))
vogelstein$X <- NULL
write.table(vogelstein[1:125,],file = "data/raw/vogelstein_driver_genes.tdv", sep = "\t",row.names = FALSE,quote = FALSE)

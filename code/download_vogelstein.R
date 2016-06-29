
##### Vogelstein Cancer Driver Gene List
library(data.table)
#install.packages("gdata")
library(gdata)
# from supplementary of http://doi.org/10.1126/science.1235122
vogelstein <- read.xls("http://science.sciencemag.org/content/sci/suppl/2013/03/27/339.6127.1546.DC1/1235122TablesS1-4.xlsx", sheet=6,pattern="Gene Symbol")
setnames(vogelstein,c("Gene.Symbol","Classification.","Ocogene.score..","Tumor.Suppressor.Gene.score..","Process","X..Mutated.Tumor.Samples.."),c("gene_name","classification","oncogene_score","tsg_score","process","tumour_sample_count"))
vogelstein$X <- NULL
write.table(vogelstein[1:125,],file = "data/raw/vogelstein_driver_genes.tdv", sep = "\t",row.names = FALSE,quote = FALSE)
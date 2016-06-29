
##### ICGC Project List and Summary Statistics

#install.packages("jsonlite")
library(jsonlite)
library(data.table)
# from https://dcc.icgc.org/projects/details, Orange button in top right, Export Table as TSV
icgc <- fromJSON("https://dcc.icgc.org/api/v1/projects?filters={}&from=1&size=100", flatten=TRUE)
icgc <- data.table(icgc$hits)[,.("Project Code"=id,"Project Name"=name,"Primary Site"=primarySite,"Donors (DCC)"=totalLiveDonorCount,"Donors (All)"=totalDonorCount,"SSM"=ssmTestedDonorCount,	"CNSM"=cnsmTestedDonorCount,	"STSM"=stsmTestedDonorCount,	"SGV"=sgvTestedDonorCount,	"METH-A"=methArrayTestedDonorCount,	"METH-S"=methSeqTestedDonorCount,	"EXP-A"=expArrayTestedDonorCount,	"EXP-S"=expSeqTestedDonorCount,	"PEXP"=pexpTestedDonorCount,	"miRNA-S"=mirnaSeqTestedDonorCount,	"JCN"=jcnTestedDonorCount)]
setnames(icgc, c("Project Code","Primary Site"), c("project_code","primary_site"))

# Filter out projects with a publication moratorium
# https://ocg.cancer.gov/programs/target/target-publication-guidelines
# http://cancergenome.nih.gov/publications/publicationguidelines
# http://docs.icgc.org/portal/publication/#current-moratorium-status-for-icgc-projects
icgc <- icgc[!project_code %in% c("SKCA-BR","COCA-CN","LUSC-CN","LIAD-FR","BOCA-FR","LIHM-FR","PACA-IT","BTCA-JP","LICA-CN","LAML-CN","WT-US")]
icgc <- icgc[SSM!=0]

write.table(icgc, file = "data/raw/ICGC_projects.tsv", sep = "\t",row.names = FALSE,quote = FALSE)

# write file of urls
urls <- paste0("https://dcc.icgc.org/api/v1/download?fn=/release_21/Projects/",icgc$project_code,"/simple_somatic_mutation.open.",icgc$project_code,".tsv.gz")

write.table(urls,"data/url-list.txt",sep = "",row.names = FALSE,col.names = FALSE,quote = FALSE)
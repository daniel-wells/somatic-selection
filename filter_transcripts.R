# create list of final transcript ID's actually used in the finial analysis so it can be used in previous steps to filter transcripts so as to not process more than required

library(data.table)
cds <- fread("data/dNdS_by_transcript.tsv", header=TRUE)
print(paste(nrow(cds),"transcripts dNdS values loaded."))

setnames(cds,"gene.id","gene")
setnames(cds,"Associated.Gene.Name","gene.name")
setnames(cds,"cds.length","cds_length")

# For each gene, find max cds length
max.cds.by.gene <- unique(cds[,.(max.cds = max(cds_length,na.rm=TRUE)),by=gene])

# For all rows in gene.maxlength, add cds values where cds_length=max(cds_length) by gene, for each gene (multiples as sometimes multiple transcript have max cds length)
setkey(max.cds.by.gene,gene,max.cds)
setkey(cds,gene,cds_length)
# cds = 57162
dNdS.by.gene <- cds[max.cds.by.gene,]
# 23,673

paste(nrow(dNdS.by.gene[duplicated(dNdS.by.gene)==TRUE,]),"duplicate entries (same gene and cds length) after getting dNdS from canonical (longest) transcript")


# set key as all columns
setkey(dNdS.by.gene,gene)

# 18,803, strickt set
unique <- unique(dNdS.by.gene[,Ensembl.Transcript.ID,by=gene])
# or remove 'by=gene' for no strict set of 23,673

unique <- unique$Ensembl.Transcript.ID

saveRDS(unique,"data/final.transcript.list.rds",compress=FALSE)
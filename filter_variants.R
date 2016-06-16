################# Filter out poor mappability reigons ##################
# Start writing to an output file
logfile <- file(paste("logs/filter_variants.R.log",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "txt", sep = "."))
sink(logfile)
sink(logfile, type="message")

source("code/functions.R")
library(data.table)

observed_variants <- readRDS("data/coding.mutations.rds")

# "gene.id" -> "gene_affected"
# "transcript.id" -> "transcript_affected"
# "variant.class" -> "consequence_type"
setnames(observed_variants,c("icgc_mutation_id","icgc_donor_id","project_code","chromosome","chromosome_start","chromosome_end","chromosome_strand","mutation_type","reference_genome_allele","mutated_from_allele","mutated_to_allele","variant.class","aa_mutation","cds_mutation","gene.id","transcript.id","sequencing_strategy"))

# columns removed in load_mutations.R creating duplicate rows
setkey(observed_variants,icgc_mutation_id,icgc_donor_id,transcript.id)
observed_variants <- unique(observed_variants)
print(paste(nrow(observed_variants),"unique observed variants loaded"))


################ Load mappability file ################
# /mnt/lustre/users/bschuster/TCGA/Coverage/mappability_24_100bp_windows_exons.bed.gz
mappability <- fread("zcat < data/raw/mappability_100bp_windows_exons.bed.gz")
setnames(mappability,c("chromosome","start","mappability"))
# convert from 0 based to 1 based
mappability[,start:=start+1]
# add stop value
mappability[,stop:=start]
# reorder columns
mappability <- mappability[,.(chromosome, start, stop, mappability)]
setkey(mappability, chromosome, start, stop)
mappability[,mappability:=as.numeric(mappability)]

# make chromosome names consistent with BEDgraph file
observed_variants[,chromosomeC:=paste0("chr",chromosome)]
observed_variants[chromosomeC=="chrMT"]$chromosomeC <- "chrM"

################ Filter by transcript ################
final.transcripts <- readRDS("data/final.transcript.list.rds")
setkey(final.transcripts,Ensembl.Transcript.ID)
setkey(observed_variants,transcript.id)

final.variants <- final.transcripts[observed_variants,nomatch=0]
print(paste(nrow(final.variants),"variants in final transcripts"))

################ Only unique variants need to be annotated with mappability
setkey(final.variants, chromosomeC, chromosome_start, chromosome_end)
uniq_obs_var <- unique(final.variants)
uniq_obs_var <- uniq_obs_var[,.(chromosomeC, chromosome_start, chromosome_end)]
print(paste(nrow(uniq_obs_var),"unique variants"))

################ Map variants to mappability scores (windowed join)
obs_var_map <- foverlaps(mappability,uniq_obs_var, type="any",nomatch=0)
# ~20 seconds
print(paste(nrow(obs_var_map),"variant-mappability pairs"))

################ remove duplicates
setkey(obs_var_map, chromosome, chromosome_start, chromosome_end)
print(paste(nrow(obs_var_map[duplicated(obs_var_map)]),"duplicated variants"))

# order by lowest mappability first
obs_var_map <- obs_var_map[order(chromosome,chromosome_start,chromosome_end, mappability),]

# choose lowest mappability score
setkey(obs_var_map,chromosome, chromosome_start, chromosome_end)
obs_var_map <- unique(obs_var_map)
print(paste(nrow(obs_var_map),"unique variant-mappability pairs"))


################ Join mappabilty values back with variants
setnames(obs_var_map,"chromosome","chromosomeC")

# for each row in observed_variants add mappability info
final.variants <- obs_var_map[final.variants]
print("rows with no mappability found:")
final.variants[is.na(mappability)]
print(paste(nrow(final.variants),"final variant-mappability pairs"))

# plot QC graph, mappability of final variants
archive.file("results/filter_variants_QC.pdf")
pdf("results/filter_variants_QC.pdf",width=9, height=9, onefile = TRUE)
hist(final.variants$mappability,breaks=2000,ylim=c(0,40000))
# more graphs below

################ Remove poor mappability reigons 
final.variants <- final.variants[mappability>0.79]
print(paste(nrow(final.variants),"variants left after removing those with mappability <0.8"))


################ Filter out common variants
exac <- fread('data/ExAC.bed')
setnames(exac,c("chromosome","chromosome.minus.1","chromosome_start","rsID","quality","ref","mutated_to_allele","freq"))
setkey(exac,chromosome,chromosome_start,mutated_to_allele)

final.variants[,chromosome:=as.character(chromosome)]
setkey(final.variants,chromosome,chromosome_start,mutated_to_allele)

# Add freqency information to each variant
final.variants <- exac[final.variants]

# plot QC graph, freq of variants in exac
hist(log(final.variants[!is.na(freq)]$freq,10))
dev.off()

# filter out variants with freq greater than 1% or not in exac
final.variants <- final.variants[freq<0.01 | is.na(freq)]
print(paste(nrow(final.variants),'variants left after removing those common in population (1%)'))

# TP53
# tp53.mutants <- unique(final.variants[Associated.Gene.Name=="TP53" & variant.class!="synonymous_variant",icgc_donor_id])
#
# final.variants <- final.variants[icgc_donor_id %in% tp53.mutants]

archive.file("data/coding.mutations.filtered.rds")
saveRDS(final.variants,"data/coding.mutations.filtered.rds")

# Save single base coding substitutions
single.base.coding.substitutions <- final.variants[mutation_type=="single base substitution"]
print(paste(nrow(single.base.coding.substitutions),"single base coding substitutions"))

# Save single base, coding mutations
archive.file("data/single.base.coding.substitutions.rds")
saveRDS(single.base.coding.substitutions,"data/single.base.coding.substitutions.rds",compress = FALSE)

sessionInfo()

sink(type="message")
sink()
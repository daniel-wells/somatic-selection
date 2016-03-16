# Start writing to an output file
logfile <- file(paste("logs/calculate_dNdS.R.log",format(Sys.time(), "%Y-%m-%d.%H-%M-%S"), "txt", sep = "."))
sink(logfile)
sink(logfile, type="message")

source("code/functions.R")
library(data.table)

# Load somatic mutation tsv
file_list <- list.files(path="data/raw/ICGC",pattern="simple_somatic_mutation.open.*.tsv.gz",full.names=TRUE)

# Filter as soon as reading or will take up too much memory (>64GB), MELA-AU alone is 34GB unzipped
read.and.filter <- function(x) {
print(paste("Reading file:",x))
data <- fread(x)
# all mutations, remove unused columns (~15% original size)
data <-  data[,.(icgc_mutation_id,icgc_donor_id,project_code,chromosome,chromosome_start,chromosome_end,chromosome_strand,mutation_type,reference_genome_allele,mutated_from_allele,mutated_to_allele,consequence_type,aa_mutation,cds_mutation,gene_affected,transcript_affected,sequencing_strategy)]
return(data)
}

# Load tsv mutations
mutations <- rbindlist(lapply(paste('zcat < ',file_list), read.and.filter))

# Set catagoricals as factors (saves 4GB)
mutations[, c("project_code","chromosome","chromosome_strand","mutation_type","consequence_type","sequencing_strategy","gene_affected","transcript_affected","icgc_donor_id") := lapply(mutations[,c("project_code","chromosome","chromosome_strand","mutation_type","consequence_type","sequencing_strategy","gene_affected","transcript_affected","icgc_donor_id"),with=FALSE],as.factor),with=FALSE]

# Save all mutations as cache
saveRDS(mutations,"data/simple.somatic.mutations.aggregated.rds",compress = FALSE)

# Save all coding substitutions
coding.mutations <- mutations[consequence_type %in% c("missense_variant", "synonymous_variant", "frameshift_variant","disruptive_inframe_deletion","disruptive_inframe_insertion","inframe_deletion","inframe_insertion","start_lost","stop_lost","stop_gained","exon_loss_variant")]

saveRDS(coding.mutations,"data/coding.mutations.rds",compress = FALSE)

# Save single base coding substitutions
single.base.coding.substitutions <- mutations[mutation_type=="single base substitution" & consequence_type %in% c("missense_variant", "synonymous_variant", "frameshift_variant","disruptive_inframe_deletion","disruptive_inframe_insertion","inframe_deletion","inframe_insertion","start_lost","stop_lost","stop_gained")]
#4,649,834

# Save single base, coding mutations
archive.file("data/single.base.coding.substitutions.rds")
saveRDS(single.base.coding.substitutions,"data/single.base.coding.substitutions.rds",compress = FALSE)

sessionInfo()

sink(type="message")
sink()
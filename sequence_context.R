# source("http://bioconductor.org/biocLite.R")
# biocLite("SomaticSignatures")
# biocLite("VariantAnnotation")
# biocLite("BSgenome")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")

library(SomaticSignatures)
library(data.table)

## Genomic sequences
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19

# Load somatic mutation tsv
file_list <- list.files(path="/mnt/lustre/users/dwells/data/raw/ICGC",pattern="simple_somatic_mutation.open.*.tsv.gz",full.names=TRUE)

# Get rid of 1.4GB (zipped!) MELA mutations
file_list <- file_list[-34]

all.simple.somatic.mutations <- rbindlist(lapply(paste('zcat < ',file_list[1:7]), fread))

all.simple.somatic.mutations[,.N,by=project_code]
all.simple.somatic.mutations[,.N,by=mutation_type]

# all single base pair, exonic mutations
ICGCraw <- all.simple.somatic.mutations[mutation_type=="single base substitution" & consequence_type %in% c("missense_variant", "synonymous_variant", "frameshift_variant","disruptive_inframe_deletion","disruptive_inframe_insertion","inframe_deletion","inframe_insertion","start_lost","stop_lost","stop_gained")]

# Remove duplicate annotations (1mut:>1annot)
setkey(ICGCraw,icgc_donor_id,icgc_mutation_id)
unique(ICGCraw)

# Make VRanges object
vr = VRanges(
	seqnames = ICGCraw$chromosome,
	ranges = IRanges(ICGCraw$chromosome_start,ICGCraw$chromosome_end),
	ref = ICGCraw$reference_genome_allele,
	alt = ICGCraw$mutated_to_allele,
	sampleNames = ICGCraw$icgc_donor_id,
	study = ICGCraw$project_code)

# add "chr" to work with UCSC.hg19
ucsc(vr)

## Annotate variants with context
vr_context <- mutationContext(vr, genome)

# for 1,191,475 variants
# user  system elapsed
# 6.637   7.223  15.405

# Way too slow ~4.5hrs
# RefGenome = FaFile("/mnt/lustre/data/GRCH37/GRCH37.fasta")
# ptm <- proc.time()
# vr_context = mutationContext(head(vr,n=100), RefGenome)
# proc.time() - ptm

motif.matrix.freq = motifMatrix(vr_context, group = "study", normalize = TRUE)
motif.matrix.count = motifMatrix(vr_context, group = "study", normalize = FALSE)

plotMutationSpectrum(vr_context, "study")
dev.off()

dput(motif.matrix.freq, file = "motif.matrix.freq.dput")
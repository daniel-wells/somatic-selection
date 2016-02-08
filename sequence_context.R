# source("http://bioconductor.org/biocLite.R")
# biocLite("SomaticSignatures")
# biocLite("VariantAnnotation")
# biocLite("BSgenome")
# biocLite("BSgenome.Hsapiens.UCSC.hg19")
# biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")

library(SomaticSignatures)
library(data.table)

## Genomic sequences
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19

# Extract trinucleotide frequency
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
cds = keepStandardChromosomes(reduce(cds(TxDb.Hsapiens.UCSC.hg19.knownGene)))
sum(width(cds)) # 35,226,001 - 35MB

k = 3
n = 1e6
k3_cds = kmerFrequency(genome, n, k, cds)

# Modified from SomaticSignatures/R/normalize.R
s = BStringSet(rownames(sca_mm))
base_motif = subseq(s, 4, 6)
subseq(base_motif, 2, 2) = subseq(s, 1, 1)
bs = as(base_motif, "character")
all(bs %in% names(k3_exons))
names(k3_exons[!names(k3_exons) %in% bs])
unique(bs) #32!
idx = match(bs, names(k3_cds))
# Modified here *2 so total frequency is 3 (as only 6/12 trinucleotides represented in mutation profile, and each is represented 3 times to make 96)
sss = 2 * as.vector(k3_cds[idx])
sum(sss) # =3
trimer.counts = unique(sss * 35226001)
names(trimer.counts) <- unique(bs)

# Load somatic mutation tsv
file_list <- list.files(path="/mnt/lustre/users/dwells/data/raw/ICGC",pattern="simple_somatic_mutation.open.*.tsv.gz",full.names=TRUE)

# Filter as soon as reading or will take up too much memory (>64GB), MELA-AU alone is 34GB unzipped
read.and.filter <- function(x) {
	print(paste("Reading file:",x))
	data <- fread(x)
	# all single base pair, exonic mutations
	data <- data[mutation_type=="single base substitution" & consequence_type %in% c("missense_variant", "synonymous_variant", "frameshift_variant","disruptive_inframe_deletion","disruptive_inframe_insertion","inframe_deletion","inframe_insertion","start_lost","stop_lost","stop_gained")]
	return(data)
}

# Load tsv mutations
single.base.coding.substitutions <- rbindlist(lapply(paste('zcat < ',file_list), read.and.filter))

# Save object
saveRDS(single.base.coding.substitutions, "data/single.base.coding.substitutions.rds")

single.base.coding.substitutions[,.N,by=project_code]
single.base.coding.substitutions[,.N,by=mutation_type]
single.base.coding.substitutions[,.N,by=consequence_type]
single.base.coding.substitutions[,.N,by=sequencing_strategy]

# Remove duplicate annotations (1mut:>1annot due to multiple transcripts)
setkey(single.base.coding.substitutions,icgc_donor_id,icgc_mutation_id)
ICGCraw <- unique(single.base.coding.substitutions)

# Make VRanges object
vr = VRanges(
	seqnames = ICGCraw$chromosome,
	ranges = IRanges(ICGCraw$chromosome_start,ICGCraw$chromosome_end),
	ref = ICGCraw$reference_genome_allele,
	alt = ICGCraw$mutated_to_allele,
	sampleNames = ICGCraw$icgc_donor_id,
	study = ICGCraw$project_code,
	sequencing_strategy = ICGCraw$sequencing_strategy)

# add "chr" to work with UCSC.hg19
vr <- ucsc(vr)

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
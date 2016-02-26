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

### LOAD ICGC MUTATION DATA ###
single.base.coding.substitutions <- readRDS("data/single.base.coding.substitutions.rds")

### CONTEXTUALISE MUTATIONS & CALCULATE MUTATION PROFILE & PROBAILITIES ###

single.base.coding.substitutions[,.N,by=project_code]
single.base.coding.substitutions[,.N,by=mutation_type]
single.base.coding.substitutions[,.N,by=consequence_type]
single.base.coding.substitutions[,.N,by=sequencing_strategy]

# Remove duplicate annotations (1mut:>1annot due to multiple transcripts)
# Warning, 1 base can mutate to 2 different bases in same patient - use mut_id not position
setkey(single.base.coding.substitutions,icgc_donor_id,icgc_mutation_id)
single.base.coding.substitutions <- unique(single.base.coding.substitutions)

# Make VRanges object
vr = VRanges(
	seqnames = single.base.coding.substitutions$chromosome,
	ranges = IRanges(single.base.coding.substitutions$chromosome_start,single.base.coding.substitutions$chromosome_end),
	ref = single.base.coding.substitutions$reference_genome_allele,
	alt = single.base.coding.substitutions$mutated_to_allele,
	sampleNames = single.base.coding.substitutions$icgc_donor_id,
	study = single.base.coding.substitutions$project_code,
	sequencing_strategy = single.base.coding.substitutions$sequencing_strategy)

# add "chr" to work with UCSC.hg19
vr <- ucsc(vr)

## Annotate variants with context
vr_context <- mutationContext(vr, genome)

motif.matrix.count = motifMatrix(vr_context, group = "study", normalize = FALSE)

# number of donors per project
setkey(single.base.coding.substitutions,icgc_donor_id)
ICGCdonors <- unique(single.base.coding.substitutions)
donor.count <- ICGCdonors[,.("donor.count"=.N),by=project_code][order(donor.count)]

CJ.dt = function(X,Y) {
  stopifnot(is.data.table(X),is.data.table(Y))
  k = NULL
  X = X[, c(k=1, .SD)]
  setkey(X, k)
  Y = Y[, c(k=1, .SD)]
  setkey(Y, NULL)
  X[Y, allow.cartesian=TRUE][, k := NULL][]
}

coding.trimer.counts <- readRDS("data/coding.trimer.counts.rds")

# cross join with donor counts
trimer.count.by.project <- CJ.dt(coding.trimer.counts,donor.count)
# overall counts of trinucleotides over all donors in a project
trimer.count.by.project$total.count <- trimer.count.by.project$coding.trimer.counts * trimer.count.by.project$donor.count


# Convert mutation motif counts to data table
motif.probabilities <- as.data.table(melt(motif.matrix.count))
setnames(motif.probabilities,c("mutation","project_code","mutation_count"))
# add base motif column
motif.probabilities$mutation <- as.character(motif.probabilities$mutation)
motif.probabilities$base_motif <- subseq(motif.probabilities$mutation, 4, 6)
subseq(motif.probabilities$base_motif, 2, 2) <- subseq(motif.probabilities$mutation, 1, 1)

# all rows in motif probabilities with timer counts added
setkey(trimer.count.by.project,project_code,base_motif)
setkey(motif.probabilities,project_code,base_motif)
motif.probabilities <- trimer.count.by.project[motif.probabilities]
# Calculate mutation "probability"
motif.probabilities$mutation.probability <- motif.probabilities$mutation_count / motif.probabilities$total.count


# number of somatic coding mutations per donor
hist(single.base.coding.substitutions[,.N,by=icgc_donor_id][order(N)]$N,breaks=3000,xlim=c(0,500))
dev.off()

w_df = melt(motif.matrix.count, varnames = c("motif", "sample"))
    w_df$alteration = sub("([ACGTN])([ACGTN]) .+", "\\1>\\2", w_df$motif)
    w_df$context = sub("[ACGTN][ACGTN] (.+)", "\\1", w_df$motif)

pdf(width=20,height=60)
ggplot(w_df) + geom_bar(aes_string(x = "context", y = "value"), stat = "identity", position = "identity") + facet_grid(sample ~ alteration,scales="free_y")
dev.off()

saveRDS(motif.probabilities, "data/motif.probabilities.rds")
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

if (file.exists("data/single.base.coding.substitutions.rds")){
	single.base.coding.substitutions <- readRDS("data/single.base.coding.substitutions.rds")
}else{

	# Load somatic mutation tsv
	file_list <- list.files(path="data/raw/ICGC",pattern="simple_somatic_mutation.open.*.tsv.gz",full.names=TRUE)

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
}

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

motif.matrix.count = motifMatrix(vr_context, group = "study", normalize = FALSE)

if (file.exists("data/coding.trimer.counts.rds")){
	coding.trimer.counts <- readRDS("data/coding.trimer.counts.rds")
}else{
	# Extract trinucleotide frequency
	library(TxDb.Hsapiens.UCSC.hg19.knownGene)
	cds = keepStandardChromosomes(reduce(cds(TxDb.Hsapiens.UCSC.hg19.knownGene)))

	k = 3
	n = 1e6
	k3_cds = kmerFrequency(genome, n, k, cds)

	# Modified from SomaticSignatures/R/normalize.R
	s = BStringSet(rownames(motif.matrix.count))
	base_motif = subseq(s, 4, 6)
	subseq(base_motif, 2, 2) = subseq(s, 1, 1)
	bs = as(base_motif, "character")
	all(bs %in% names(k3_cds))
	names(k3_cds[!names(k3_cds) %in% bs])
	unique(bs) #32!
	idx = match(bs, names(k3_cds))
	# Modified here *2 so total frequency is 3 (as only 6/12 trinucleotides represented in mutation profile, and each is represented 3 times to make 96)
	sss = 2 * as.vector(k3_cds[idx])
	sum(sss) # =3
	coding.trimer.counts = unique(sss * sum(width(cds)))  # 35,226,001 - 35Mbp coding exome length
	names(coding.trimer.counts) <- unique(bs)
	saveRDS(coding.trimer.counts, "data/coding.trimer.counts.rds")	
}

# number of donors per project
setkey(ICGCraw,icgc_donor_id)
ICGCdonors <- unique(ICGCraw)
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


# Convert coding.trimer.counts to data table and cross join with donor counts
trimer.count.by.project <- CJ.dt(data.table("base_motif"=names(coding.trimer.counts),coding.trimer.counts),donor.count)
# overall counts of trinucleotides over all donors in a project
trimer.count.by.project$total.count <- trimer.count.by.project$coding.trimer.counts * trimer.count.by.project$donor.count


# Convert mutation motif counts to data table
motif.probabilities <- as.data.table(melt(motif.matrix.count))
setnames(motif.probabilities,c("mutation","project_code","mutation_count"))
# add base motif column
motif.probabilities$mutation <- as.character(motif.probabilities$mutation)
motif.probabilities$base_motif <- subseq(motif.probabilities$mutation, 4, 6)
subseq(motif.probabilities$base_motif, 2, 2) <- subseq(motif.probabilities$mutation, 1, 1)

# all rows in motif matrix with timer counts added
setkey(trimer.count.by.project,project_code,base_motif)
setkey(motif.probabilities,project_code,base_motif)
motif.probabilities <- trimer.count.by.project[motif.probabilities]
# Calculate mutation "probability"
motif.probabilities$mutation.probability <- motif.probabilities$mutation_count / motif.probabilities$total.count


# number of somatic coding mutations per donor
hist(ICGCraw[,.N,by=icgc_donor_id][order(N)]$N,breaks=3000,xlim=c(0,500))
dev.off()

w_df = melt(motif.matrix.count, varnames = c("motif", "sample"))
    w_df$alteration = sub("([ACGTN])([ACGTN]) .+", "\\1>\\2", w_df$motif)
    w_df$context = sub("[ACGTN][ACGTN] (.+)", "\\1", w_df$motif)

pdf(width=20,height=60)
ggplot(w_df) + geom_bar(aes_string(x = "context", y = "value"), stat = "identity", position = "identity") + facet_grid(sample ~ alteration,scales="free_y")
dev.off()

saveRDS(motif.probabilities, "data/motif.probabilities.rds")
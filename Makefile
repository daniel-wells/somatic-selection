.PHONY : results all
results: results/figures.pdf results/dNdS_by_gene.tsv
all: raw_data results

# To print a variable type e.g. "make -n print-TGCA_RNAseq_data"
print-%  : ; @echo $* = $($*)

#################################
####### Download Raw Data #######
#################################

data/raw/ICGC_projects.tsv:
	# Manual Download list of projects -O data/raw/ICGC_projects.tsv
	echo https://dcc.icgc.org/projects/details https://dcc.icgc.org/78e07ae9-2c8c-4511-a2b0-77c71475688d

data/url-list.txt: data/raw/ICGC_projects.tsv
	for OUTPUT in $$(cut -f 1,7 data/raw/ICGC_projects.tsv | grep -Pv '(\t0|Project)' | cut -f 1) ; do \
		echo "https://dcc.icgc.org/api/v1/download?fn=/release_20/Projects/$$OUTPUT/simple_somatic_mutation.open.$$OUTPUT.tsv.gz" >> data/url-list.txt \ ; \
	done

URLS=$(shell awk '{printf "%s\n", $$1}' data/url-list.txt)

ICGC_project_mutation_files=$(addprefix data/raw/ICGC/,$(notdir $(URLS)))

$(ICGC_project_mutation_files): data/url-list.txt
	mkdir -p $(dir $@)
	curl "$(filter $(addprefix %/,$(notdir $@)),$(URLS))" -o $@
	# 2.6G in total
	#? find data/raw/ICGC -name "simple_somatic_mutation.open.*.tsv.gz" -exec touch {} \;

data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz:
	# Download ensembl cds nt sequence - referenced as ftp://ftp.ensembl.org/pub/release-75/mysql/homo_sapiens_core_75_37/ at https://docs.icgc.org/methods
	curl ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.75.cds.all.fa.gz -o data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz
	curl ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cds/CHECKSUMS -o data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz.CHECKSUMS
	curl ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cds/README -o data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz.README3
	sum data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz
	head -1 data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz.CHECKSUMS
	sha256sum data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz >> data/raw/SHA256SUMS


data/raw/mart_export.txt.gz:
	# Requires manual download
	http://www.ensembl.org/biomart/martview/a9103936cd932b57917528550c7c9a2b?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.external_gene_name|hsapiens_gene_ensembl.default.feature_page.gene_biotype|hsapiens_gene_ensembl.default.feature_page.transcript_biotype|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id&FILTERS=&VISIBLEPANEL=resultspanel
	sha256sum data/raw/mart_export.txt.gz >> data/raw/SHA256SUMS

data/raw/cancer_gene_census.csv:
	# Download cancer gene census for annotation in analysis
	# Requires Password
	sftp "daniel.wells@well.ox.ac.uk"@sftp-cancer.sanger.ac.uk
	get /files/grch38/cosmic/v75/cancer_gene_census.csv data/raw/cancer_gene_census.csv
	#572 genes

data/raw/vogelstein_driver_genes.tdv:
	# Requires manual download

$(TGCA_RNAseq_data):
	# Requires manual download

data/raw/HGNC.tsv:
	# Requires manual download

data/raw/exons.hg19.mappability100.bed.gz:
	# Requires manual download

raw_data: data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz data/raw/mart_export.txt.gz $(ICGC_project_mutation_files) data/raw/cancer_gene_census.csv data/raw/vogelstein_driver_genes.tdv data/raw/ICGC_projects.tsv data/raw/exons.hg19.mappability100.bed.gz data/raw/HGNC.tsv $(TGCA_RNAseq_data)


#############################
####### Main Analysis #######
#############################

#######
####### RNAseq
#######

TGCA_RNAseq_data=$(wildcard /users/bschuster/sharedscratch/Splicing/GeneExpRaw/*_hiseq_rnaseqv2_gene_exp-raw.tdv.gz)

data/tgca_RNAseq.total.tsv: code/tgca_expression.R $(TGCA_RNAseq_data)
	Rscript code/tgca_expression.R

data/RNAseq.by.gene.tsv: code/expression_analysis.R data/tgca_RNAseq.total.tsv
	Rscript code/expression_analysis.R

#######
####### Calculate Expected N:S
#######

# Load/filter Mutations
data/coding.mutations%rds data/observed.transcripts%rds: code/load_mutations.R
	Rscript code/load_mutations.R

# Choose which transcripts to use
data/final.transcript.list.rds: code/filter_transcripts.R data/observed.transcripts.rds data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz data/raw/mart_export.txt.gz
	Rscript filter_transcripts.R

# Filter variants
data/coding.mutations.filtered%rds data/single.base.coding.substitutions%rds: code/filter_variants.R data/coding.mutations.rds data/raw/mappability_100bp_windows_exons.bed.gz data/final.transcript.list.rds
	Rscript code/filter_variants.R

# Count trimers
data/coding.trimer.counts.rds: code/trimerise_codome.R data/final.transcript.list.rds data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz
	Rscript trimerise_codome.R

# Calculate substitution matrix / mutation profiles
data/motif.probabilities.rds: code/calculate_mutation_profiles.R data/coding.trimer.counts.rds data/single.base.coding.substitutions.rds
	Rscript code/calculate_mutation_profiles.R

# Calculate expected N:S
data/expected_variants_per_transcript.tsv: code/calculate_expected_variants.R data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz data/motif.probabilities.rds data/final.transcript.list.rds
	Rscript code/calculate_expected_variants.R

#######
####### Randomisation Control
#######

# To Do!

#######
####### Expected & Actual Integration
#######

# Join expected and actual N & S, calculate odds ratio and p-values
data/dNdS_byproject%tsv data/dNdS_bysite%tsv data/dNdS_pancancer%tsv: code/calculate_dNdS.R data/expected_variants_per_transcript.tsv data/coding.mutations.filtered.rds
	Rscript code/calculate_dNdS.R

# Plot results
# % are to ensure rule is not run twice for both targets when running make in parallel
results/figures%pdf results/dNdS_by_gene%tsv: code/analyse_dNdS.R data/raw/cancer_gene_census.csv data/RNAseq.by.gene.tsv data/dNdS_byproject.tsv data/dNdS_bysite.tsv data/dNdS_pancancer.tsv data/raw/vogelstein_driver_genes.tdv data/raw/HGNC.tsv
	# Generate histogram and distribution of dNdS
	Rscript code/analyse_dNdS.R
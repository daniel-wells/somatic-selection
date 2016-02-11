.PHONY : results all
results: results/figures.pdf results/dNdS_by_gene.tsv
all: raw_data results

# To print a variable type e.g. "make -n print-TGCA_RNAseq_data"
print-%  : ; @echo $* = $($*)

#################################
####### Download Raw Data #######
#################################

ICGC_project_mutation_files=$(wildcard data/raw/ICGC/simple_somatic_mutation.open.*.tsv.gz)

raw_data: data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz data/raw/mart_export.txt.gz data/raw/simple_somatic_mutation.aggregated.vcf.gz $(ICGC_project_mutation_files) data/raw/cancer_gene_census.csv data/raw/vogelstein_driver_genes.tdv data/raw/ICGC_projects.tsv

data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz:
	# Download ensembl cds nt sequence - referenced as ftp://ftp.ensembl.org/pub/release-75/mysql/homo_sapiens_core_75_37/ at https://docs.icgc.org/methods
	wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.75.cds.all.fa.gz -O data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz
	wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cds/CHECKSUMS -O data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz.CHECKSUMS
	wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cds/README -O data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz.README
	sum data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz
	head -1 data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz.CHECKSUMS
	sha256sum data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz >> data/raw/SHA256SUMS


data/raw/simple_somatic_mutation.aggregated.vcf.gz:
	# Download ICGC aggregated varients
	wget https://dcc.icgc.org/api/v1/download?fn=/release_20/Summary/simple_somatic_mutation.aggregated.vcf.gz -O data/raw/simple_somatic_mutation.aggregated.vcf.gz
	# No hash avaliable
	wget https://dcc.icgc.org/api/v1/download?fn=/release_20/README.txt -O data/raw/simple_somatic_mutation.aggregated.vcf.gz.README
	sha256sum data/raw/simple_somatic_mutation.aggregated.vcf.gz >> data/raw/SHA256SUMS
	sha256sum -c data/raw/SHA256SUMS

data/raw/mart_export.txt.gz:
	# Requires manual download
	http://www.ensembl.org/biomart/martview/a9103936cd932b57917528550c7c9a2b?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.external_gene_name|hsapiens_gene_ensembl.default.feature_page.gene_biotype|hsapiens_gene_ensembl.default.feature_page.transcript_biotype|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id&FILTERS=&VISIBLEPANEL=resultspanel
	sha256sum data/raw/mart_export.txt.gz >> data/raw/SHA256SUMS

data/raw/ICGC_projects.tsv:
	# Manual Download list of projects -O data/raw/ICGC_projects.tsv
	https://dcc.icgc.org/78e07ae9-2c8c-4511-a2b0-77c71475688d

$(ICGC_project_mutation_files): data/raw/ICGC_projects.tsv
	# Download each project (which has somatic mutation data)
	for OUTPUT in $(cut -f 1,7 data/raw/ICGC_projects.tsv | grep -Pv '(\t0|Project)' | cut -f 1)
	do
		echo $OUTPUT
		wget https://dcc.icgc.org/api/v1/download?fn=/release_20/Projects/$OUTPUT/simple_somatic_mutation.open.$OUTPUT.tsv.gz -O data/raw/ICGC/simple_somatic_mutation.open.$OUTPUT.tsv.gz
	done
	# 2.6G in total

data/raw/cancer_gene_census.csv:
	# Download cancer gene census for annotation in analysis
	# Requires Password
	sftp "daniel.wells@well.ox.ac.uk"@sftp-cancer.sanger.ac.uk
	get /files/grch38/cosmic/v75/cancer_gene_census.csv data/raw/cancer_gene_census.csv
	#572 genes

data/raw/vogelstein_driver_genes.tdv:
	# Requires manual download

$(TGCA_RNAseq_data):
	# ?????

#################################
####### Main Analysis #######
#################################

#######
####### RNAseq
#######

TGCA_RNAseq_data=$(wildcard /users/bschuster/sharedscratch/Splicing/GeneExpRaw/*_hiseq_rnaseqv2_gene_exp-raw.tdv.gz)

data/tgca_RNAseq.total.tsv: code/tgca_expression.R $(TGCA_RNAseq_data)
	Rscript code/tgca_expression.R

data/RNAseq.by.gene.tsv: code/expression_analysis.R data/tgca_RNAseq.total.tsv
	Rscript code/expression_analysis.R

#######
####### Expected Variants 
#######

data/single.base.coding.substitutions.rds data/coding.trimer.counts.rds: code/sequence_context.R
	Rscript code/sequence_context.R

data/motif.probabilities.rds: code/sequence_context.R $(ICGC_project_mutation_files)
	Rscript code/sequence_context.R

data/expected_variants_per_transcript.tsv: code/calculate_nonsynon_sites.R data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz data/motif.probabilities.rds
	Rscript code/calculate_nonsynon_sites.R

#######
####### Actual Variants
#######

data/simple_somatic_mutation.aggregated.filtered.vcf: data/raw/simple_somatic_mutation.aggregated.vcf.gz
	# To keep headers (modify grep to append), or modify below by egrep '(^#|missense_variant)'
	# head -13 simple_somatic_mutation.aggregated.vcf > simple_somatic_mutation.aggregated.coding.vcf

	# Remove lines which are irrelevant (e.g. only intron, intergenic)
	zgrep 'missense_variant\|synonymous_variant\|frameshift_variant\|disruptive_inframe_deletion\|disruptive_inframe_insertion\|inframe_deletion\|inframe_insertion\|start_lost\|stop_lost\|exon_variant' data/raw/simple_somatic_mutation.aggregated.vcf.gz > data/simple_somatic_mutation.aggregated.filtered.vcf
	# 2min40s
	# 1,942,374 lines
	#? find data/raw -name "Homo_sapiens.GRCh37.75.cds.all.fa.gz*" -exec touch {} \;

data/observed_variants_by_transcript.tsv: code/flatten_vcf.py data/simple_somatic_mutation.aggregated.filtered.vcf
	# Create table with each variant annotation pair on seperate line
	python code/flatten_vcf.py data/simple_somatic_mutation.aggregated.filtered.vcf > data/observed_variants_by_transcript.tsv
	# 1 min
	# 12,857,317 lines

data/observed_variants_by_transcript.filtered.tsv: data/observed_variants_by_transcript.tsv
	# Subset only synon or non-synon for dnds calculation
	# NB this removes genes which were sequenced but have non measured synon or non-synon variations
	zgrep 'missense_variant\|synonymous_variant\|frameshift_variant\|disruptive_inframe_deletion\|disruptive_inframe_insertion\|inframe_deletion\|inframe_insertion\|start_lost\|stop_lost\|exon_variant' data/observed_variants_by_transcript.tsv > data/observed_variants_by_transcript.filtered.tsv

#######
####### Expected & Actual Integration
#######

data/dNdS_by_transcript.tsv: code/calculate_dNdS.R data/observed_variants_by_transcript.filtered.tsv data/expected_variants_per_transcript.tsv
	Rscript code/calculate_dNdS.R

# % are to ensure rule is not run twice for both targets when running make in parallel
results/figures%pdf results/dNdS_by_gene%tsv: code/analyse_dNdS.R data/raw/cancer_gene_census.csv data/RNAseq.by.gene.tsv data/dNdS_by_transcript.tsv data/raw/vogelstein_driver_genes.tdv data/raw/mart_export.txt.gz
	# Generate histogram and distribution of dNdS
	Rscript code/analyse_dNdS.R
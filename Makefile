.PHONY : results all
results: results/volcano.pdf results/dNdS_by_gene.tsv
all: raw_data results


# To print a variable type e.g. "make -n print-TGCA_RNAseq_data"
print-%  : ; @echo $* = $($*)




#################################
####### Download Raw Data #######
#################################

data/raw/ICGC_projects.tsv:
	# Manual Download list of projects -O data/raw/ICGC_projects.tsv
	echo https://dcc.icgc.org/projects/details
	mkdir data
	mkdir data/raw
	mkdir logs
	mkdir code
	mkdir archive
	mkdir results
	mkdir data/raw/ICGC/

data/url-list.txt: data/raw/ICGC_projects.tsv
	for OUTPUT in $$(cut -f 1,7 data/raw/ICGC_projects.tsv | grep -Pv '(\t0|Project)' | cut -f 1) ; do \
		echo "https://dcc.icgc.org/api/v1/download?fn=/release_20/Projects/$$OUTPUT/simple_somatic_mutation.open.$$OUTPUT.tsv.gz" >> data/url-list.txt \ ; \
	done

URLS=$(shell awk '{printf "%s\n", $$1}' data/url-list.txt)

ICGC_project_mutation_files=$(addprefix data/raw/ICGC/,$(notdir $(URLS)))

# Download simple somatic mutation files
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


data/raw/data/raw/ExAC.r0.3.1.sites.vep.vcf.gz:
	# Download ExAC for SNP frequencies ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/
	# see also /mnt/lustre/data/ExAC/
	curl ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz -o data/raw/ExAC.r0.3.1.sites.vep.vcf.gz
	curl ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz.tbi -o data/raw/ExAC.r0.3.1.sites.vep.vcf.gz.tbi
	curl ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/md5sum.txt -o data/raw/ExAC.md5sum.txt
	sum data/raw/ExAC.r0.3.1.sites.vep.vcf.gz
	head -1 data/raw/ExAC.md5sum.txt
	sha256sum data/raw/ExAC.r0.3.1.sites.vep.vcf.gz >> data/raw/SHA256SUMS

data/raw/cancer_gene_census.csv:
	# Download cancer gene census for annotation in analysis
	# Requires Password
	sftp "daniel.wells@well.ox.ac.uk"@sftp-cancer.sanger.ac.uk
	get /files/grch38/cosmic/v77/cancer_gene_census.csv data/raw/cancer_gene_census.csv
	#595 genes

data/raw/vogelstein_driver_genes.tdv:
	# Requires manual download

$(TGCA_RNAseq_data):
	# Requires manual download

data/raw/HGNC.tsv:
	curl --data 'col=gd_hgnc_id&col=gd_app_sym&col=gd_locus_group&col=gd_pub_eg_id&col=gd_pub_ensembl_id&col=md_eg_id&col=md_ensembl_id&status=Approved&status_opt=2&where=&order_by=gd_app_sym_sort&format=text&limit=&submit=submit' http://www.genenames.org/cgi-bin/download -o data/raw/HGNC.tsv


data/raw/mappability_100bp_windows_exons.bed.gz:
	# Requires manual download
	cp /mnt/lustre/users/bschuster/TCGA/Coverage/mappability_100bp_windows_exons.bed.gz /mnt/lustre/users/dwells/v21/data/raw/mappability_100bp_windows_exons.bed.gz

ssm: $(ICGC_project_mutation_files)

raw_data: data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz $(ICGC_project_mutation_files) data/raw/cancer_gene_census.csv data/raw/vogelstein_driver_genes.tdv data/raw/ICGC_projects.tsv data/raw/mappability_100bp_windows_exons.bed.gz data/raw/HGNC.tsv data/raw/data/raw/ExAC.r0.3.1.sites.vep.vcf.gz $(TGCA_RNAseq_data)


#############################
####### Main Analysis #######
#############################

#######
####### ExAC
#######
data/raw/ExAC.bed: data/raw/ExAC.r0.3.1.sites.vep.vcf.gz data/raw/ExAC.r0.3.1.sites.vep.vcf.gz.tbi
	module load apps/bcftools/1.0/gcc-4.4.7
	bcftools query -f "%CHROM\t%POS\t%ID\t%QUAL\t%REF\t%ALT\t%AF\n" /mnt/lustre/data/ExAC/ExAC.r0.3.1.sites.vep.vcf.gz | perl -ne '@l=split; if (/,/){@a=split /,/, $l[5]; @af=split /,/, $l[6]; for (0..$#a){print join("\t", $l[0], $l[1]-1, $l[1], @l[2..4], $a[$_], $af[$_])."\n"}}else{print join("\t", $l[0], $l[1]-1, $l[1], @l[2..6])."\n"}' > data/ExAC.bed

#######
####### RNAseq
#######

# TGCA_RNAseq_data=$(wildcard /users/bschuster/sharedscratch/Splicing/GeneExpRaw/*_hiseq_rnaseqv2_gene_exp-raw.tdv.gz)
#
# data/tgca_RNAseq.total.tsv: code/tgca_expression.R $(TGCA_RNAseq_data)
# 	Rscript code/tgca_expression.R
#
# data/RNAseq.by.gene.tsv: code/expression_analysis.R data/tgca_RNAseq.total.tsv
# 	Rscript code/expression_analysis.R

#######
####### Calculate Expected N:S
#######

# Load/filter Mutations
data/coding.mutations%rds data/observed.transcripts%rds: code/load_mutations.R
	Rscript code/load_mutations.R

# Choose which transcripts to use
data/final.transcript.list.rds: code/filter_transcripts.R data/observed.transcripts.rds data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz
	Rscript filter_transcripts.R

# Filter variants
data/coding.mutations.filtered%rds data/single.base.coding.substitutions%rds: code/filter_variants.R data/coding.mutations.rds data/ExAC.bed data/raw/mappability_100bp_windows_exons.bed.gz data/final.transcript.list.rds
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
results/volcano%pdf results/dNdS_by_gene%tsv: code/analyse_dNdS.R data/raw/cancer_gene_census.csv data/dNdS_byproject.tsv data/dNdS_bysite.tsv data/dNdS_pancancer.tsv data/raw/vogelstein_driver_genes.tdv data/raw/HGNC.tsv
	# possibly add data/RNAseq.by.gene.tsv
	# Generate histogram and distribution of dNdS
	Rscript code/analyse_dNdS.R
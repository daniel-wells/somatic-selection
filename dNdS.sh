#!/bin/bash
echo "***************************"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "***************************"

PROJECT_DIRECTORY="/users/dwells/sharedscratch/"

############################################################################################################
######## Calculate total Nonsynonymous and Synonymous sites with homemade script from Ensembl fasta ########
############################################################################################################

# Download ensembl cds nt sequence - referenced as ftp://ftp.ensembl.org/pub/release-75/mysql/homo_sapiens_core_75_37/ at https://docs.icgc.org/methods
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.75.cds.all.fa.gz -O data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cds/CHECKSUMS -O data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz.CHECKSUMS
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cds/README -O data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz.README
sum data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz
head -1 data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz.CHECKSUMS
sha256sum data/raw/Homo_sapiens.GRCh37.75.cds.all.fa.gz >> data/raw/SHA256SUMS


http://www.ensembl.org/biomart/martview/a9103936cd932b57917528550c7c9a2b?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.external_gene_name|hsapiens_gene_ensembl.default.feature_page.gene_biotype|hsapiens_gene_ensembl.default.feature_page.transcript_biotype|hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id&FILTERS=&VISIBLEPANEL=resultspanel
sha256sum data/raw/mart_export.txt.gz >> data/raw/SHA256SUMS


Rscript calculate_nonsynon_sites.R
# 15 mins

###############################################################################
####### Calculate (non)synonymous mutations by transcript from ICGC VCF #######
###############################################################################

# Manual Download list of projects -O data/raw/ICGC_projects.tsv
https://dcc.icgc.org/78e07ae9-2c8c-4511-a2b0-77c71475688d

# Download each project (which has somatic mutation data)
for OUTPUT in $(cut -f 1,7 data/raw/ICGC_projects.tsv | grep -Pv '(\t0|Project)' | cut -f 1)
do
	echo $OUTPUT
	wget https://dcc.icgc.org/api/v1/download?fn=/release_20/Projects/$OUTPUT/simple_somatic_mutation.open.$OUTPUT.tsv.gz -O data/raw/ICGC/simple_somatic_mutation.open.$OUTPUT.tsv.gz
done
# 2.6G in total

# Download ICGC aggregated varients
wget https://dcc.icgc.org/api/v1/download?fn=/release_20/Summary/simple_somatic_mutation.aggregated.vcf.gz -O data/raw/simple_somatic_mutation.aggregated.vcf.gz
# No hash avaliable
wget https://dcc.icgc.org/api/v1/download?fn=/release_20/README.txt -O data/raw/simple_somatic_mutation.aggregated.vcf.gz.README
sha256sum data/raw/simple_somatic_mutation.aggregated.vcf.gz >> data/raw/SHA256SUMS
sha256sum -c data/raw/SHA256SUMS

# To keep headers (modify grep to append), or modify below by egrep '(^#|missense_variant)'
# head -13 simple_somatic_mutation.aggregated.vcf > simple_somatic_mutation.aggregated.coding.vcf

# Remove lines which are irrelevant (e.g. only intron, intergenic)
zgrep 'missense_variant\|synonymous_variant\|frameshift_variant\|disruptive_inframe_deletion\|disruptive_inframe_insertion\|inframe_deletion\|inframe_insertion\|start_lost\|stop_lost\|exon_variant' data/raw/simple_somatic_mutation.aggregated.vcf.gz > data/simple_somatic_mutation.aggregated.filtered.vcf
# 2min40s
# 1,942,374 lines

# Create table with each variant annotation pair on seperate line
python code/flatten_vcf.py data/simple_somatic_mutation.aggregated.filtered.vcf > data/observed_variants_by_transcript.tsv
# 1 min
# 12,857,317 lines

# Subset only synon or non-synon for dnds calculation
# NB this removes genes which were sequenced but have non measured synon or non-synon variations
zgrep 'missense_variant\|synonymous_variant\|frameshift_variant\|disruptive_inframe_deletion\|disruptive_inframe_insertion\|inframe_deletion\|inframe_insertion\|start_lost\|stop_lost\|exon_variant' data/observed_variants_by_transcript.tsv > data/observed_variants_by_transcript.filtered.tsv

module load apps/R/3.2.2/gcc-4.4.7+lapack-3.5.0+blas-20110419

# Map counted actual N and S from ICGC VCF to total possible N and S from ensembl cds fasta
Rscript calculate_dNdS.R

# Download cancer gene census for annotation in analysis
sftp "daniel.wells@well.ox.ac.uk"@sftp-cancer.sanger.ac.uk
get /files/grch38/cosmic/v75/cancer_gene_census.csv data/raw/cancer_gene_census.csv
#572 genes

# Generate histogram and distribution of dNdS
Rscript analyse_dNdS.R


echo "***************************"
echo "Finished at: "`date`
echo "***************************"
#!/bin/bash -l

echo "***************************"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "***************************"

#########################
# SGE SUBMISSION SCRIPT #
#########################

# Run in current working directory
#$ -cwd
# Use bash
#$ -S /bin/bash
# Export all variables present at job submission time
#$ -V
# Merge STDOUT and STDERR
#$ -j y
# Set the output file for STDOUT
mkdir -p "$JOB_NAME"
# no SGE_ prefix here for some reason
#@ MONITOR $JOB_NAME/$JOB_NAME.$JOB_ID.$TASK_ID.output
#$ -o $JOB_NAME/$JOB_NAME.$JOB_ID.$TASK_ID.output

# What is the name of your job?
#
#$ -N genome_randomiser

# Set the destination email address
##$ -M EMAILHERE
# Set the conditions under which you wish to receive an email (b:begin, e:end, a:aborted, s:suspended)
##$ -m 

#####################################
# SGE DIRECTIVES - Always set these
#####################################

# Expected RUNTIME
# Format: hours:minutes:seconds, e.g.: 72:0:0
#$ -l h_rt=48:00:00

# Expected HARD MEMORY LIMIT (Per Slot)
#$ -l h_vmem=4G

#=======================
#  Array configuration
#-----------------------
# How many tasks?
#
# For example 1-10 will run 10 jobs, numbered from 1 to 10.
#
#$ -t 1-10000

# How many tasks to schedule at once?
#
# Please note that the total may be limited within your compute environment.
#
# #$ -tc 5


############################
# LOAD ENVIRONMENT MODULES
############################
module load apps/R/3.2.2/gcc-4.4.7+lapack-3.5.0+blas-20110419

mkdir -p "data/randomised_genomes"

###############################
# APPLICATION LAUNCH COMMANDS
###############################

Rscript code/randomise_genome.R $SGE_TASK_ID


echo "***************************"
echo "Finished at: "`date`
echo "***************************"

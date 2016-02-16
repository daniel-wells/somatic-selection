#!/bin/bash
echo "***************************"
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "***************************"

PROJECT_DIRECTORY="/mnt/lustre/users/dwells"

cd $PROJECT_DIRECTORY

make

module load apps/R/3.2.2/gcc-4.4.7+lapack-3.5.0+blas-20110419

echo "***************************"
echo "Finished at: "`date`
echo "***************************"
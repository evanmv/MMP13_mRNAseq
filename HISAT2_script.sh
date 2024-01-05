#!/bin/sh
#SBATCH --job-name=HISAT2
#SBATCH --account=ec12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=10G

#setup
set -o errexit
set -o nounset

module --quiet purge

module load HISAT2/2.2.1-gompi-2022a
module load SAMtools/1.17-GCC-12.2.0

#Argument
ref=$1
fq1=$2
fq2=$3
samout=$4

#New directory 
mkdir -p $SUBMITDIR/alignment_output

#Run Bowtie2
hisat2 -x $SUBMITDIR/grch38/$ref -1 $SUBMITDIR/Raw_reads/$fq1 -2 $SUBMITDIR/Raw_reads/$fq2 \
-S $SUBMITDIR/alignment_output/$samout --summary-file hisat_summary


#!/bin/sh
#SBATCH --job-name=TrimGalore
#SBATCH --account=ec12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:30:00
#SBATCH --mem=1G

#setup
set -o errexit
set -o nounset

module --quiet purge

module load TrimGalore/0.6.10-GCCcore-11.3.0

#Argument
fq1=$1
fq2=$2

#Run TrimGalore
trim_galore --paired $SUBMITDIR/../Raw_reads/$fq1 $SUBMITDIR/../Raw_reads/$fq2 

#Finish message
echo "Done"

#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=QC_s2

#Note: I found this script from here:
#https://www.biostars.org/p/317385/

cd $SLURM_SUBMIT_DIR

#not necessary to submit as batch script

#location of raw reads to be merged
rawreads="/work/bs66/dasanthera_novaseq/raw_reads"
cd $rawreads


#for loop to find matching file names and merge
#change syntax of shared file names if necessary, though should work as-is

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)

    do echo "Merging R1"

cat "$i"_L00*_R1_001.fastq.gz > "$i"_merged_L001_R1_001.fastq.gz

       echo "Merging R2"

cat "$i"_L00*_R2_001.fastq.gz > "$i"_merged_L001_R2_001.fastq.gz

done;
#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=pullCDS
#SBATCH --output=slurm-generate-CDS-fasta_%j.out


#This is an arrayed batch script, requiring special syntax to submit. Submit with:
#sbatch --array [0-17] VCFcall_4b.ARRAY_generate_CDS_fastas.sh


cd $SLURM_SUBMIT_DIR

#source appropriate environments to enable use of conda installs through batch submission
source /home/bs66/.bashrc
source /home/bs66/.bash_profile


#activate conda environment with gffread installed
conda activate mapping_etc



#path to CDS.bed file from the annotation, and desired outfile
bedfile="/work/bs66/project_compare_genomes/annot_Pdavidsonii_1mb.gffread.genes.bed"
outdir="/work/bs66/dasanthera_novaseq/consensus_alignments/individual_CDS_fastas"


#path to directory with genomes (fasta format) for each sample
genomedir="/work/bs66/dasanthera_novaseq/consensus_alignments/individual_fullgenome_fastas"


#set up variables array
filelist=($genomedir/*.fa)
infile="${filelist[$SLURM_ARRAY_TASK_ID]}"
samplename="${infile##*/}"
samplename="${samplename/consensus_fullgenome_/}"



#gffread command to extract CDS for all mapped individuals
gffread -x $outdir/CDS_$samplename -C -M -K -Y -E --sort-alpha -g $infile $bedfile



#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=generate_genetrees

#this is a masterscript to generate sliding window trees 
#This is an array script. Array is 0 to n-1, where n is number of scaffolds in scaflist


cd $SLURM_SUBMIT_DIR


#source bash profile and activate conda environment with python3 as base python env
source /home/bs66/.bashrc
source /home/bs66/.bash_profile


#paths to important scripts and directories
pythonscript="/work/bs66/dasanthera_novaseq/analysis/TREES.create_fasta_window_alignments.py"
scaffold_dir="/work/bs66/dasanthera_novaseq/consensus_alignments/scaffold_fullgenome_fastas"
outdir="/work/bs66/dasanthera_novaseq/analysis/genetree_infiles"


#list of scaffold names
scaflist=("1085" "1086" "1087" "2151" "2446" "2531" "2532" "2533" "2684" "2685" "2686" "2687")
inscaf="${scaflist[$SLURM_ARRAY_TASK_ID]}"

#Use the python script in the directory to filter accordingly
#syntax: python3 scriptname input_fasta windowsize prefix
python3 $pythonscript $scaffold_dir/allseqs_consensus_scaffold_$inscaf.fa 10000 0.75 $outdir/scaf_$inscaf/scaffold_$inscaf

#Alternative code, for submission without an ARRAY:
#run for loop for python script to generate gene tree files for each scaffold
#for i in "${scaflist[@]}"
#do
#	python3 $pythonscript $scaffold_dir/allseqs_consensus_scaffold_$i.fa 10000 0.75 $outdir/scaf_$i/scaffold_$i
#done


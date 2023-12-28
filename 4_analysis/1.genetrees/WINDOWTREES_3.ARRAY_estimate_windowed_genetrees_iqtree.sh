#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p wessinger-48core
#SBATCH --job-name=iqtree_genetrees

##NOTE: this is an arrayed batch script. The array is 0-n, where n is the number of scaffolds-1. For P. davidsonii genome, this is sbatch --array [0-11]


cd $SLURM_SUBMIT_DIR

#load conda env. Need iqtree installed.
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc


#path to directory holding subdirectories -- one for each scaffold
aligndir="/work/bs66/dasanthera_novaseq/analysis/genetree_infiles"

#set up variable to infer trees in each scaffold directory
indirs=($aligndir/*/)
inscaf="${indirs[$SLURM_ARRAY_TASK_ID]}"


#change to the input directory
cd $inscaf

#for loop to go through each alignment, estimate subs. model, and est. tree
#the ARRAY number changes the scaffold for which this occurs
#"rooting" at P. montanus, but this is still an unrooted tree
for i in *.fa
do
	iqtree -s $i \
	 --seqtype DNA \
	 -m MFP \
	 -o mon_61-7_S440 \
	 -T 4 \
	 --prefix outfile_$i
done

#list outdir file and move files
outdir=${inscaf/genetree_infiles/genetree_outfiles}
mv outfile* $outdir



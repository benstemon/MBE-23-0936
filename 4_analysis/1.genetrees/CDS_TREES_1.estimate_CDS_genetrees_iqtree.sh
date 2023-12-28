#!/bin/sh

#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p wessinger-48core
#SBATCH --job-name=iqtree_CDS_genetrees
#SBATCH --output=slurm-iqtree_CDS_genetrees_%j_%a.out


cd $SLURM_SUBMIT_DIR

#load conda env. Need iqtree installed.
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc


#path to directory holding subdirectories -- one for each scaffold
aligndir="/work/bs66/dasanthera_novaseq/analysis/CDS_genetree_infiles"
outdir="/work/bs66/dasanthera_novaseq/analysis/CDS_genetree_outfiles"


#for loop to go through each alignment, estimate subs. model, and est. tree
cd $aligndir
for i in *.fa
do
	iqtree -s $i \
	 --seqtype DNA \
	 -m MFP \
	 -T 4 \
	 --prefix outfile_$i
done


#list outdir file and move files
mv outfile* $outdir



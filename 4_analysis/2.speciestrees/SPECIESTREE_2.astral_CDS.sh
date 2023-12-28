#!/bin/sh

#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p wessinger-48core
#SBATCH --job-name=astral_tree
#SBATCH --output=slurm-astral-CDStree_%j.out


cd $SLURM_SUBMIT_DIR

#path to ASTRAL
astral="/work/bs66/software/ASTRAL/Astral/astral.5.7.8.jar"

#path to window trees, catted into combined treefile
treefile="/work/bs66/dasanthera_novaseq/analysis/treemetrics/combined_CDStrees.tre"

#path to output file
outpath="/work/bs66/dasanthera_novaseq/analysis/astral_trees"


#run astral. Using single ML estimates of gene trees rather than bootstrap resampling.
#the branch annotations are for full annotations
java -jar $astral --input $treefile --branch-annotate 2 --output $outpath/astral_CDS_annotations.tre


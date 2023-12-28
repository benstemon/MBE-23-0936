#!/bin/sh

#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p wessinger-48core
#SBATCH --job-name=concat_tree
#SBATCH --output=iqtree_speciestree_%j.out


cd $SLURM_SUBMIT_DIR

#load conda env. Need numpy, matplotlib, biopython, installed.
source /home/bs66/.bashrc
source /home/bs66/.bash_profile
conda activate mapping_etc

#path to scaffold fasta alignments
scaffolds="/work/bs66/dasanthera_novaseq/consensus_alignments/scaffold_fullgenome_fastas"


#path to output file
outpath="/work/bs66/dasanthera_novaseq/analysis/iqtree_speciestree"



#infer species tree in IQtree
iqtree -p $scaffolds \
 -m GTR+I+R \
 -o mon_61-7_S440 \
 -T 32 \
 -B 1000 \
 --prefix $outpath/concat_speciestree_scaffolds


## Analysis
This README file contains code and links to scripts/directories that contain all the information necessary to conduct the analyses in the manuscript. Because certain analyses require input from upstream analysis (e.g., constructing the species tree requires sliding window trees or gene trees as input), analyses should be conducted sequentially as they are described in this README.

The main sections are as follows:
1. [Gene Trees (gene, sliding window)](README.md#1-gene-trees-gene-sliding-window)
2. [Species trees (from ASTRAL with sliding windows and CDS, and a concatenated ML approach)](README.md#2-species-trees-from-astral-with-sliding-windows-and-cds-and-a-concatenated-ml-approach)
3. [Tree metrics](README.md#3-tree-metrics)
4. [D-statistics](README.md#4-d-statistics)
5. [TWISST](README.md#5-twisst)
6. [PIXY](README.md#6-pixy)
7. [f3p5ph](README.md#7-f3p5ph)



### 1. [Gene Trees (gene, sliding window)](1.genetrees/)

#### 1a. Windowed gene trees

##### Prepare fasta files for windowed gene tree inference
The fasta files generated from consensus have newlines every 50bp. Additionally, they are sorted by species. We want to have a fasta file for each scaffold, with the sequence for that scaffold from each species. There is a script, [`WINDOWTREES_1.rearrange_consensus_sequences.sh`](1.genetrees/WINDOWTREES_1.rearrange_consensus_sequences.sh), which does this. Note that this functions on Linux OS, and will need modified slightly if using on MacOS (read script for more details). Script will need edited to match variable naming scheme. Use:

`bash WINDOWTREES_1.rearrange_consensus_sequences.sh`


##### Generate windowed gene tree input files
This is done in two parts, with three files:
1. Make output directories for each scaffold. 
2. Generate gene trees with [`WINDOWTREES_2.ARRAY_generate_windowed_genetree_infiles_masterscript.sh`](1.genetrees/WINDOWTREES_2.ARRAY_generate_windowed_genetree_infiles_masterscript.sh), which uses the custom python script [`WINDOWTREES.create_fasta_window_alignments.py`](1.genetrees/WINDOWTREES.create_fasta_window_alignments.py). It is currently set up as an array script. Parameters are defined in the script, including:
	* Input fasta file
	* Window size
	* Missing data threshold (only generates windows in which all species pass missing data threshold)
	* prefix to append to outfiles (can function as outdir)


##### Estimate windowed gene trees
This pipeline uses [IQtree](http://www.iqtree.org/) for gene tree inference. Again, this is done in two parts: 
1. Make output directories for gene trees estimated along each scaffold.
2. Estimate gene trees (in array batch submission), setting outgroup (here it is P. montanus), and specifying substitution model inference and out prefix.
* See [`WINDOWTREES_3.ARRAY_estimate_windowed_genetrees_iqtree.sh`](1.genetrees/WINDOWTREES_3.ARRAY_estimate_genetrees_iqtree.sh)



#### 1b. CDS gene trees

##### Prepare CDS fasta files for gene tree inference
The fasta files generated for CDS regions are initially sorted by sample. To estimate CDS gene trees, we need a fasta for each CDS, with each sample's sequence included. First we will want to reformat our fasta files, so they are single line fasta rather than interleaved. Use this code on the CDS fastas.
```shell

for i in *.fa;
do
    awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < $i > "${i/.fa/.fixed.fa}"
    rm $i
done

```

Next, we will want to generate the infiles for gene tree inference. To generate these, see [`CDS_TREES.concat_CDS_fastas.py`](1.genetrees/CDS_TREES.concat_CDS_fastas.py). This python script takes input fasta, output directory, and missing data threshold, and appends to an output fasta for each CDS, naming the output after the scaffold and region the CDS corresponds to. A simple shell for loop can be run with the python script to add all samples to the output. The missing data threshold here filters for individuals, rather than windows (i.e., if a sample has more missing data than desired, that individual is not added to the output fasta, rather than the output fasta not being generated).

```shell
for i in individual_CDS_fastas/CDS_*.fa;
do
    python3 CDS_TREES.concat_CDS_fastas.py -i $i -o /work/bs66/dasanthera_novaseq/analysis/CDS_genetree_infiles -m 0.5
done
```

Here I generate aligned fasta for each CDS, excluding individuals with > 50% missing data. The script is fast, and possible to use on an interactive node with this number of samples (18). But for many samples or many CDS regions, it would be wise to submit as a batch/array script.


I also generated individual fasta files for each CDS for the reference genome. This will be useful for later inferences when looking at function of specific genomic regions.

```
#give the 1mb genome and desired output directory

refgenome="/work/bs66/project_compare_genomes/annot_Pdavidsonii_genome.1mb.fasta"
bedfile="/work/bs66/project_compare_genomes/annot_Pdavidsonii_1mb.gffread.genes.bed"
outdir="/work/bs66/project_compare_genomes/davidsonii_CDS"


gffread -x $outdir/concat_CDS_davidsonii_refgenome.fa -C -M -K -Y -E --sort-alpha -g $refgenome $bedfile


#change to single line per fasta
for i in *.fa;
do
    awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < $i > "${i/.fa/.fixed.fa}"
    rm $i
done


#use python script to generate individual fastas for each CDS
python3 CDS_TREES.concat_CDS_fastas.py -i concat_CDS_davidsonii_refgenome.fixed.fa -o $outdir

#remove files no longer needed
#rm CDS_TREES.concat_CDS_fastas.py
#rm concat_CDS_davidsonii_refgenome.fa
```



##### Estimate CDS gene trees
Still using IQtree for gene tree inference here. This used to be done on a scaffold-by-scaffold basis, but because that information has been lost in the header, instead we are just submitting all files for inference. They could be split up into sections to speed this part up.
* To infer gene trees, see [`CDS_TREES_1.estimate_CDS_genetrees_iqtree.sh`](1.genetrees/CDS_TREES_1.estimate_CDS_genetrees_iqtree.sh)





### 2. [Species trees (from ASTRAL with sliding windows and CDS, and a concatenated ML approach)](2.speciestrees/)

#### 2a. ASTRAL species tree: windowed analysis

First, cat trees into one large file. Change to directory with subdirectories, each containing gene trees estimated for each scaffold.

```shell
cd /work/bs66/dasanthera_novaseq/analysis/genetree_outfiles
for i in scaf_*;
do
    cat $i/*.treefile >> combined_10kbwindowtrees.tre
done
```

Also, keep a log of which trees were catted to the treefile, and in what order. This will come in handy later.

```shell
for i in scaf_*;
do
    echo $i/*.treefile | tr " " "\n" >> tmpout.txt
done

cat --number tmpout.txt > numbered_10kbtreepaths.txt
rm tmpout.txt
```

I put both of these outputs in a new directory, `treemetrics`. They will be used for several downstream analyses. Next, use the combined windowtrees file to estimate the species tree in ASTRAL.

* See [`SPECIESTREE_1.astral_windowtrees.sh`](2.speciestrees/SPECIESTREE_1.astral_windowtrees.sh)



#### 2b. ASTRAL species tree: CDS analysis

Perform the same commands as described for the windowed analysis; cat trees to a combined tree file and keep a log of the order in which this was done. 

```shell
cd /work/bs66/dasanthera_novaseq/analysis/CDS_genetree_outfiles
for i in *.treefile;
do
    cat $i >> combined_CDStrees.tre
done

for i in *.treefile;
do
    echo $i | tr " " "\n" >> tmpout.txt
done

cat --number tmpout.txt > numbered_CDStreepaths.txt
rm tmpout.txt

mv combined_CDStrees.tre ../treemetrics
mv numbered_CDStreepaths.txt ../treemetrics
```

Then, estimate the species tree in ASTRAL.
* See [`SPECIESTREE_2.astral_CDS.sh`](2.speciestrees/SPECIESTREE_2.astral_CDS.sh)



#### 2c. Concatenated ML species tree
The concatenated ML species tree is estimated in IQtree. I specified the GTR+I+R model, and estimated rates for each scaffold, performing 1000 ultra-fast bootstrap replicates, and specifying _P. montanus_ as the outgroup.
* See [`SPECIESTREE_3.IQtree_concat.sh`](2.speciestrees/SPECIESTREE_3.IQtree_concat.sh)





### 3. [Tree metrics](3.treemetrics/)
These are scripts to calculate RF distance and concordance factors, assuming gene trees, window trees, and species trees have been generated.


#### 3a. Calculate gene content in sliding windows (prep)
* see [`TREEMETRICS_1.calculate_percentage_CDS_slidingwindows.sh`](3.treemetrics/TREEMETRICS_1.calculate_percentage_CDS_slidingwindows.sh). Using this script will require a .txt file with the names and sizes of each chromosome (see [`genomesize_scaffolds_davidsonii_1mb.txt`](3.treemetrics/genomesize_scaffolds_davidsonii_1mb.txt)) and a .bed file with the genomic coordinates of CDS. This could also be pulled from the fasta.fai file.
* To enable use of the script as intended, it must be made executable with `chmod +x TREEMETRICS_1.calculate_percentage_CDS_slidingwindows.sh` and run with `./TREEMETRICS_1.calculate_percentage_CDS_slidingwindows.sh`



#### 3b. Robinson-Foulds distance
Estimate RF distance of 10kb sliding window trees compared to the species tree. This was plotted along with genic content to visualize whether there is a relationship between RF distance and genic content along scaffolds. For the script to calculate RF distance and plot results, see [`TREEMETRICS_2.calculate-RF_plot-genecontent.R`](3.treemetrics/TREEMETRICS_2.calculate-RF_plot-genecontent.R)



#### 3c. Concordance factors
Estimated site concordance factors and gene concordance factors in IQtree. See [`TREEMETRICS_3.calculate_concordance_factors.sh`](3.treemetrics/TREEMETRICS_3.calculate_concordance_factors.sh). Prior to running for CDS, needed to filter out sites with fewer than 3 samples, using provided python script [`TREEMETRICS_remove_fewsequence_CDS.py`](3.treemetrics/TREEMETRICS_remove_fewsequence_CDS.py)

From this, used the .log file to match locus IDs (10kb trees) to their site concordance factors score. First, I copied the output from the IQtree.log file containing identification into a new .txt file, then used `grep -v '^WARNING' 10kbsource_10kbref.log > 10kbsource_cf_site_ids.txt` to generate an output text file with only this information.

* This information will be plotted later, in the same script as plots for TWISST internal branch plots. (To jump ahead, see [`plot_twisst_IB-test_and_CFs_IB.R`](twisst_internal_branches/plot_twisst_IB-test_and_CFs_IB.R)) 





### 4. [D-statistics](4.Dstats/)
We used Dsuite to calculate all introgression statistics, including f-branch tests and fdM. These were performed in three main groups of analyses.


#### 4a. Dstats-fbranch (whole-genome signatures of introgression)
First, perform Dsuite Dtrios to generate D statistics for all possible triplets of taxa. These tests do not include P. lyallii and use the Astral 10kb topology as the input species tree. Samples here were treated as individuals.
* See [`DSTATS_1.ARRAY_dsuite_dtrios-fbranch.sh`](4.Dstats/DSTATS_1.ARRAY_dsuite_dtrios-fbranch.sh) for this first batch script, as well as [`intree_dtrios-fbranch.tre`](4.Dstats/intree_dtrios-fbranch.tre) for the input species tree and [`popset_dtrios-fbranch.txt`](4.Dstats/popset_dtrios-fbranch.txt) for the population identification set.

Next, combine Dtrios output for each scaffold into a genome-wide analysis. This will also produce an f-branch plot for the whole genome.
* See [`DSTATS_2.combineDtrios-fbranch.sh`](Dstats_fbranch/DSTATS_2.combineDtrios-fbranch.sh)



#### 4b. Dstats-fullspecies-sliding (signatures of introgression across genomic regions)
We used Dsuite Dinvestigate to calculate introgression metrics in sliding windows across the genome. These scripts estimate D, fd, fdM, and Df in overlapping windows of 10kb SNPs, sliding every 2500 SNPs. Rather than as individuals, samples were assigned to species identity.

We calculated these introgression statistics for every possible rooted triplet of taxa (ten total rooted triplets), specifying relationships as inferred by the species tree.
* For the full sample analysis, see [`DSTATS_3.sliding_10kb-fullspecies.sh`](4.Dstats/DSTATS_3.sliding_10kb-fullspecies.sh) for the shell script,  [`popset-fullspecies.txt`](4.Dstats/popset-fullspecies.txt) for the population set, and [`trioset-fullspecies.txt`](4.Dstats/trioset-fullspecies.txt) for specifying each of the rooted triplets.

We were also interested in the potential relationship that introgression metrics have with gene density. To explore this relationship, we calculated the proportion of coding sites for each of the 10kb SNP windows using a custom script and bedtools. See script [`DSTATS_4.generate_dstat_genic_plotfiles.sh`](4.Dstats/DSTATS_4.generate_dstat_genic_plotfiles.sh) for the script. It is used as follows:
```bash
# Input for the script includes:
# -d : the d-statistics outfile produced from Dinvestigate
# -c : a .bed file with CDS coordinates
# -o : the name for the tab-separated output file


chmod +x generate_dstat_genic_plotfiles.sh

./DSTATS_4.generate_dstat_genic_plotfiles.sh -d dstat_file.txt -c CDS_regions.bed -o outfile.bed

# the script will generate a "plotting" file which can be used to compare Dinvestigate output with genic fraction.
```

* Finally, see [`DSTATS_5.plot_fullspecies_sliding_genic.R`](4.Dstats/DSTATS_5.plot_fullspecies_sliding_genic.R) for plotting options.



#### 4c. Dstats-new-rup (signatures of introgression specific to the newberryi-cardwellii-rupicola branch)

We used Dsuite to calculate introgression metrics for the *P. newberryi* - *P. cardwellii* - *P. rupicola* branch of the species tree. This was done with and without two samples from Castle Lake (one each of *P. newberryi* and *P. rupicola*) to see whether (1) there was a prevailing significant signature of introgression between new-rup, (2) whether removing those samples changed the result, and (3) to assess whether these values were elevated specifically around the f3'5'h locus.

* The main bash script [`DSTATS_5.new-rup.sh`](4.Dstats/DSTATS_6.new-rup.sh) runs all the tests, and uses two sets of population files (one each for the [`full scheme`](4.Dstats/popset_new-rup-full.txt) and [`no Castle lake scheme`](4.Dstats/popset_new-rup-nohybrid.txt)) and one [`shared trioset`](4.Dstats/trioset_new-rup.txt).

* See [`DSTATS_7.plot_new-rup_fdm.R`](4.Dstats/DSTATS_7.plot_new-rup_fdm.R) for plotting options and analysis.





### 5. [TWISST](5.TWISST/)
We ran TWISST using the 10kb sliding window trees as input.


#### TWISST internal branch test
We considered all possible tree topologies consistent wth the three main internal branches of the species tree: 1: new+car+rup, 2: new+car, and 3: dav+fru.

* [`TWISST_1.run_twisst_IB-test.sh`](5.TWISST/TWISST_1.run_twisst_IB-test.sh) is the main script to run TWISST. Note that this uses a file named `combined_10kbwindowtrees.tre`, which was generated all the way back during the ATRAL species tree step. It also uses a groupsfile, which connects tips on the window trees with species ID: [`groupsfile_IB-test.txt`](5.TWISST/groupsfile_IB-test.txt).

* This output (topology weights) is then plotted for each 5 taxon tree consistent with each of the three internal branches (along with concordance factors) with [`TWISST_2.plot_twisst_IB-test_and_CFs_IB.R`](5.TWISST/TWISST_2.plot_twisst_IB-test_and_CFs_IB.R).





### 6. [PIXY](6.pixy/)
We ran pixy to calculate dxy, Fst, and pi in 10kb and 50kb non-overlapping windows.

#### 6a. Generate .bed file in matching increments.

The goal of these analyses are to visualize how these diversity and differentiation metrics vary across the genome, and whether there are shared patterns across species and/or with particular genomic features. One key metric is gene density. You should already have a .bed file with genic fraction in sliding 10kb windows, but if not (or if you need to make it for a different window size), see the [TREEMETRICS](3.treemetrics/) section, and specifically the script [`TREEMETRICS_1.calculate_percentage_CDS_slidingwindows.sh`](3.treemetrics/TREEMETRICS_1.calculate_percentage_CDS_slidingwindows.sh).



#### 6b. Run and plot pixy results across the genome

* To run the main batch script, see [`PIXY_1.10kb-50kb_allpops-fullspecies.sh`](6.pixy/PIXY_1.10kb-50kb_allpops-fullspecies.sh). You will also need the population file, which is provided here: [`popfile_allpops-fullspecies.txt`](6.pixy/popfile_allpops-fullspecies.txt). These results are plotted with [`PIXY_2.plot_pixy-fullspecies.R`](6.pixy/PIXY_2.plot_pixy-fullspecies.R)



#### 6c. Find dxy outliers between species that form hybrid zones

* Here we use an R script [`PIXY_3.dxy_outliers_hybridzone.R`](6.pixy/PIXY_3.dxy_outliers_hybridzone.R) to identify dxy outliers that are in the same genomic windows as another species. Specifically, we assess which regions of the genome have especially high levels of dxy in both *P. davidsonii* x *P. newberryi* and *P. davidsonii* x *P. rupicola* contrasts.





### 7. [f3p5ph](7.f3p5ph)
We ran a series of analyses to re-call variants surrounding the f3'5'h locus, generate a gene tree for f4'5'h, and compare protein and nucleotide sequences for each sample in the study.


#### Re-call variants around the f3'5'h locus and estimate a gene tree
This script re-calls variants and additionally calls indels. It then uses IQtree to estimate a gene tree, after reconstructing CDS from the three f3'5'h exons. This produces all of the files needed to manually compare samples. The script can be found here: [`F3P5PH_1.call_variants-filter-genetree.sh`](7.f3p5ph/F3P5PH_1.call_variants-filter-genetree.sh)


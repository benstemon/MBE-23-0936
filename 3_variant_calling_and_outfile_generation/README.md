## Genotyping and variant calling


### Prepare .bed file with annotated regions of interest
In this case, the .bed file produces genic regions, which we will use downstream to generate consensus sequences for genes of interest

* See [`VCFcall_1.gffread_convert_gff-to-bed.sh`](VCFcall_1.gffread_convert_gff-to-bed.sh)
* Requires installation of [gffread](https://github.com/gpertea/gffread)



### Call genotypes to produce an "allsites" .vcf
The VCF includes variants in addition to invariant sites.

* see [`VCFcall_2.call_genotypes_mpileup.sh`](VCFcall_2.call_genotypes_mpileup.sh)
* Filters on base quality <20 and mapping quality <20, and calls invariant sites in such a way that you can filter on read depth and other parameters for invariant sites if you choose.
* uses bcftools >= v1.15.1



### Filter vcfs (goal: generate consensus sequence)
These filters are for the production of consensus sequences for use in a phylogenetic context. The filters used here are therefore appropriate for this purpose, but are likely not suitable for generating VCFs in a population genomic context.

Invariant and variant sites are filtered differently, because we want stricter filters on sites which are putatively variable. This means there are three main steps: (1) filtering invariant sites, (2) filtering variant sites, and (3) recombining post-filtered data into a single VCF. For this reason there are three main steps conducted in a single script, [`VCFcall_3.filter_vcf_consensus.sh`](VCFcall_3.filter_vcf_consensus.sh).

* Step 1 filters invariant sites. Currently filters are placed on invariant sites, but one could reasonably apply depth and/or missingness filters.

* Step 2 filters variant sites. Filters implemented include genotype quality <20, minimum mean depth 3x, maximum mean depth 60x, minimum depth (for individual genotype call) 2x, maximum missing data allowed = 50%. In addition, genotypes at sites with low quality variants (QUAL <20) are changed to match the reference allele.

* Step 3 combines the two separate filtered VCFs into one large consensus VCF. From here one can move on to the production of consensus sequence alignments (fasta files) or perform certain analyses straight from the consensus VCF.



### Generate consensus whole genome and CDS fastas
First we generate consensus genome sequences. This script uses the reference genome in conjunction with the consensus VCF file to generate consensus sequences for each sample. The command `samtools faidx $referencegenome` will need to be run prior to the use of this script. This is an ARRAY script, meaning that to run, you will need to implement some form of `sbatch --array [0-n] scriptname.sh`, where n is the number of samples -1.
* See [`VCFcall_4a.ARRAY_generate_consensus_fullgenome.sh`](VCFcall_4a.ARRAY_generate_consensus_fullgenome.sh)
<br />

Next, we will make use of the file generated in the first script to pull CDS from the fastas. This only functions correctly if consensus fastas are aligned perfectly to the reference fasta (NO INDELS); it just pulls the same CDS regions that were annotated in the reference, based on position.

First, .bed file needs filtered to contain only scaffolds that are actually in the genome. I made a file named [`scaflist.txt`](scaflist.txt) which contained the names of each scaffold I wanted to retain, one name per line. Then, I filtered with grep:
`grep -f scaflist.txt annot_Pdavidsonii_gffread.genes.bed > annot_Pdavidsonii_1mb.gffread.genes.bed`

Then, use this newly created filtered .bed file to generate CDS sequences for each sample.
* See [`VCFcall_4b.ARRAY_generate_CDS_fastas.sh`](VCFcall_4b.ARRAY_generate_CDS_fastas.sh)


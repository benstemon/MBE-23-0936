# Ecological diversification in an adaptive radiation of plants: the role of de novo mutation and introgression

This data set corresponds to Stone and Wessinger 2023: "Ecological diversification in an adaptive radiation of plants: the role of de novo mutation and introgression". [`DOI: 10.1093/molbev/msae007`]([https://www.biorxiv.org/content/10.1101/2023.11.01.565185v1.full.pdf+html](https://academic.oup.com/mbe/article/41/1/msae007/7564791))

The link to supplemental data on FigShare is here: [`Dataset DOI: 10.6084/m9.figshare.24480499`](https://figshare.com/articles/dataset/Data_from_Ecological_diversification_in_an_adaptive_radiation_of_plants_the_role_of_de_novo_mutation_and_introgression/24480499)

The raw reads associated with this study are stored in the SRA repository under project number PRJNA1057825: [`https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1057825`](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1057825)

Note that this study uses the *Penstemon davidsonii* genome as a reference for mapping. Information about that genome can be found in [`Ostevik et al. 2023`](https://www.biorxiv.org/content/10.1101/2023.08.31.555769v1), and the genome and associated materials are stored in the SRA here: [`https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1010203`](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1010203)

## Directory setup
### [Quality control](1_QC/)
Scripts for assessing read quality and filtering raw reads.


### [Mapping and Filtering](2_mapping_and_filtering/)
Scripts for mapping reads to the *P. davidsonii* reference genome, quality filtering, and generating summary statistics.


### [Variant calling, .vcf filtering, and generation of analysis files (.vcfs, multiple sequence alignments, etc.)](3_variant_calling_and_outfile_generation)
Scripts for variant calling, filtering .vcf files, and generating files for downstream analysis.


### [Analysis](4_analysis/)
Contains scripts and code to conduct all further analysis on files generated thus far.


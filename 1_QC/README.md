## NovaSeq data for 18 accessions of *Penstemon* subgenus *Dasanthera*

### Note: this pipeline requires the use of conda installs of various softwares. Ensure the proper conda environments are installed prior to use.

### Quality control pipeline
1. Run QC of raw sequencing reads (fastqc) and summarize (multiqc)
2. Merge Illumina lanes by forward and reverse reads
3. 
    a. Trim adapters, quality filter, and enable base correction in overlapped regions (fastp)
    b. Run QC on trimmed, filtered data (fastqc)
    c. Summarize results (multiqc)


#### 1. QC on raw reads with fastqc, summarized with multiqc
* see [`QC_s1_rawqc-summarize.sh`](QC_s1_rawqc-summarize.sh)


#### 2. Merge data across lanes
* See [`QC_s2_merge-illumina-lanes.sh`](QC_s2_merge-illumina-lanes.sh)
* After using this script, move the merged reads to a new directory


#### 3a-c. Trimming and quality filtering (fastp), check QC on filtered reads (fastqc) and summarize (multiqc)
* See [`QC_s3_fastp-filterqc-summarize.sh`](QC_s3_fastp-filterqc-summarize.sh)
* Default options for quality filtering
* Base correction for overlapping reads enabled
* poly-x trimming on 3' ends enabled
* Limit read length to 30 bp
* Enable auto-detection of adapters

### Mapping, deduplication, overlap clipping, and summary stats

* Reads are mapped to the *P. davidsonii* reference genome with bwa-mem
* Duplicates marked and removed with samtools markdup
* Reads with low mapping quality (q<20)
* Overlapping paired end reads clipped with bamutil clipOverlap
* Some basic summary statistics generated with samtools coverage and samtools stats

All of the main commands are piped to avoid the creation of many intermediate files. See [`mapping_pipeline_array.sh`](mapping_pipeline_array.sh)
  * Uses bwa, samtools v1.15.1, and bamUtil v1.0.15

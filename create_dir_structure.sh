#! /bin/bash
# Create directory structure, move files
# and delete unwanted file after the pipeline is run

mkdir fastqc_results bam qualimap_results read_counts raw_fastq trimmed_fastq
mv *html *zip fastqc_results
mv *trim.sorted.bam *.bai bam/
mv *.counts *counts.summary read_counts/
mv *trim.fastq trimmed_fastq/
mv *sorted_stats qualimap_results
mv *.fastq raw_fastq
rm *.sam *.unpaired

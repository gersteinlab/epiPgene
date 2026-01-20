#!/bin/bash

## Download files from ENCODE (DNAme)
sh bin/download.metadata.sh "https://www.encodeproject.org/metadata/?status=released&internal_tags=ENTEx&assay_title=DNAme+array&type=Experiment&files.processed=true" metadata.tsv

## list of bed_bed3+ files
grep bed_bed3+ metadata.tsv | grep -v hg19 |  cut -f1,7-8,11,35,46 | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$6}'  > bed.files.tsv

## download bed_bed3+ files
cut -f1 bed.files.tsv | while read file; do echo -e "wget -P bed.files https://www.encodeproject.org/files/"$file"/@@download/"$file".bed.gz; md5sum bed.files/"$file".bed.gz >> bed.gz.md5sum"; done > bed.download.job.list.txt

module load dSQ
dsq --job-file bed.download.job.list.txt --mem-per-cpu 8g -t 24:00:00 --mail-type ALL --partition pi_gerstein,day,scavenge
sbatch dsq-bed.download.job.list-2023-07-28.sh

## check md5sum
bin/join.py -a <(sed 's/bed.files\///g;s/.bigWig//g' bed.gz.md5sum | awk 'BEGIN{OFS="\t"}{print $2, $1}') -b bed.files.tsv | awk '$(NF-1)!=$NF'

## prepare bed exp list
module load R
Rscript bin/bed.exp.R

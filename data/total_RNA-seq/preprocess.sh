#!/bin/bash

## Download total_RNA-Seq metadata files from ENCODE
sh bin/download.metadata.sh "https://www.encodeproject.org/metadata/?status=released&internal_tags=ENTEx&assay_title=total+RNA-seq&type=Experiment" metadata.tsv

## Retrieve list of transcript quantification (v29)
grep tsv /gpfs/gibbs/pi/gerstein/yj329/epi/total_RNA-seq/preprocessing/metadata.tsv | grep -v hg19 | awk -F '\t' '($8=="total_RNA-seq" && $5=="transcript_quantifications")' | grep -i v29 | cut -f1,7-8,11,46 > /gpfs/gibbs/pi/gerstein/yj329/epi/total_RNA-seq/preprocessing/transcript_quantification.files.tsv

## Download transcript-level files (rsem and kallisto applied for transcript-level quantification)
cut -f1 /gpfs/gibbs/pi/gerstein/yj329/epi/total_RNA-seq/preprocessing/transcript_quantification.files.tsv | while read file; do echo -e "wget -P transcript_quantification https://www.encodeproject.org/files/"$file"/@@download/"$file".tsv; md5sum /gpfs/gibbs/pi/gerstein/yj329/epi/total_RNA-seq/preprocessing/transcript_quantification/"$file".tsv >> /gpfs/gibbs/pi/gerstein/yj329/epi/total_RNA-seq/preprocessing/transcript_quantification.md5sum"; done > transcript_quantification.download.job.list.txt

module load dSQ
dsq --job-file transcript_quantification.download.job.list.txt --mem-per-cpu 4g -t 1:00:00 --mail-type ALL

## Check md5sum
bin/join.py -a <(sed 's/\/.*\/transcript_quantification\///g;s/.tsv//g' transcript_quantification.md5sum | awk 'BEGIN{OFS="\t"}{print $2, $1}') -b transcript_quantification.files.tsv| awk '$(NF-1)!=$NF'

## Rename
Rscript bin/reformat.R

## Merge expression by tissue
python /gpfs/gibbs/pi/gerstein/yj329/epi/total_RNA-seq/merge_transcript_quantification.py
#!/bin/bash

## Download bigbed and bidwig from ENCODE (histone marks + DNAse + ATAC)
sh bin/download.metadata.sh "https://www.encodeproject.org/metadata/?status=released&internal_tags=ENTEx&assay_title=Histone+ChIP-seq&assay_title=DNase-seq&assay_title=ATAC-seq&type=Experiment&files.analyses.status=released&files.preferred_default=true" metadata.tsv
## list of bed_narrowPeak files
grep bed_narrowPeak metadata.tsv | grep -v hg19 | awk '($8=="DNase-seq" && $5=="peaks") || ($5=="pseudoreplicated_peaks" || $5=="replicated_peaks")' | grep -i encode4 | cut -f1,7-8,11,23,35,46 | awk 'BEGIN{FS=OFS="\t"}{if ($(NF-1)==1){$(NF-1)="b_1"} else if ($(NF-1)==2){$(NF-1)="c_2"} else {$(NF-1)="a_1,2"}; print $0}' | sort -k2,2 -k6,6 | sort -u -k2,2 | awk 'BEGIN{FS=OFS="\t"}{if ($3=="ATAC-seq" || $3=="DNase-seq"){$5="openChr"}; print $1, $2, $3, $4, $5, $7}' > bed.files.tsv
## download bed_narrowPeak files
cut -f1 bed.files.tsv | while read file; do echo -e "wget -P bed.files https://www.encodeproject.org/files/"$file"/@@download/"$file".bed.gz; md5sum bed.files/"$file".bed.gz >> bed.gz.md5sum"; done > bed.download.job.list.txt

module load dSQ
dsq --job-file bed.download.job.list.txt --mem-per-cpu 8g -t 24:00:00 --mail-type ALL --partition pi_gerstein,day,scavenge
sbatch dsq-bed.download.job.list-2024-05-07.sh
## check md5sum
bin/join.py -a <(sed 's/bed.files\///g;s/.bigWig//g' bed.gz.md5sum | awk 'BEGIN{OFS="\t"}{print $2, $1}') -b bed.files.tsv | awk '$(NF-1)!=$NF'

## prepare bigbed exp list
module load R
Rscript bin/bed.exp.R

## list of bigwig files
grep bigWig metadata.tsv | grep -v hg19 | awk 'BEGIN{FS=OFS="\t"}$5=="fold_change_over_control" || $5=="read-depth_normalized_signal"' |grep -i encode4 | cut -f1,7-8,11,23,35,46 | awk 'BEGIN{FS=OFS="\t"}{if ($(NF-1)==1){$(NF-1)="b_1"} else if ($(NF-1)==2){$(NF-1)="c_2"} else {$(NF-1)="a_1,2"}; print $0}' | sort -k2,2 -k6,6 | sort -u -k2,2 | awk 'BEGIN{FS=OFS="\t"}{if ($3=="ATAC-seq" || $3=="DNase-seq"){$5="openChr"}; print $1, $2, $3, $4, $5, $7}' > bigwig.files.tsv
# download bigwig files
cut -f1 bigwig.files.tsv | while read file; do echo -e "wget -P bigwig.files https://www.encodeproject.org/files/"$file"/@@download/"$file".bigWig; md5sum bigwig.files/"$file".bigWig >> bigwig.md5sum"; done > bigwig.download.job.list.txt

module load dSQ
dsq --job-file bigwig.download.job.list.txt --mem-per-cpu 10g -t 24:00:00 --mail-type ALL
sbatch dsq-bigwig.download.job.list-2024-02-29.sh

bin/join.py -a <(sed 's/bigwig.files\///g;s/.bigWig//g' bigwig.md5sum | awk 'BEGIN{OFS="\t"}{print $2, $1}') -b bigwig.files.tsv | awk '$(NF-1)!=$NF'


#!/bin/bash

## Download metadata files from ENCODE
sh bin/download.metadata.sh "https://www.encodeproject.org/metadata/?status=released&internal_tags=ENTEx&biosample_ontology.classification=tissue&assay_title=in+situ+Hi-C&type=Experiment&files.analyses.status=released&files.preferred_default=true" metadata.tsv

## list of hic, TAD, and loop files
grep hic metadata.tsv | grep -v hg19 |cut -f1,7-8,11,35,46 | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$6}' > hic.files.tsv
grep contact_domains metadata.tsv | grep -v hg19 |cut -f1,7-8,11,35,46 | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$6}' > TAD.files.tsv
grep loops metadata.tsv | grep -v hg19 |cut -f1,7-8,11,35,46 | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$6}' > loop.files.tsv
grep bigWig metadata.tsv | grep -v hg19 |cut -f1,7-8,11,35,46 | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$6}' > bigwig.files.tsv
## download bed_bed3+ files
cut -f1 hic.files.tsv | while read file; do echo -e "wget -P hic.files https://www.encodeproject.org/files/"$file"/@@download/"$file".hic; md5sum hic.files/"$file".hic >> hic.md5sum"; done > hic.download.job.list.txt
cut -f1 TAD.files.tsv | while read file; do echo -e "wget -P TAD.files https://www.encodeproject.org/files/"$file"/@@download/"$file".bedpe.gz; md5sum TAD.files/"$file".bedpe.gz >> TAD.gz.md5sum"; done > TAD.download.job.list.txt
cut -f1 loop.files.tsv | while read file; do echo -e "wget -P loop.files https://www.encodeproject.org/files/"$file"/@@download/"$file".bedpe.gz; md5sum loop.files/"$file".bedpe.gz >> loop.gz.md5sum"; done > loop.download.job.list.txt

cut -f1 bigwig.files.tsv | while read file; do echo -e "wget -P bigwig.files https://www.encodeproject.org/files/"$file"/@@download/"$file".bigWig; md5sum bigwig.files/"$file".bigWig >> bigwig.md5sum"; done > bigwig.download.job.list.txt



module load dSQ
dsq --job-file hic.download.job.list.txt --mem-per-cpu 8g -t 3:00:00 --mail-type ALL --partition pi_gerstein,day
dsq --job-file TAD.download.job.list.txt --mem-per-cpu 8g -t 3:00:00 --mail-type ALL --partition pi_gerstein,day
dsq --job-file loop.download.job.list.txt --mem-per-cpu 8g -t 3:00:00 --mail-type ALL --partition pi_gerstein,day
dsq --job-file bigwig.download.job.list.txt --mem-per-cpu 8g -t 3:00:00 --mail-type ALL --partition pi_gerstein,day

sbatch dsq-hic.download.job.list-2024-03-15.sh
sbatch dsq-TAD.download.job.list-2024-03-15.sh
sbatch dsq-loop.download.job.list-2024-03-15.sh
sbatch dsq-bigwig.download.job.list-2024-03-17.sh

## check md5sum
bin/join.py -a <(sed 's/bed.files\///g;s/.bigWig//g' bed.gz.md5sum | awk 'BEGIN{OFS="\t"}{print $2, $1}') -b bed.files.tsv | awk '$(NF-1)!=$NF'

## prepare big exp list
module load R
Rscript bin/bigwig.exp.R

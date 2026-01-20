#!/bin/bash
# Download files from ENCODE (HiC)
sh bin/download.metadata.sh "https://www.encodeproject.org/metadata/?control_type%21=%2A&status=released&perturbed=false&assay_title=intact+Hi-C&assay_title=in+situ+Hi-C&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&type=Experiment&files.analyses.status=released&files.preferred_default=true" metadata.tsv
# list of bedpe and HiC files (loops, not contact domains)
grep bedpe metadata.tsv | grep GRCh38 | awk -F '\t' '($5=="loops" && $12=="tissue")' | cut -f1,7-8,11,46 > bedpe.files.tsv
cut -f1 bedpe.files.tsv | while read file; do echo "wget -P loop.files https://www.encodeproject.org/files/"$file"/@@download/"$file".bedpe.gz; md5sum loop.files/"$file".bedpe.gz >> bedpe.gz.md5sum"; done > bedpe.download.job.list.txt

module load dSQ
dsq --job-file bedpe.download.job.list.txt --mem-per-cpu 8g -t 24:00:00 --mail-type ALL --partition pi_gerstein,day,scavenge
sbatch dsq-bedpe.download.job.list-2024-06-10.sh

## check md5sum
bin/join.py -a <(sed 's/bedpe.files\///g;s/.bedpe.gz//g' bedpe.gz.md5sum | awk 'BEGIN{OFS="\t"}{print $2, $1}') -b bedpe.files.tsv | awk '$(NF-1)!=$NF'
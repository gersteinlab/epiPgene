# generate promoters (1000 bp around TSS) of protein coding genes and lncRNAs from the GTF file
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' gencode.v29.annotation.gtf | gtf2bed - > gencode.v29.annotation.bed

cat gencode.v29.annotation.bed | awk 'BEGIN{OFS=FS="\t"}{split($NF, a, ";"); if($8=="transcript") {if($6=="+") {start=$2-1000; end=$2+1000;} else {if($6=="-") start=$3-1000; end=$3+1000; } if(start<0) start=0; print $1,start,end,$4,$5,$6, a[2], a[3]}}'| sed -E 's/transcript_id|gene_type//g' | sed 's/[" ]//g' > GRCh38_tss_1kb.v29.bed 
cat gencode.v29.annotation.gtf | grep -v '#' | sed 's/"/\t/g' | awk -F '\t' '{if($3=="transcript") print $10"\t"$12"\t"$14}' > geneType.v29.txt
cat gencode.v29.annotation.bed | awk 'BEGIN{OFS=FS="\t"}{split($NF, a, ";"); if($8=="transcript") print $1,$2,$3,$4,$5,$6,a[2],a[3]}' | sed -E 's/transcript_id|gene_type//g' | sed 's/[" ]//g' > GRCh38_transcripts.v29.bed

# finding promoters that overlap with others
bedtools intersect -a GRCh38_tss_1kb.v29.bed -b GRCh38_tss_1kb.v29.bed -wao | awk -F '\t' '{if($4!=$12) print $0}' | python findOverlappingPromoters.py
# transcripts whose promoters overlap with other promoters
awk '{print $1; print $2}' overlapping_promoters.txt | sort | uniq > blacklist.txt
# upstream TSS only
cat gencode.v29.annotation.bed | awk 'BEGIN{OFS=FS="\t"}{split($NF, a, ";"); if($8=="transcript") {if($6=="+") {start=$2-1000; end=$2;} else {if($6=="-") start=$3; end=$3+1000; } if(start<0) start=0;print $1,start,end,$4,$5,$6, a[2], a[3]}}'| sed -E 's/transcript_id|gene_type//g' | sed 's/[" ]//g' > GRCh38_tss_up1kb_only.v29.bed
#!/bin/bash

bedtools nuc -fi hg38.fa -bed GRCh38_tss_upstream.v29.bed > GRCh38_tss_upstream.v29.GC.txt
bedtools nuc -fi hg38.fa -bed GRCh38_tss_downstream.v29.bed > GRCh38_tss_downstream.v29.GC.txt
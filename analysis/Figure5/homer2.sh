#!/bin/bash
#SBATCH --job-name=motif
#SBATCH --partition=bigmem,day,scavenge
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --requeue
#SBATCH --mem=128g
#SBATCH --mail-type=ALL

module load miniconda
conda activate py3

# processed
findMotifsGenome.pl processed_pgene_promoter_regions.bed hg38 processed_random -p 16 -size 200 -S 10
findMotifsGenome.pl processed_pgene_promoter_regions.bed hg38 processed_pcg -p 16 -size 200 -bg protein_coding_promoter_regions.bed -S 10
# unprocessed
findMotifsGenome.pl unprocessed_pgene_promoter_regions.bed hg38 unprocessed_random -p 16 -size 200 -S 10
findMotifsGenome.pl unprocessed_pgene_promoter_regions.bed hg38 unprocessed_pcg -p 16 -size 200 -bg protein_coding_promoter_regions.bed -S 10
# PCG
findMotifsGenome.pl protein_coding_promoter_regions.bed hg38 pcg_random -p 16 -size 200 -S 10
# processed vs unprocessed
findMotifsGenome.pl processed_pgene_promoter_regions.bed hg38 processed_unprocessed -p 16 -size 200 -bg unprocessed_pgene_promoter_regions.bed -S 10 -h

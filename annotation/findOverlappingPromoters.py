import pandas as pd
import sys

def parse_bed_line(line):
    fields = line.strip().split()
    chrom, start, end, gene_id, score, strand, transcript_id = fields[:7]
    return chrom, int(start), int(end), gene_id, score, strand, transcript_id

def get_5prime(promoter):
    _, start, end, _, _, strand, _ = promoter
    return start if strand == '+' else end

def get_3prime(promoter):
    _, start, end, _, _, strand, _ = promoter
    return end if strand == '+' else start

def distance(coord1, coord2):
    return abs(coord1 - coord2)

def is_unidirectional(promoter1, promoter2):
    return promoter1[5] == promoter2[5] 

def is_head_to_head(promoter1, promoter2):
    return promoter1[5] != promoter2[5] and distance(get_5prime(promoter1), get_5prime(promoter2)) < distance(get_3prime(promoter1), get_3prime(promoter2))

def is_tail_to_tail(promoter1, promoter2):
    return promoter1[5] != promoter2[5] and distance(get_5prime(promoter1), get_5prime(promoter2)) > distance(get_3prime(promoter1), get_3prime(promoter2))

res = []
flag = set()

## BEDtools result from the stdin
for line in sys.stdin:
    parts = line.strip().split("\t")
    promoter1 = parse_bed_line("\t".join(parts[0:7]))
    promoter2 = parse_bed_line("\t".join(parts[8:15]))
    
    if is_unidirectional(promoter1, promoter2):
        overlap = 'unidirectional'
    elif is_head_to_head(promoter1, promoter2):
        overlap = 'head-to-head'
    elif is_tail_to_tail(promoter1, promoter2):
        overlap = 'tail-to-tail'

    item = tuple(sorted([promoter1[-1], promoter2[-1]]))
    if item not in flag:
        flag.add(item)
        res.append([promoter1[-1], promoter2[-1], overlap])

data = pd.DataFrame(res)
data.to_csv('/gpfs/gibbs/pi/gerstein/yj329/epi/gencode_v29/overlapping_promoters.txt', index = False, header = False, sep = '\t')
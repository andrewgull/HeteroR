#!/usr/bin/env python3
from Bio import SeqIO


with open("/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_perfect_repeats_GRF_test/perfect.spacer.id") as f:
    lines = [line.rstrip() for line in f.readlines()]

# a line in lines looks like this: >1:0-190543:951:190482:14m'
# >1:0-190543 - fasta id
# 951:190482:14m - start, end and length of a repeat

# read GBK to look inside
gbk = [rec for rec in SeqIO.parse("/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_genomic.gbk", "genbank")]

gbk[0].features[0]

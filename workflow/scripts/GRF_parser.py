#!/usr/bin/env python3
from Bio import SeqIO
import pandas as pd

with open("/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_perfect_repeats_GRF_test/perfect.spacer.id") as f:
    lines = [line.rstrip() for line in f.readlines()]

# a line in lines looks like this: >1:0-190543:951:190482:14m'
# >1:0-190543 - fasta id
# 951:190482:14m - start of the 1st match, end of the 2nd match and length of the match

# read GBK to look inside
gbk = [rec for rec in SeqIO.parse("/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_genomic.gbk", "genbank")]

gbk[0].features[0]

in_rgi = "/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_rgi_table.txt"
rgi = pd.read_csv(in_rgi, sep="\t")

# parse gff or genebank

# you have info on genes
# see the copybook

#!/usr/bin/env python3

"""
renames IDs in prokka-made GBK annotation files
so that these GBK will be compatible with assembly's IDs and genomic browsers
"""

from Bio import SeqIO
import os


# in_gbk = "/home/andrei/Data/HeteroR/results/annotations/DA62886/prokka/DA62886_genomic.gbk"
# out_gbk = "/home/andrei/Data/HeteroR/results/annotations/DA62886/prokka/DA62886_genomic_renamed.gbk"
# input is a directory
in_gbk_dir = snakemake.input[0]
out_gbk = snakemake.output[0]
filename = snakemake.params[0]

# make a filename
in_gbk = os.path.join(in_gbk_dir, filename)

gbk = [rec for rec in SeqIO.parse(in_gbk, "genbank")]
for rec in gbk:
    # ID name pattern is like this: LOOOHHPE_1
    # I need number after the underscore - this is chromosome's name in assembly file
    new_id = rec.id.split("_")[-1]
    rec.id = new_id

SeqIO.write(gbk, out_gbk, "genbank")

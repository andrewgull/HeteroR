#!/usr/bin/env python3

"""
to covert rgi output table to bed/gff file
"""

import pandas as pd
from BCBio import GFF


in_rgi = "/home/andrei/Data/HeteroR/results/resistance_genes/DA62886/rgi_table.txt"
# no coordinates of genes in this table
# take them from GFF?
in_gff = "/home/andrei/Data/HeteroR/results/annotations/DA62886/prokka/DA62886_genomic.gff"
in_gbk = "/home/andrei/Data/HeteroR/results/annotations/DA62886/prokka/DA62886_genomic.gbk"

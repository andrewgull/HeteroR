"""
sequence of commands/a script prototype
to join RGI and RGF results
requires: bcbio-gff and gffutils
"""


from Bio import SeqIO
import pandas as pd
from BCBio import GFF

# cd /home/andrei/Data/HeteroR/test_dir/GRF

# rgi results
rgi = pd.read_csv("DA62886_rgi_table.tsv", sep="\t")

# parse gbk file
gbk = [rec for rec in SeqIO.parse("DA62886_genomic.gbk", "genbank")]
gbk_dict = SeqIO.to_dict(SeqIO.parse("DA62886_genomic.gbk", "gb"))

# it contains 3 records: a chromosome and two plasmids, they are accessible through .features
# get chromosomal genes

# chrom_genes = [item for item in gbk[0].features if item.type == 'gene']  # 4947

# try GFF file instead
in_file = "DA62886_genomic.gff"

with open(in_file) as f:
    gff = [rec for rec in GFF.parse(f)]

chromosome_genes = [feature for feature in gff[0].features if feature.type == "gene"]  # here we have IDs and positions

# filter the genes with resistance - both in RGI and GBK
rgi["ORF_ID"]

# get their coordinates - from GBK
# get ± 100 kb region for each gene


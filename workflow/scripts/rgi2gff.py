#!/usr/bin/env python3

"""
to covert rgi output table to bed/gff file
"""

import pandas as pd
from Bio import SeqIO


# in_rgi = "/home/andrei/Data/HeteroR/results/resistance_genes/DA62886/rgi_table.txt"
# in_gbk = "/home/andrei/Data/HeteroR/results/annotations/DA62886/prokka/DA62886_genomic.gbk"
# out_gbk = "/home/andrei/Data/HeteroR/results/annotations/DA62886/prokka/DA62886_resistance_genes.gbk"
# filter_condition = "Loose"

in_rgi = snakemake.input[0]
in_gbk = snakemke.input[1]
out_gbk = snakemake.output[0]
filter_condition = snakemake.params[0]

# no coordinates of genes in this table
rgi = pd.read_csv(in_rgi, delimiter="\t")

# take them from GFF or GBK?
gbk = [rec for rec in SeqIO.parse(in_gbk, "genbank")]
# ORF_IDs from rgi are present in GBK file, find them in GBK object

# iterate through rows and collect genes with features
rgi_not_loose = rgi[rgi.Cut_Off != filter_condition]
rgi_not_loose_ids = [row[1][0].split(" ")[0] for row in rgi_not_loose.iterrows()]
# keep more information to add to qualifiers
# ORF_ID(0), SNP(12), Drug Class(14), AMR Gene Family(16)
rgi_not_loose_dict = dict()
for row in rgi_not_loose.iterrows():
    id = row[1][0].split(" ")[0]
    if id in rgi_not_loose_ids:
        rgi_not_loose_dict[id] = [row[1][12], row[1][14], row[1][16]]

# rgi_not_loose_tuples = [(row[1][0].split(" ")[0], row[1][12], row[1][14], row[1][16]) for row in rgi_not_loose.iterrows()]

for record in gbk:
    # change record.id to compatible with assembly's records names
    record.id = record.id.split("_")[-1]
    rgi_from_gbk = [f for f in record.features if f.type == 'gene' and f.qualifiers['locus_tag'][0] in rgi_not_loose_ids]
    # rewrite record's features
    record.features = rgi_from_gbk

# remove records with no features (i.e. no resistance genes)
rgi_gbk = [record for record in gbk if len(record.features) > 0]

# add more info to qualifiers
for record in rgi_gbk:
    for feature in record.features:
        qualifiers_to_add = rgi_not_loose_dict[feature.qualifiers['locus_tag'][0]]
        feature.qualifiers['SNP'] = qualifiers_to_add[0]
        feature.qualifiers['drug_class'] = qualifiers_to_add[1]
        feature.qualifiers['AMR_gene_family'] = qualifiers_to_add[2]

# write to a file
SeqIO.write(rgi_gbk, out_gbk, "genbank")

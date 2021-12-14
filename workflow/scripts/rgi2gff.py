#!/usr/bin/env python3

"""
to covert rgi output table to bed/gff file
"""

import pandas as pd
from BCBio import GFF
from Bio import SeqIO


in_rgi = "/home/andrei/Data/HeteroR/results/resistance_genes/DA62886/rgi_table.txt"

# take them from GFF?
in_gff = "/home/andrei/Data/HeteroR/results/annotations/DA62886/prokka/DA62886_genomic.gff"
in_gbk = "/home/andrei/Data/HeteroR/results/annotations/DA62886/prokka/DA62886_genomic.gbk"

# no coordinates of genes in this table
rgi = pd.read_csv(in_rgi, delimiter="\t")

# take them from GFF or GBK?
gbk = [rec for rec in SeqIO.parse(in_gbk, "genbank")]
# ORF_IDs from rgi are present in GBK file, find them in GBK object
for rec in gbk:
    print(rec.id)
    for f in rec.features:
        print(f)

# type: mRNA
# location: [2092185:2093025](+)
# qualifiers:
#     Key: gene, Value: ['hdfR_2']
#     Key: locus_tag, Value: ['LOOOHHPE_01937']
#     Key: product, Value: ['hypothetical protein']
# type: CDS
# location: [2092185:2093025](+)
# qualifiers:
#     Key: codon_start, Value: ['1']
#     Key: gene, Value: ['hdfR_2']
#     Key: inference, Value: ['ab initio prediction:Prodigal:002006', 'protein motif:HAMAP:MF_01233']
#     Key: locus_tag, Value: ['LOOOHHPE_01937']
#     Key: product, Value: ['HTH-type transcriptional regulator HdfR']
#     Key: protein_id, Value: ['UU:LOOOHHPE_01937']
#     Key: transl_table, Value: ['11']
#     Key: translation, Value: ['MDTELLKTFLEVSRTRHFGRAA...

qual = [f.qualifiers for f in gbk[0].features]
# qual[0]
# first qualifier is entire chromosome\plasmid qualifier (type='source'), has no locus_tags
# OrderedDict([('organism', ['Escherichia coli']),
#              ('mol_type', ['genomic DNA']),
#              ('strain', ['DA62886']),
#              ('db_xref', ['taxon:562'])])

feat = [f for f in gbk[0].features]
# a single gene has 3 features: 1 type='gene', 1 type'mRNA', 1 type='CDS'
genes = [f for f in feat if f.type == 'gene']
tags = [gene.qualifiers['locus_tag'][0] for gene in genes]

# iterate through rows and collect genes with features
rgi_not_loose = rgi[rgi.Cut_Off != "Loose"]
rgi_from_gbk = list()
for row in rgi_not_loose.iterrows():
    orf_id = row[1][0].split(" ")[0]
    for gene in genes:
        if orf_id == gene.qualifiers['locus_tag'][0]:
            rgi_from_gbk.append(gene)

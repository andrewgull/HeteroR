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
rgi = pd.read_csv("DA62886_rgi_table.tsv", sep="\t")  # all genes are in this list???

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
# naive looping approach
resistance_genes_coords = list()
for orf in rgi["ORF_ID"]:
    for gene in chromosome_genes:
        if orf.split(" ")[0] in gene.id:
            resistance_genes_coords.append(gene)

# get their coordinates - from GBK
# get ± 100 kb region for each gene

genome_end = 5131220
span = 100000
for gene in resistance_genes_coords:
    start = int(gene.location.start)
    end = int(gene.location.end)
    span_start = start - span
    # check if it's negative
    #if span_start < 0:
    #    span_start = genome_end + span_start
    span_end = end + span
    row = [gene.id, start, end, span_start, span_end, int(gene.location.strand)]
    print(row)
genes_and_spans = pd.DataFrame()
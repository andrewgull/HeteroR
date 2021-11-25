"""
sequence of commands/a script prototype
to join RGI and RGF results
requires: bcbio-gff and gffutils
"""


from Bio import SeqIO
import pandas as pd
from BCBio import GFF

# cd /home/andrei/Data/HeteroR/test_dir/GRF
# VARIABLES
in_rgi = "DA62886_rgi_table.tsv"
in_gff = "DA62886_genomic.gff"
genome_len = 5131220  # must be determined for each particular genome
span_len = 100000

# rgi results - resistance information
rgi = pd.read_csv(in_rgi, sep="\t")
# rgi["Cut_Off"].unique()
# some genes are 'Loose', leave 'Strict' and 'Perfect' only
rgi_notLoose = rgi[rgi["Cut_Off"] != "Loose"]

# parse GFF annotation file
# to retrieve chromosomal genes' coordinates and strand
with open(in_gff) as f:
    gff = [rec for rec in GFF.parse(f)]

chromosome_genes = [feature for feature in gff[0].features if feature.type == "gene"]  # here we have IDs and positions
chromosome_id = [gff[0].id]

# and plasmid genes' coordinates and strand
plasmid_genes = list()
plasmid_id = list()
for i in range(1, len(gff)):
    plasmid_genes.append([feature for feature in gff[i].features if feature.type == "gene"])
    plasmid_id.append(gff[i].id)

# find resistance genes' coords: RGI - resistance, GBK - coords
# naive looping approach: 61*4947 comparisons
resistance_genes_coords = list()
for orf in rgi_notLoose["ORF_ID"]:
    for gene in chromosome_genes:
        if orf.split(" ")[0] in gene.id:
            resistance_genes_coords.append(gene)

print("Coordinates of %i resistance genes not found" % (len(rgi_notLoose) - len(resistance_genes_coords)))

# get their coordinates - from GBK
# get ± 100 kb region for each gene

genes_lol = list()

for gene in resistance_genes_coords:
    start = int(gene.location.start)
    end = int(gene.location.end)
    span_start = start - span_len
    # check if it's negative
    #if span_start < 0:
    #    span_start = genome_len + span_start
    span_end = end + span_len
    row = [gene.id, start, end, span_start, span_end, int(gene.location.strand)]
    genes_lol.append(row)

# a table with span and resistance gene coordinates
rg_spans_and_coords = pd.DataFrame(columns=["gene_id", "gene_start", "gene_end", "span_start", "span_end", "strand"],
                           data=genes_lol)

# turn this data frame to a bed file - but how to deal with negative coordinates?
# https://bacteria.ensembl.org/info/website/upload/bed.html
# Use bed tools - create a bed file
# something should be done with the negative coordinates

rg_spans_and_coords_positive = rg_spans_and_coords[rg_spans_and_coords["span_start"] >= 0]
rg_spans_and_coords_negative = rg_spans_and_coords[rg_spans_and_coords["span_start"] < 0]

# making a bed file for ranges not crossing oriC
bed_file = pd.DataFrame()

bed_file["range_start"] = rg_spans_and_coords_positive["span_start"]
bed_file["range_end"] = rg_spans_and_coords_positive["span_end"]
bed_file["name"] = rg_spans_and_coords_positive["gene_id"]
bed_file["score"] = "0"
bed_file["strand"] = rg_spans_and_coords_positive["strand"]



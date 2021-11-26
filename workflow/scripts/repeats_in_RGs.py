"""
sequence of commands/a script prototype
to join RGI and RGF results
requires: bcbio-gff and gffutils
"""


from Bio import SeqIO
import pandas as pd
from BCBio import GFF


def make_bed_file(gff_record, rgi_dataframe):
    """
    makes bed file for genomic ranges with resistance genes
    """
    genes = [feature for feature in gff_record.features if feature.type == "gene"]  # here we have IDs and positions
    item_id = gff_record.id

    # find resistance genes' coords: RGI - resistance, GBK - coords
    # naive looping approach: 61*4947 comparisons
    resistance_genes_coords = list()
    for orf in rgi_dataframe["ORF_ID"]:
        for gene in genes:
            if orf.split(" ")[0] in gene.id:
                resistance_genes_coords.append(gene)

    message = "In record %s %i resistance genes not found" % (item_id, len(rgi_dataframe) - len(resistance_genes_coords))

    # get their coordinates - from GBK
    # get ± 100 kb region for each gene

    rg_collection = list()

    for gene in resistance_genes_coords:
        start = int(gene.location.start)
        end = int(gene.location.end)
        span_start = start - span_len
        # check if it's negative
        #if span_start < 0:
        #    span_start = genome_len + span_start
        span_end = end + span_len
        row = [gene.id, start, end, span_start, span_end, int(gene.location.strand)]
        rg_collection.append(row)

    # a table with span and resistance gene coordinates
    rg_ranges = pd.DataFrame(columns=["gene_id", "gene_start", "gene_end", "span_start", "span_end", "strand"],
                             data=rg_collection)

    # turn this data frame to a bed file - but how to deal with negative coordinates?
    # https://bacteria.ensembl.org/info/website/upload/bed.html
    # Use bed tools - create a bed file
    # something should be done with the negative coordinates
    # split ranges not crossing oriC and not crossing it
    rg_ranges_pos = rg_ranges[rg_ranges["span_start"] >= 0]
    rg_ranges_neg = rg_ranges[rg_ranges["span_start"] < 0]

    # making a bed file for ranges not crossing oriC
    bed_file = pd.DataFrame()
    bed_file["range_start"] = rg_ranges_pos["span_start"]
    bed_file["range_end"] = rg_ranges_pos["span_end"]
    bed_file["name"] = rg_ranges_pos["gene_id"]
    bed_file["score"] = "0"
    bed_file["strand"] = rg_ranges_pos["strand"]
    return bed_file, message


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

# iterate through chromosome and plasmids
for item in gff:
    ranges_bed, bed_message = make_bed_file(item, rgi_notLoose)

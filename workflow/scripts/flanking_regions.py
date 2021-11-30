"""
sequence of commands/a script prototype
to join RGI and RGF results
requires: bcbio-gff and gffutils
"""
import pybedtools
from Bio import SeqIO
import pandas as pd
from BCBio import GFF
from pybedtools import BedTool


def make_bed_file(gff_record, rgi_dataframe, dna_len, span_len, circular):
    """
    makes bed file for genomic ranges with resistance genes
    it takes into account the length of dna record and range
    :param gff_record: GFF-object, list of gff records of a particular assembly
    :param rgi_dataframe: pandas DataFrame, pd DataFrame of RGI results
    :param dna_len: integer, length of a current record (chromosome or plasmid)
    :param span_len: integer, length of up- and downstream region flanking a resistance gene
    :param circular: boolean, is current dna record circular or not?
    :return: bed formatted pandas DataFrame, negative span pandas DataFrame, message with number of found genes
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

    message = "In record %s %i of %i resistance genes found" % (item_id, len(resistance_genes_coords), len(rgi_dataframe))

    # get their coordinates - from GBK
    # get ± 100 kb region for each gene

    rg_collection = list()
    if not circular:
        # then don't care about ranges crossing oriC
        # range coordinate will always be less than span length
        for gene in resistance_genes_coords:
            chrom_name = item_id.split("_")[-1]  # goes to the first column of abed file
            start = int(gene.location.start)
            end = int(gene.location.end)
            # check left end
            if start < span_len:
                # to include leftmost letter subtract 1
                span_start = 0
            else:
                span_start = start - span_len - 1
            # check right end: bedtools includes it
            if end + span_len > dna_len:
                span_end = dna_len
            else:
                span_end = end + span_len
            # make a future table row
            row = [chrom_name, gene.id, start, end, span_start, span_end, int(gene.location.strand)]
            rg_collection.append(row)

        # a table with span and resistance gene coordinates
        rg_ranges = pd.DataFrame(columns=["chrom", "gene_id", "gene_start", "gene_end", "span_start", "span_end", "strand"],
                                 data=rg_collection)

        # split ranges not crossing oriC and not crossing it
        rg_ranges_pos = rg_ranges[rg_ranges["span_start"] >= 0]
        rg_ranges_neg = rg_ranges[rg_ranges["span_start"] < 0]

        # making a bed file for ranges not crossing oriC
        bed_file = pd.DataFrame()
        bed_file["chrom"] = rg_ranges_pos["chrom"]
        bed_file["range_start"] = rg_ranges_pos["span_start"]
        bed_file["range_end"] = rg_ranges_pos["span_end"]
        bed_file["name"] = rg_ranges_pos["gene_id"]
        bed_file["score"] = "0"
        bed_file["strand"] = rg_ranges_pos["strand"]
    else:
        # this line doesn't do anything
        bed_file, rg_ranges_neg, message = handle_negative_coords()

    return bed_file, rg_ranges_neg, message


def handle_negative_coords():
    """
    do something with negative coordinates
    """
    return 1, 1, 1


# cd /home/andrei/Data/HeteroR/test_dir/GRF
# VARIABLES NON CIRCULAR CHROMOSOME
in_rgi = "DA62886_rgi_table.txt"
in_gff = "DA62886_genomic.gff"
in_assembly = "DA62886_assembly.fasta"
regions_bed_output = "regions_output.bed"
regions_fasta_output = "regions_output.fasta"

# VARIABLES CIRCULAR EVERYTHING
in_rgi_circ = "DA63004_rgi_table.txt"
in_gff_circ = "DA63004_genomic.gff"
in_assembly_circ = "DA63004_assembly.fasta"

# genome_len = 5131220  # must be determined for each particular genome
range_len = 100000

# rgi results - resistance information
rgi = pd.read_csv(in_rgi, sep="\t")
# rgi["Cut_Off"].unique()
# some genes are 'Loose', leave 'Strict' and 'Perfect' only
rgi_notLoose = rgi[rgi["Cut_Off"] != "Loose"]

# parse GFF annotation file
# to retrieve chromosomal genes' coordinates and strand
with open(in_gff) as f:
    gff = [rec for rec in GFF.parse(f)]

# read joined assembly file
assembly = [rec for rec in SeqIO.parse(in_assembly, "fasta")]

# iterate through chromosome and plasmids
for i in range(len(gff)):
    # TODO: negative coordinates?
    record_len = len(assembly[i].seq)
    # find is it circular
    if "circular=true" in assembly[i].description:
        circ = True
    else:
        circ = False

    ranges_bed, negative_coords, bed_message = make_bed_file(gff_record=gff[i], rgi_dataframe=rgi_notLoose,
                                                             dna_len=record_len, span_len=range_len, circular=circ)
    # no merge needed here
    ranges_bed.to_csv(regions_bed_output, sep="\t", index=False, header=False)
    # cut regions using bedtools
    bed_file = BedTool(regions_bed_output)
    # write fasta regions to a file
    bedtool_write = pybedtools.bedtool.BedTool.sequence(bed_file, fi=in_assembly, fo=regions_fasta_output)


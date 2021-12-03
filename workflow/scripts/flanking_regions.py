#!/usr/bin/env python3
"""
script to create a bed files for a given assembly and get regions from bed in fasta format
requires: bcbio-gff and gffutils
"""

from Bio import SeqIO
import pandas as pd
from BCBio import GFF


def make_bed(collection, score):
    """
    :param collection: pd DataFrame of ranges for a bed file
    :param score: score column value for bed file
    :return: bed formatted pd DataFrame
    """
    bed = pd.DataFrame()
    bed["chrom"] = collection["chrom"]
    bed["range_start"] = collection["span_start"]
    bed["range_end"] = collection["span_end"]
    bed["name"] = collection["gene_id"]
    bed["score"] = str(score)
    bed["strand"] = collection["strand"]
    return bed


def make_bed_file_for_rg(gff_record, rgi_dataframe, dna_len, span_len, circular):
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

    msg = "In record %s %i of %i resistance genes found\n" % (item_id, len(resistance_genes_coords), len(rgi_dataframe))

    # get their coordinates - from GBK
    # get ± 100 kb region for each gene

    rg_collection = list()
    # if not circular:  # TODO: SOLVE THE CASE WHEN IT'S CIRCULAR!
    # then don't care about ranges crossing oriC
    # range coordinate will always be less than span length
    for gene in resistance_genes_coords:
        chrom_name = item_id.split("_")[-1]  # goes to the first column of a bed file
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
    bed_dataframe = make_bed(rg_ranges_pos, score=0)

    return bed_dataframe, rg_ranges_neg, msg


def handle_circular_records():
    """
    do something with negative coordinates
    """
    return 1, 1, 1


def join_bed_files(bed_files_list):
    """
    join bed files
    """
    return pd.concat(bed_files_list)

# cd /home/andrei/Data/HeteroR/test_dir/GRF
# VARIABLES TEST NON CIRCULAR CHROMOSOME
# in_rgi = "DA62886_rgi_table.txt"
# in_gff = "DA62886_genomic.gff"
# in_assembly = "DA62886_assembly.fasta"
# regions_bed_output = "regions_output.bed"
# regions_fasta_output = "regions_output.fasta"
#
# # VARIABLES TEST CIRCULAR EVERYTHING
# in_rgi_circ = "DA63004_rgi_table.txt"
# in_gff_circ = "DA63004_genomic.gff"
# in_assembly_circ = "DA63004_assembly.fasta"


# VARIABLES FOR SNAKEMAKE
in_assembly = snakemake.input[0]
in_gff = snakemake.input[1]
in_rgi = snakemake.input[2]
range_len = int(snakemake.params[0])
min_plasmid_size = int(snakemake.params[1])
regions_bed_output = snakemake.output[0]

# write things to log
with open(snakemake.log[0], "w") as log:
    # rgi results - resistance information
    rgi = pd.read_csv(in_rgi, sep="\t")
    # rgi["Cut_Off"].unique()
    # some genes are 'Loose', leave 'Strict' and 'Perfect' only
    rgi_notLoose = rgi[rgi["Cut_Off"] != "Loose"]

    # parse GFF annotation file
    # to retrieve chromosomal genes' coordinates and strand
    with open(in_gff) as f:
        gff = [rec for rec in GFF.parse(f) if len(rec.seq) >= min_plasmid_size]
    gff_ids = [rec.id.split("_")[-1] for rec in gff]
    # read joined assembly file as dict
    assembly = SeqIO.to_dict(SeqIO.parse(in_assembly, "fasta"))
    # filter it because not all of the records present in GFF
    assembly_filtered = [assembly[key] for key in gff_ids]

    # iterate through chromosome and plasmids
    bed_list = list()
    messages = list()
    for i in range(len(gff)):
        # TODO: negative coordinates?
        # i in assembly does not correspond to i in gff!
        record_len = len(assembly_filtered[i].seq)
        record_id = assembly_filtered[i].id
        # find is it circular
        if "circular=true" in assembly_filtered[i].description:
            circ = True
        else:
            circ = False

        ranges_bed, negative_coords, bed_message = make_bed_file_for_rg(gff_record=gff[i], rgi_dataframe=rgi_notLoose,
                                                                        dna_len=record_len, span_len=range_len,
                                                                        circular=circ)
        bed_list.append(ranges_bed)
        messages.append(bed_message)

    joined_bed_dataframe = join_bed_files(bed_list)
    # write dataframe to a BED file
    joined_bed_dataframe.to_csv(regions_bed_output, sep="\t", index=False, header=False)
    # write messages to log
    log.writelines(messages)

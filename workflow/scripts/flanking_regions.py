#!/usr/bin/env python3
"""
script to create a bed files for a given assembly and get regions from bed in fasta format
requires: bcbio-gff and gffutils
"""

from Bio import SeqIO
import pandas as pd
from BCBio import GFF
import os
import numpy as np


def make_bed(collection, score):
    """
    :param collection: pd DataFrame of ranges for a bed_normal file
    :param score: score column value for bed_normal file
    :return: bed_normal formatted pd DataFrame
    """
    # make bed_normal for normal coordinates
    bed_normal = pd.DataFrame()
    bed_normal["chrom"] = collection["chrom"]
    bed_normal["range_start"] = collection["span_start"]
    bed_normal["range_end"] = collection["span_end"]
    bed_normal["name"] = collection["gene_id"]
    bed_normal["score"] = str(score)
    bed_normal["strand"] = collection["strand"]
    # make bed_5 for spans crossing 5-end
    bed_5 = pd.DataFrame()
    bed_5["chrom"] = collection["chrom"]
    bed_5["range_start"] = collection["span_over_5_start"]
    bed_5["range_end"] = collection["span_over_5_end"]
    bed_5["name"] = collection["gene_id"]
    bed_5["score"] = str(score)
    bed_5["strand"] = collection["strand"]
    # drop NaNs in range_start (i.e. span_over_5_start)
    bed_5 = bed_5[bed_5.range_start.notnull()]
    # make bed_3 for spans crossing 3-end
    bed_3 = pd.DataFrame()
    bed_3["chrom"] = collection["chrom"]
    bed_3["range_start"] = collection["span_over_3_start"]
    bed_3["range_end"] = collection["span_over_3_end"]
    bed_3["name"] = collection["gene_id"]
    bed_3["score"] = str(score)
    bed_3["strand"] = collection["strand"]
    # drop NaNs
    bed_3 = bed_3[bed_3.range_end.notnull()]
    # join them
    # bed = pd.concat([bed_normal, bed_5, bed_3])
    return bed_normal, bed_5, bed_3


def make_bed_file_for_rg(gff_record, rgi_dataframe, dna_len, span_len):
    """
    makes bed files for genomic ranges with resistance genes
    it takes into account the length of dna record and range
    :param gff_record: GFF-object, list of gff records of a particular assembly
    :param rgi_dataframe: pandas DataFrame, pd DataFrame of RGI results
    :param dna_len: integer, length of a current record (chromosome or plasmid)
    :param span_len: integer, length of up- and downstream region flanking a resistance gene
    :return: a list of bed formatted pandas DataFrame, message with number of found genes
    """
    genes = [feature for feature in gff_record.features if feature.type == "gene"]  # here we have IDs and positions
    item_id = gff_record.id

    # find median gene length
    genes_lengths = [gene.location.end.position - gene.location.start.position for gene in genes]
    if len(genes) > 0:
        median_gene = pd.Series(genes_lengths).median()
    else:
        # it can be 0 genes in record
        median_gene = 0
    # adjust span_len according to the record length
    # span_len = (record_len - median_gene)/2, if 2*span_len > record_len
    record_length = len(gff_record.seq)
    if record_length < span_len * 2:
        span_len = round((record_length - median_gene) / 2)
        msg_len = "span length was adjusted to %i in record %s\n" % (span_len, item_id)
    else:
        msg_len = ""

    # find resistance genes' coords: RGI - resistance, GBK - coords
    # naive looping approach: 61*4947 comparisons
    resistance_genes_coords = list()
    for orf in rgi_dataframe["ORF_ID"]:
        for gene in genes:
            if orf.split(" ")[0] in gene.id:
                resistance_genes_coords.append(gene)

    msg_count = "In record %s %i of %i resistance genes found\n" % (item_id, len(resistance_genes_coords),
                                                                    len(rgi_dataframe))

    # get their coordinates - from GBK
    # get ± 100 kb region for each gene
    if len(resistance_genes_coords) > 0:
        rg_collection = list()
        # if not circular:
        # then don't care about ranges crossing oriC
        # range coordinate will always be less than span length
        for gene in resistance_genes_coords:
            chrom_name = item_id.split("_")[-1]  # goes to the first column of a bed file
            start = int(gene.location.start)
            end = int(gene.location.end)
            span_start = start - span_len - 1
            span_end = end + span_len
            # make a future table row
            row = [chrom_name, gene.id, start, end, span_start, span_end, int(gene.location.strand)]
            rg_collection.append(row)

        # a table with span and resistance gene coordinates
        rg_ranges = pd.DataFrame(
            columns=["chrom", "gene_id", "gene_start", "gene_end", "span_start", "span_end", "strand"],
            data=rg_collection)
        # add columns for ranges crossing ends of a chromosome
        rg_ranges["span_over_5_start"] = np.where(rg_ranges["span_start"] < 0, rg_ranges["span_start"] + dna_len,
                                                  np.nan)
        rg_ranges["span_over_5_end"] = np.where(np.isnan(rg_ranges["span_over_5_start"]), np.nan, dna_len)
        rg_ranges["span_over_3_end"] = np.where(rg_ranges["span_end"] > dna_len, rg_ranges["span_end"] - dna_len,
                                                np.nan)
        rg_ranges["span_over_3_start"] = np.where(np.isnan(rg_ranges["span_over_3_end"]), np.nan, 0)

        # making a bed file for all ranges
        bed_dataframes_list = make_bed(rg_ranges, score=0)
        # change float to int for BEDtools
        for bed_df in bed_dataframes_list:
            bed_df["range_start"] = bed_df.range_start.astype("int64")
            bed_df["range_end"] = bed_df.range_end.astype("int64")

    else:
        empty_columns = ["chrom", "range_start", "range_end", "name", "score", "strand"]
        bed_dataframes_list = [pd.DataFrame(columns=empty_columns)]

    # the output bed data frame contains negative and overly positive coordinates
    return bed_dataframes_list, msg_len + msg_count


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
strain = in_assembly.split("/")[2]
in_gff = os.path.join(snakemake.input[1], strain + "_genomic.gff")
in_rgi = snakemake.input[2]
range_len = int(snakemake.params[0])
min_plasmid_size = int(snakemake.params[1])

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
    # gff_ids = [rec.id.split("_")[-1] for rec in gff]
    # read joined assembly file as dict
    # assembly = SeqIO.to_dict(SeqIO.parse(in_assembly, "fasta"))
    # filter it because not all of the records present in GFF
    # but you can not filter by ID because IDs coming from SPAdes (if it finished successfully) are different from
    # IDs coming from Unicycler
    # filter assembly using lengths
    assembly_filtered = list()
    for assembly_item in SeqIO.parse(in_assembly, "fasta"):
        for gff_item in gff:
            if len(assembly_item) == len(gff_item):
                assembly_filtered.append(assembly_item)

    # sort both assembly and gff in order to iteration below worked
    assembly_filtered = sorted(assembly_filtered, reverse=True, key=len)
    gff = sorted(gff, reverse=True, key=len)

    # iterate through chromosome and plasmids
    bed_lol = [list(), list(), list()]
    messages = list()
    for i in range(len(gff)):
        # i in assembly corresponds to i in gff
        record_len = len(assembly_filtered[i].seq)
        record_id = assembly_filtered[i].id

        ranges_bed_list, bed_message = make_bed_file_for_rg(gff_record=gff[i], rgi_dataframe=rgi_notLoose,
                                                            dna_len=record_len, span_len=range_len)
        # turn 5-end crossing ranges' starts to zeros and 3-end crossing ranges to chromosome length
        for x in range(len(ranges_bed_list)):
            ranges_bed_list[x]["range_start"] = np.where(ranges_bed_list[x]["range_start"] < 0,
                                                         0, ranges_bed_list[x]["range_start"])
            ranges_bed_list[x]["range_end"] = np.where(ranges_bed_list[x]["range_end"] > record_len,
                                                       record_len, ranges_bed_list[x]["range_end"])
        # add bed data frames to the corresponding lists
        for j in range(len(ranges_bed_list)):
            bed_lol[j].append(ranges_bed_list[j])

        messages.append(bed_message)
    # here you have three types of bed files as lists of data frames
    # join bed files for all chromosomes/plasmids in current GFF
    bed_df_list = [pd.concat(lst) for lst in bed_lol]

    # write dataframes to three BED files
    for i in range(3):
        bed_df_list[i].to_csv(snakemake.output[i], sep="\t", index=False, header=False)
    # write messages to log
    log.writelines(messages)

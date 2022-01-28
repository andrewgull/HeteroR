from Bio import SeqIO
import pandas as pd
from BCBio import GFF
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
    bed = pd.concat([bed_normal, bed_5, bed_3])
    return bed


def make_bed_file_for_rg(gff_record, rgi_dataframe, dna_len, span_len):
    """
    makes bed file for genomic ranges with resistance genes
    it takes into account the length of dna record and range
    :param gff_record: GFF-object, list of gff records of a particular assembly
    :param rgi_dataframe: pandas DataFrame, pd DataFrame of RGI results
    :param dna_len: integer, length of a current record (chromosome or plasmid)
    :param span_len: integer, length of up- and downstream region flanking a resistance gene
    :return: bed formatted pandas DataFrame, negative span pandas DataFrame, message with number of found genes
    """
    genes = [feature for feature in gff_record.features if feature.type == "gene"]  # here we have IDs and positions
    item_id = gff_record.id

    # find median gene length
    genes_lengths = [gene.location.end.position - gene.location.start.position for gene in genes]
    median_gene = pd.Series(genes_lengths).median()
    # adjust span_len according to the record length
    # span_len = (record_len - median_gene)/2, if 2*span_len > record_len
    record_length = len(gff[1].seq)
    if record_length < span_len*2:
        span_len = round((record_length - median_gene)/2)
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
        rg_ranges = pd.DataFrame(columns=["chrom", "gene_id", "gene_start", "gene_end", "span_start", "span_end", "strand"],
                                 data=rg_collection)
        # add columns for ranges crossing ends of a chromosome
        rg_ranges["span_over_5_start"] = np.where(rg_ranges["span_start"] < 0, rg_ranges["span_start"] + dna_len, np.nan)
        rg_ranges["span_over_5_end"] = np.where(np.isnan(rg_ranges["span_over_5_start"]), np.nan, dna_len)
        rg_ranges["span_over_3_end"] = np.where(rg_ranges["span_end"] > dna_len, rg_ranges["span_end"] - dna_len, np.nan)
        rg_ranges["span_over_3_start"] = np.where(np.isnan(rg_ranges["span_over_3_end"]), np.nan, 0)

        # making a bed file for all ranges
        bed_dataframe = make_bed(rg_ranges, score=0)
        # change float to int for BEDtools
        bed_dataframe["range_start"] = bed_dataframe.range_start.astype("int64")
        bed_dataframe["range_end"] = bed_dataframe.range_end.astype("int64")
    else:
        empty_columns = ["chrom", "range_start", "range_end", "name", "score", "strand"]
        bed_dataframe = pd.DataFrame(columns=empty_columns)
    # the output bed data frame must be further transformed
    return bed_dataframe, msg_len + msg_count


in_assembly = "results/assemblies_joined/DA63746/assembly.fasta"

in_gff = "results/annotations/DA63746/prokka/DA63746_genomic.gff"
in_rgi = "results/resistance_genes/DA63746/rgi_table.txt"

rgi = pd.read_csv(in_rgi, sep="\t")
rgi_notLoose = rgi[rgi["Cut_Off"] != "Loose"]

i = 0
min_plasmid_size = 1000
range_len = 100000

with open(in_gff) as f:
    gff = [rec for rec in GFF.parse(f) if len(rec.seq) >= min_plasmid_size]
gff_ids = [rec.id.split("_")[-1] for rec in gff]

assembly = SeqIO.to_dict(SeqIO.parse(in_assembly, "fasta"))
assembly_filtered = [assembly[key] for key in gff_ids]
strain = in_assembly.split("/")[2]

bed_list = list()
messages = list()

for i in range(len(gff)):
    record_len = len(assembly_filtered[i].seq)
    record_id = assembly_filtered[i].id

    ranges_bed, bed_message = make_bed_file_for_rg(gff_record=gff[i], rgi_dataframe=rgi_notLoose,
                                                   dna_len=record_len, span_len=range_len)
    # turn negative range starts to zeros and 3-end crossing ranges to chromosome length
    ranges_bed["range_start"] = np.where(ranges_bed["range_start"] < 0, 0, ranges_bed["range_start"])
    ranges_bed["range_end"] = np.where(ranges_bed["range_end"] > record_len, record_len, ranges_bed["range_end"])
    bed_list.append(ranges_bed)
    messages.append(bed_message)

joined_bed_dataframe = pd.concat(bed_list)
joined_bed_dataframe.to_csv("test_bed_output.tsv", sep="\t", index=False, header=False)


# to join sequences from two files with the same IDs use seqkit concat
# but you will need to split records with repeated IDs to separate files :(

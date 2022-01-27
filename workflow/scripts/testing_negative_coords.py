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
    # the output bed data frame must be further transformed
    return bed_dataframe, msg


in_assembly = "results/assemblies_joined/DA63746/assembly.fasta"

in_gff = "results/annotations/DA63746/prokka/DA63746_genomic.gff"
in_rgi = "results/resistance_genes/DA63746/rgi_table.txt"

rgi = pd.read_csv(in_rgi, sep="\t")
rgi_notLoose = rgi[rgi["Cut_Off"] != "Loose"]

i = 0
min_plasmid_size = 1000
circ = False
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
                                                                    dna_len=record_len, span_len=range_len,
                                                                    circular=circ)
    # turn negative range starts to zeros and 3-end crossing ranges to chromosome length
    ranges_bed["range_start"] = np.where(ranges_bed["range_start"] < 0, 0, ranges_bed["range_start"])
    ranges_bed["range_end"] = np.where(ranges_bed["range_end"] > record_len, record_len, ranges_bed["range_end"])
    bed_list.append(ranges_bed)
    messages.append(bed_message)

joined_bed_dataframe = pd.concat(bed_list)
joined_bed_dataframe.to_csv("test_bed_output.tsv", sep="\t", index=False, header=False)

# look at spans crossing left end of the DNA
neg_coords_left = all_coords[all_coords["span_start"] < 0]
neg_coords_left.head()  # there's just one line

# look at spans crossing right end of the DNA
neg_coords_right = all_coords[all_coords["span_end"] > record_len]  # there are 0 such lines
# add one
df_add = pd.DataFrame({"chrom": ["1"], "gene_id": ["GENE"], "gene_start": [4900000], "gene_end": [4901353],
              "span_start": [4800000], "span_end": [5001353], "strand": [1]})
all_coords = all_coords.append(df_add)
# add two columns: crossing_left_end & crossing_right_end
all_coords["span_over_5_start"] = np.where(all_coords["span_start"] < 0, all_coords["span_start"] + record_len, np.nan)
all_coords["span_over_5_end"] = np.where(np.isnan(all_coords["span_over_5_start"]), np.nan, record_len)
all_coords["span_over_3_end"] = np.where(all_coords["span_end"] > record_len, all_coords["span_end"] - record_len, np.nan)
all_coords["span_over_3_start"] = np.where(np.isnan(all_coords["span_over_3_end"]), np.nan, 0)

# you can make a single bed file (i.e. DataFrame) for everything

# to join sequences from two files with the same IDs use seqkit concat

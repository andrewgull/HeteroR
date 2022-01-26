from Bio import SeqIO
import pandas as pd
from BCBio import GFF
import numpy as np


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
        span_start = start - span_len - 1
        span_end = end + span_len
        # check left end
        # if start < span_len:
        #     # to include leftmost letter subtract 1
        #     span_start = 0
        # else:
        #     span_start = start - span_len - 1
        # # check right end: bedtools includes it
        # if end + span_len > dna_len:
        #     span_end = dna_len
        # else:
        #     span_end = end + span_len
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
    # TODO: make two bed files for rg_ranges_neg: one with up to oriC and the other for trans-oriC coords, see pseudocode.txt on Argos

    return bed_dataframe, rg_ranges, msg


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

# for i in range(len(gff)):
record_len = len(assembly_filtered[i].seq)
record_id = assembly_filtered[i].id

ranges_bed, all_coords, bed_message = make_bed_file_for_rg(gff_record=gff[i], rgi_dataframe=rgi_notLoose,
                                                                dna_len=record_len, span_len=range_len,
                                                                circular=circ)
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
all_coords["span_over_5_end"] = np.where(all_coords["span_start"] < 0, all_coords["span_start"] + record_len, np.nan)
all_coords["span_over_3_end"] = np.where(all_coords["span_end"] > record_len, all_coords["span_end"] - record_len, np.nan)

#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
import pandas as pd
from BCBio import GFF
import os


def parse_spacer(header):
    """
    :param header: a char string like '>1:0-190543:951:190482:14m'
    :return: a list of values from parsed header
    """
    first_split = header.split(":")
    record_id = first_split[0][1:]
    repeat_len = int(first_split[-1][:-1])
    range_ = first_split[1].split("-")
    range_start, range_end = int(range_[0]), int(range_[-1])
    repeat_1_start_in_range = int(first_split[2])
    repeat_2_end_in_range = int(first_split[3])

    repeat_1_start_in_chrom = range_start + repeat_1_start_in_range
    repeat_1_end_in_chrom = repeat_1_start_in_chrom + repeat_len

    repeat_2_end_in_chrom = range_start + repeat_2_end_in_range
    repeat_2_start_in_chrom = repeat_2_end_in_chrom - repeat_len

    return [record_id, repeat_1_start_in_chrom, repeat_1_end_in_chrom, repeat_2_start_in_chrom,
            repeat_2_end_in_chrom, repeat_len]


def gff_object(features_df, record, strand=1, feature_type="direct_repeat"):
    """
    :param features_df: a pd.DataFrame with features for gff object
    :param record: SeqRecord object from assembly
    :param strand: strand for GFF files, 1 or -1
    :param feature_type: name of the feature
    :return: SeqRecord with features for GFF.write()
    """
    # take features found in provided record
    features_in_record = features_df[features_df.record_id == record.id]
    # init repeat IDs
    repeat_id = 0
    for row in features_in_record.iterrows():
        repeat_id += 1
        # use one ID for each pair of repeats (1,2,3,4 elements of each row)
        record.features.append(SeqFeature(FeatureLocation(ExactPosition(row[1][1]),
                                                          ExactPosition(row[1][2]), strand=strand),
                                          type=feature_type, id=str(repeat_id)))
        record.features.append(SeqFeature(FeatureLocation(ExactPosition(row[1][3]),
                                                          ExactPosition(row[1][4]), strand=strand),
                                          type=feature_type, id=str(repeat_id)))

    return record


# DEFINE INPUTS AND OUTPUTS
# in_spacers = "/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_perfect_repeats_GRF_test/perfect.spacer.id"
# in_assembly = "/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_assembly.fasta"
# out_gff = "/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_repeats.gff"
grf_out_filename = "perfect.spacer.id"  # IT CAN BE CHANGED
in_spacers = os.path.join(snakemake.input[0], grf_out_filename)
in_assembly = snakemake.input[1]
out_gff = snakemake.output[0]
# set minimal repeat length
min_repeat_length = int(snakemake.params[0])

# read spacer IDs
with open(in_spacers) as f:
    spacer_ids = [line.rstrip() for line in f.readlines()]

# read assembly
assembly = [rec for rec in SeqIO.parse(in_assembly, "fasta")]

# make features rows from spacer IDs
gff_rows = [parse_spacer(line) for line in spacer_ids]
gff_df = pd.DataFrame(columns=["record_id", "start_1", "end_1", "start_2", "end_2", "length"], data=gff_rows)

# filter out too short repeats
gff_df = gff_df[gff_df.length > min_repeat_length]

# make a gff object from this filtered data frame
# one SeqRecord with features per record in assembly
gff_records = [gff_object(gff_df, record) for record in assembly]

# write to a GFF file
with open(out_gff, 'w') as out:
    GFF.write(recs=gff_records, out_handle=out, include_fasta=False)

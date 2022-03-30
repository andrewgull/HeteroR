#!/usr/bin/env python3
"""
parses GRF output (grf-main) and makes GFF file for visualization of repeats (TDR) in a genomic browser
for each strain it runs two times: perfect.spacer.id and imperfect.id (for repeats with no mismatches)
Second file in each run corresponds to
type="imperfect" in parse_grf_output()
"""
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
import pandas as pd
from BCBio import GFF
import os


def parse_grf_output_w_ranges(header):
    """
    MUST BE RE-WRITTEN (OR NOT) FOR A CASE WHEN A GENE NAME IS IN THE HEADER
    :param header: a char string like '>1:0-190543:951:190482:14m' or '>1:155171-355997:1372:198472:10m1D3m'
    :return: a list of values from parsed header
    """
    first_split = header.split(":")
    record_id = first_split[0][1:]
    # this works for both mismatch or no mismatch types of output
    repeat_len = int(first_split[-1].split("m")[0])
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


def parse_grf_output_no_ranges(header):
    """
    for newer version of grf files: '>IPFHMEHC_00036_gene:12520:113509:13m1I2m2I5m' or
    '>IPFHMEHC_00036_gene:95122:101284:29m'
    where headers do not contain coordinates of regions on chromosome (i.e. 'ranges')
    :param header: a string from grf output, see above
    :return: a list of coordinates from this header ready to be turned into a gff-record
    """
    first_split = header.split(":")
    record_id = first_split[0][1:]
    repeat_len = int(first_split[-1].split("m")[0])
    repeat_1_start_in_range = int(first_split[1])
    repeat_1_end_in_range = repeat_1_start_in_range + repeat_len
    repeat_2_end_in_range = int(first_split[2])
    repeat_2_start_in_range = repeat_2_end_in_range - repeat_len
    return [record_id, repeat_1_start_in_range, repeat_1_end_in_range, repeat_2_start_in_range, repeat_2_end_in_range,
            repeat_len]


def gff_object(features_df, record, strand=1, feature_type="direct_repeat"):
    """
    :param features_df: a pd.DataFrame with features for gff object
    :param record: SeqRecord object from assembly
    :param strand: strand for GFF files, 1 or -1
    :param feature_type: name of the feature
    :return: SeqRecord with features for GFF.write()
    """
    # take features found in provided record
    # record_id in features_df contains "_gene"
    features_in_record = features_df[features_df.record_id == record.id + "_gene"]
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


def main(input_assembly, input_grf, min_len):
    # read spacer IDs
    with open(input_grf) as f:
        repeat_ids = [line.rstrip() for line in f.readlines()]

    # read assembly - fasta with annotated proteins/genes
    assembly = [rec for rec in SeqIO.parse(input_assembly, "fasta")]

    # make features rows from spacer IDs
    gff_rows = [parse_grf_output_no_ranges(line) for line in repeat_ids]
    gff_df = pd.DataFrame(columns=["record_id", "start_1", "end_1", "start_2", "end_2", "length"], data=gff_rows)
    gff_df.drop_duplicates(inplace=True)

    # filter out too short repeats
    gff_df = gff_df[gff_df.length >= min_len]

    # make a gff object from this filtered data frame
    # one SeqRecord with features per record in assembly
    gff_records = [gff_object(gff_df, record) for record in assembly]
    return gff_records


if __name__ == '__main__':
    # DEFINE INPUTS AND OUTPUTS
    # in_spacers = "/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_perfect_repeats_GRF_test/perfect.spacer.id"
    # in_assembly = "/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_assembly.fasta"
    # out_gff = "/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_repeats.gff"
    # grf_out_filename = "perfect.spacer.id"  # IT CAN BE CHANGED
    in_perfect = os.path.join(snakemake.input[0], "perfect.spacer.id")
    in_imperfect = os.path.join(snakemake.input[0], "imperfect.id")
    in_assembly = snakemake.input[1]
    out_gff_perfect = snakemake.output[0]
    out_gff_imperfect = snakemake.output[1]
    # set minimal repeat length
    min_repeat_length = int(snakemake.params[0])

    # create gff records from perfect.spacer.id
    gff_records_perfect = main(in_assembly, in_perfect, min_len=min_repeat_length)
    # remove records (genes) with 0 repeats
    gff_records_perfect = [rec for rec in gff_records_perfect if len(rec.features) > 0]
    # write to a GFF file
    with open(out_gff_perfect, 'w') as out:
        GFF.write(recs=gff_records_perfect, out_handle=out, include_fasta=False)

    # create gff records from imperfect.id
    gff_records_imperfect = main(in_assembly, in_imperfect, min_len=min_repeat_length)
    # remove records (genes) with 0 repeats
    gff_records_imperfect = [rec for rec in gff_records_imperfect if len(rec.features) > 0]
    # write them to a file
    with open(out_gff_imperfect, 'w') as out:
        GFF.write(recs=gff_records_imperfect, out_handle=out, include_fasta=False)

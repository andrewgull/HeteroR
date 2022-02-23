# script for making the feature table - should it be SQLite?
# list of features:
# 1. number of RG
# 2. number of RG on plasmids
# 3. nearest distance to oriC
# 4. median distance to oriC
# 5. total number of repeats
# 6. median number of repeats
# 7. median repeat length
# 8. longest repeat length
# 9. match concetration (=number of mismatches)
# 10. amplifiable region (AR) min length
# 11. AR median length
# 12. BLACK BOX

import pandas as pd

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


def make_repeat_df(input_grf, min_len):
    # read spacer IDs
    with open(input_grf) as f:
        repeat_ids = [line.rstrip() for line in f.readlines()]

    # make features rows from spacer IDs
    gff_rows = [parse_grf_output_no_ranges(line) for line in repeat_ids]
    gff_df = pd.DataFrame(columns=["record_id", "start_1", "end_1", "start_2", "end_2", "length"], data=gff_rows)
    gff_df.drop_duplicates(inplace=True)

    # filter out too short repeats
    gff_df = gff_df[gff_df.length > min_len]

    return gff_df


strains = ["DA63084", "DA63186", "DA63322", "DA63946", "DA64026"]
antibiotics = []

path_to_perfect = "/home/andrei/Data/HeteroR/results/direct_repeats/%s/repeats_no_mismatch/perfect.spacer.id"
path_to_rgi = "/home/andrei/Data/HeteroR/results/resistance_genes/%s/rgi_table.txt" % strains[0]

df_lst = list()
for strain in strains:
    file_path = path_to_perfect %strain
    df = make_repeat_df(input_grf=file_path, min_len=20)
    df["strain"] = strain
    df_lst.append(df)

repeat_df = pd.concat(df_lst)
repeat_df.head()

rgi_df = pd.read_csv(path_to_rgi, delimiter="\t")

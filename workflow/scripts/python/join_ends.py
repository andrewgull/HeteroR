from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys


def join_ends(side_dict, normal_dict, left=True):
    """
    joins seqs in 'left+normal' way or in 'normal+right' way
    :param side_dict: dictionary of SeqRecords overlapping chromosomal ends; should not be empty
    :param normal_dict: dictionary of SeqRecords non overlapping chromosomal ends
    :param left: boolean, if True then side_dict has sequences that came from left; else side_dict contains
    sequences from right
    :return: a list of joined SeqRecords
    """
    joined_ends_list = list()
    for key in side_dict.keys():
        if left:
            # start joining from the side dict
            joined_seq = side_dict[key].seq + normal_dict[key].seq
        else:
            joined_seq = normal_dict[key].seq + side_dict[key].seq
        # get IDs and make a SeqRecord
        joined_id = side_dict[key].id
        joined_record = SeqRecord(seq=joined_seq, id=joined_id, name=joined_id, description=joined_id)
        joined_ends_list.append(joined_record)
    return joined_ends_list

# input files
input_normal = snakemake.input[0]
input_left = snakemake.input[1]
input_right = snakemake.input[2]

# output files
output_normal = snakemake.output[0]
output_5_end = snakemake.output[1]
output_3_end = snakemake.output[2]

# open log file
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    # read sequences
    normal_records = SeqIO.to_dict(SeqIO.parse(input_normal, "fasta"))
    left_records = SeqIO.to_dict(SeqIO.parse(input_left, "fasta"))
    right_records = SeqIO.to_dict(SeqIO.parse(input_right, "fasta"))

    # CASE 1 (the most probable): normal is non-empty, left and right are empty
    if len(normal_records) > 0 & len(left_records) == 0 & len(right_records) == 0:
        # write essentially the same file, snakemake will detect these outputs
        joined_5_end = list()
        joined_3_end = list()
        # NORMAL STAYS THE SAME
    # CASE 2 (less probable): normal and left are non-empty, right is empty
    elif len(normal_records) > 0 & len(left_records) > 0 & len(right_records) == 0:
        joined_5_end = join_ends(left_records, normal_records, left=True)
        joined_3_end = list()
        for key in left_records.keys():
            dropped = normal_records.drop(key, None)
    # CASE 3 (same probability): normal and right are non-empty, left is empty
    elif len(normal_records) > 0 & len(right_records) > 0 & len(left_records) == 0:
        joined_3_end = join_ends(right_records, normal_records, left=False)
        joined_5_end = list()
        for key in right_records.keys():
            dropped = normal_records.drop(key, None)
    # CASE 4 (unlikely to happen): all three dicts are non-empty
    elif len(normal_records) > 0 & len(right_records) > 0 & len(left_records) > 0:
        joined_5_end = join_ends(left_records, normal_records, left=True)
        joined_3_end = join_ends(right_records, normal_records, left=False)
        keys_to_drop = set(list(left_records.keys()) + list(right_records.keys()))
        for key in keys_to_drop:
            dropped = normal_records.drop(key, None)
    else:  # should not happen
        joined_5_end = list()
        joined_3_end = list()

    # write outputs
    SeqIO.write(list(normal_records.values()), output_normal, "fasta")
    SeqIO.write(joined_5_end, output_5_end, "fasta")
    SeqIO.write(joined_3_end, output_3_end, "fasta")

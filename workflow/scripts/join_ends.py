from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


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


# in test directory
# /home/andrei/Data/HeteroR/test_dir/bed_seqkit
# input files
input_normal = "regions_normal.fasta"
input_left = "regions_left.fasta"
input_right = "regions_right.fasta"

# output files
output_5_end = ""
output_3_end = ""
output_normal = ""

# read sequences
normal_records = SeqIO.to_dict(SeqIO.parse(input_normal, "fasta"))
left_records = SeqIO.to_dict(SeqIO.parse(input_left, "fasta"))
right_records = SeqIO.to_dict(SeqIO.parse(input_right, "fasta"))

# join overlapping 5-end & ranges within a chromosome
joined_5_end = join_ends(left_records, normal_records, left=True)

# do the same with 3-end
joined_3_end = join_ends(right_records, normal_records, left=False)

# collect the keys from both left and right dicts and remove them from normal
keys_to_drop = set(list(left_records.keys()) + list(right_records.keys()))
for key in keys_to_drop:
    dropped = normal_records.drop(key, None)

# write new normal_records (without dropped keys of course)
SeqIO.write(joined_5_end, output_5_end, "fasta")
SeqIO.write(joined_3_end, output_3_end, "fasta")
SeqIO.write(normal_records, output_normal, "fasta")

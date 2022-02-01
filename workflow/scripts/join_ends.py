from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# in test directory
# /home/andrei/Data/HeteroR/test_dir/bed_seqkit
input_normal = "regions_normal.fasta"
input_left = "regions_left.fasta"
input_right = "regions_right.fasta"

normal = SeqIO.to_dict(SeqIO.parse(input_normal, "fasta"))
left = SeqIO.to_dict(SeqIO.parse(input_left, "fasta"))
right = SeqIO.to_dict(SeqIO.parse(input_right, "fasta"))

# join 5-end overlapping
joined_5_end = list()
for key in left.keys():
    joined_seq = left[key].seq + normal[key].seq
    joined_id = left[key].id
    joined_record = SeqRecord(seq=joined_seq, id=joined_id, name=joined_id, description=joined_id)
    joined_5_end.append(joined_record)
    # drop key from normal
    normal.pop(key, None)

# write joined_5_end
# do the same with 3-end
# write new normal (without dropped keys of course)

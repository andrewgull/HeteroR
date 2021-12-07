#!/usr/bin/env python3
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.Seq import Seq
import pandas as pd
from BCBio import GFF


def parse_spacer(header):
    # header '>1:0-190543:951:190482:14m'
    first_split = header.split(":")
    record_id = first_split[0][1:]
    repeat_len = int(first_split[-1][:-1])
    range_ = first_split[1].split("-")
    range_start, range_end = int(range_[0]), int(range_[-1])
    repeat_1_start_inRange = int(first_split[2])
    repeat_2_end_inRange = int(first_split[3])

    repeat_1_start_inChrom = range_start + repeat_1_start_inRange
    repeat_1_end_inChrom = repeat_1_start_inChrom + repeat_len

    repeat_2_end_inChrom = range_start + repeat_2_end_inRange
    repeat_2_start_inChrom = repeat_2_end_inChrom - repeat_len

    return [record_id, repeat_1_start_inChrom, repeat_1_end_inChrom, repeat_2_start_inChrom, repeat_2_end_inChrom]


min_repeat_length = 10
with open("/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_perfect_repeats_GRF_test/perfect.spacer.id") as f:
    spacer_ids = [line.rstrip() for line in f.readlines() if str(min_repeat_length) + "m" not in line]

assembly = [rec for rec in SeqIO.parse("/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_assembly.fasta", "fasta")]

# a line in lines looks like this: >1:0-190543:951:190482:14m'
# >1:0-190543 - fasta id
# 951:190482:14m - start of the 1st match, end of the 2nd match and length of the match

# read GBK to look inside
gbk = [rec for rec in SeqIO.parse("/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_genomic.gbk", "genbank")]

gbk[0].features[0]

in_rgi = "/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_rgi_table.txt"
rgi = pd.read_csv(in_rgi, sep="\t")

# parse gff or genebank
gff = [rec for rec in GFF.parse("/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_genomic.gff")]
# same structure inside
gff[0].features[0]

# you have info on genes
# see the copybook
# gff files look easier to create
# each record in assembly should have separate gff file with repeats
# create a SeqRecord (chromosome or plasmid) where each repeat is a feature with location
# to get location features parse and calculate coordinates from spacer IDs
strand = 1  # is always 1
n = 0
for line in spacer_ids[:20]:
    n += 1
    repeat_features = parse_spacer(line)
    # use 1st record from the assembly as an example
    assembly[0].features.append(SeqFeature(FeatureLocation(ExactPosition(repeat_features[1]),
                                                           ExactPosition(repeat_features[2]), strand=strand),
                                           type='direct_repeat', id=str(n)))
    assembly[0].features.append(SeqFeature(FeatureLocation(ExactPosition(repeat_features[3]),
                                                           ExactPosition(repeat_features[4]), strand=strand),
                                           type='direct_repeat', id=str(n)))
with open("/home/andrei/Data/HeteroR/test_dir/GRF/DA62886_repeats.gff", 'w') as out_gff:
    GFF.write(recs=[assembly[0]], out_handle=out_gff, include_fasta=False)

# TODO: expand to all records in assembly

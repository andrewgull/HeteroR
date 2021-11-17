#!/usr/bin/env python3

"""a script for genomes and repeats generation"""


import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse


def get_args():
    """
    Get command line arguments
    """

    parser = argparse.ArgumentParser(
        description="This scripts generates a random genome with insertions.\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-gl", "--genome_length", type=int, metavar="<integer>",
                        help="length of a genome to generate")
    parser.add_argument("-rl", "--repeat_length",  type=int, metavar="<integer>",
                        help="length of a repeat unit to insert")
    parser.add_argument("-p", "--positions", type=str, metavar="<pos1,pos2>",
                        help="positions to insert repeats to, coma separated, no whitespace")
    parser.add_argument("-o", "--output", type=str, metavar="<output file>",
                        help="output file, fasta format")

    return parser.parse_args()


def random_seq(length):
    return ''.join(random.choice('ACTG') for char in range(length))


def insert_repeats(dna, positions, repeat_unit):
    """
    param genome: DNA string
    param positions: a list with positions to insert repeats to
    param n_copies: number of times to insert a repeat
    param length: a repeat unit's length
    return: genome with inserted repeats
    """
    # create a random repeat
    # repeat = random_seq(rep_length)

    for pos in positions:
        dna = dna[:int(pos)] + repeat_unit + dna[int(pos):]

    return dna


args = get_args()

# parse positions string
positions_list = args.positions.split(",")

# create a repeat unit
repeat = random_seq(args.repeat_length)
# create a genomic string with inserts
genome = insert_repeats(dna=random_seq(args.genome_length), positions=positions_list, repeat_unit=repeat)

# make SeqRecords
genome_record = SeqRecord(Seq(genome), id="genomic_string", description="genomic_string_pos_%s" % "-".join(positions_list))
repeat_record = SeqRecord(Seq(repeat), id="repeat_unit", description="repeat_unit_length_%i" % len(repeat))
# write them to a file
SeqIO.write([genome_record, repeat_record], args.output, "fasta")

print("A DNA string with %i repeat units has been generated" % len(positions_list))

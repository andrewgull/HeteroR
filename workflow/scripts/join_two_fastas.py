#!/usr/bin/env python3
"""Joins Unicycler assembly/prokka annotation and SPAdes plasmid assembly/tRNAScan-SE annotation. """
import argparse
import os
import sys
from Bio import SeqIO


def get_args():
    """
    Get command line arguments
    requires directories not files as inputs because of the outputs produced by rule unicycler and rule plasmid_assembly
    """

    parser = argparse.ArgumentParser(
        description="Joins Unicycler assembly and SPAdes plasmid assembly.\n"
                    "If no plasmid assembly exists it just copies Unicycler assembly to a new file.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("first_file", metavar="<path to a dir with the 1st fasta file (e.g. unicycler assembly)>",
                        help="file with genome assembly produced by Unicycler")
    parser.add_argument("second_file", metavar="<path to a directory with the 2nd fasta file (e.g. plasmids)>",
                        help="file with plasmids assembly produced by SPAdes")
    parser.add_argument("output", metavar="<filename>",
                        help="new joined assembly file")
    parser.add_argument("mode", metavar="<assembly or annotation>", help="which mode of joining to use")

    return parser.parse_args()


def joiner(file1, file2):
    """
    Reads in provided files and joins them together as lists
    :param file1: 1st file to join (unicycler assembly or any other fasta)
    :param file2: 2nd file to join (plasmid assembly or any other fasta)
    :return: joined fasta as list
    """
    # read the 1st file
    fasta1 = [seq for seq in SeqIO.parse(file1, 'fasta')]
    # check if the 2nd file exists and its size is greater than 0
    # useful in 'assembly' mode when file2 = plasmid
    if os.path.isfile(file2) and os.path.getsize(file2) > 0:
        # file exists, read it and join
        fasta2 = [seq for seq in SeqIO.parse(file2, 'fasta')]
        fasta_joined = fasta1 + fasta2
        print("Joined assembly contains %i sequences" % len(fasta_joined))
    else:
        # return the 1st file
        # in 'assembly' mode it is an unicycler assembly
        fasta_joined = fasta1
        print("No 2nd file found (finished plasmid assembly doesn't exist)")

    return fasta_joined


# execute joiner() and send its stderr/stdout to the log file
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    # make proper filenames (not directories!)
    filename1 = snakemake.input[0] + "/assembly.fasta"
    filename2 = snakemake.input[1] + "/scaffolds.fasta"

    # elif snakemake.input[2] == "annotation":
    # It might not work because it's a mix of NUC and AA
    #     filename1 = snakemake.input[0] + "/annotation_prokka.faa"
    #     filename2 = snakemake.input[1]
    output_file = joiner(file1=filename1, file2=filename2)
    SeqIO.write(output_file, snakemake.output[0], 'fasta')

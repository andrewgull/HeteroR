#!/usr/bin/env python3
"""Joins Unicycler assembly/prokka annotation and SPAdes plasmid assembly/tRNAScan-SE annotation. """
import argparse
import os
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
    Joins provided files
    param: file1 - 1st file to join (unicycler assembly or any other fasta)
    param: file2 - 2nd file to join (plasmid assembly or any other fasta)
    """
    # read the 1st file
    fasta1 = [seq for seq in SeqIO.parse(file1, 'fasta')]
    # check if plasmid assembly exists and its size
    if os.path.isfile(file2) and os.path.getsize(file2) > 0:
        # file exists, read it and join
        fasta2 = [seq for seq in SeqIO.parse(file2, 'fasta')]
        fasta_joined = fasta1 + fasta2
    else:
        # return the 1st file (unicycler assembly)
        fasta_joined = fasta1
        # it will go to a log file
        print("no finished plasmid assembly file found")

    return fasta_joined


if __name__ == '__main__':
    args = get_args()
    if args.mode == "assembly":
        # make proper filenames (not directories!)
        filename1 = args.first_file + "/assembly.fasta"
        filename2 = args.second_file + "/scaffolds.fasta"

    elif args.mode == "annotation":
        # It might not work because it's a mix of NUC and AA
        filename1 = args.first_file + "/annotation_prokka.faa"
        filename2 = args.second_file
    output_file = joiner(file1=filename1, file2=filename2)
    SeqIO.write(output_file, args.output, 'fasta')

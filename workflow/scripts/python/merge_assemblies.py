#!/usr/bin/env python3
"""Joins Unicycler assembly and SPAdes plasmid assembly. """
import os
import sys
from Bio import SeqIO


def joiner(file1, file2):
    """
    Reads in provided files and joins them together as lists
    :param file1: 1st file to join (unicycler assembly or any other fasta)
    :param file2: 2nd file to join (plasmid assembly or any other fasta)
    :return: joined fasta as list
    """
    # read the 1st file
    fasta1 = list(SeqIO.parse(file1, 'fasta'))
    
    # check if the 2nd file exists and its size is greater than 0
    if file2 and os.path.isfile(file2) and os.path.getsize(file2) > 0:
        # file exists, read it and join
        fasta2 = list(SeqIO.parse(file2, 'fasta'))
        fasta_joined = fasta1 + fasta2
        print(f"Joined assembly contains {len(fasta_joined)} sequences")
    else:
        # return the 1st file
        fasta_joined = fasta1
        print("No 2nd file found (finished plasmid assembly doesn't exist or is empty)")

    return fasta_joined

def run_snakemake():
    # execute joiner() and send its stderr/stdout to the log file
    with open(snakemake.log[0], "w") as f:
        sys.stderr = sys.stdout = f
        # make proper filenames (not directories!)
        filename1 = os.path.join(snakemake.input[0], "assembly.fasta")
        filename2 = os.path.join(snakemake.input[1], "scaffolds.fasta")

        output_records = joiner(file1=filename1, file2=filename2)
        SeqIO.write(output_records, snakemake.output[0], 'fasta')

if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake()
    else:
        # For importing in tests
        pass

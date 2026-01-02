#!/usr/bin/env python3

"""
renames IDs in prokka-made GBK annotation files
so that these GBK will be compatible with assembly's IDs and genomic browsers
"""

from Bio import SeqIO
import os
import sys


def rename_gbk_records(records):
    for rec in records:
        # ID name pattern is like this: LOOOHHPE_1
        # I need number after the underscore - this is chromosome's name in assembly file
        new_id = rec.id.split("_")[-1]
        rec.id = new_id
    return records

def run_snakemake():
    with open(snakemake.log[0], "w") as f:
        sys.stderr = sys.stdout = f
        # input is a directory
        in_gbk_dir = snakemake.input[0]
        out_gbk = snakemake.output[0]
        filename = snakemake.params[0]

        # make a filename
        in_gbk = os.path.join(in_gbk_dir, filename)

        records = list(SeqIO.parse(in_gbk, "genbank"))
        renamed_records = rename_gbk_records(records)

        SeqIO.write(renamed_records, out_gbk, "genbank")

if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake()
    else:
        # For importing in tests
        pass

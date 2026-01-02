#!/usr/bin/env python3

"""
to covert rgi output table to bed/gff file
"""

import pandas as pd
from Bio import SeqIO
import os
import sys


def parse_rgi_output(rgi_file, filter_condition):
    rgi = pd.read_csv(rgi_file, delimiter="\t")
    # iterate through rows and collect genes with features
    rgi_not_loose = rgi[rgi.Cut_Off != filter_condition]
    rgi_not_loose_ids = [row[1][0].split(" ")[0] for row in rgi_not_loose.iterrows()]
    
    # keep more information to add to qualifiers
    # ORF_ID(0), SNP(12), Drug Class(14), AMR Gene Family(16)
    rgi_not_loose_dict = dict()
    for row in rgi_not_loose.iterrows():
        id = row[1][0].split(" ")[0]
        if id in rgi_not_loose_ids:
            rgi_not_loose_dict[id] = [row[1][12], row[1][14], row[1][16]]
    
    return rgi_not_loose_ids, rgi_not_loose_dict

def process_gbk_records(gbk_records, rgi_not_loose_ids, rgi_not_loose_dict):
    for record in gbk_records:
        # change record.id to compatible with assembly's records names
        record.id = record.id.split("_")[-1]
        rgi_from_gbk = [f for f in record.features if f.type == 'gene' and f.qualifiers.get('locus_tag', [None])[0] in rgi_not_loose_ids]
        # rewrite record's features
        record.features = rgi_from_gbk

    # remove records with no features (i.e. no resistance genes)
    rgi_gbk = [record for record in gbk_records if len(record.features) > 0]

    # add more info to qualifiers
    for record in rgi_gbk:
        for feature in record.features:
            locus_tag = feature.qualifiers['locus_tag'][0]
            qualifiers_to_add = rgi_not_loose_dict[locus_tag]
            feature.qualifiers['SNP'] = qualifiers_to_add[0]
            feature.qualifiers['drug_class'] = qualifiers_to_add[1]
            feature.qualifiers['AMR_gene_family'] = qualifiers_to_add[2]
            
    return rgi_gbk

def run_snakemake():
    # inputs
    in_rgi = snakemake.input[0]
    strain = in_rgi.split("/")[2]
    in_gbk = os.path.join(snakemake.input[1], strain + "_genomic.gbk")
    out_gbk = snakemake.output[0]
    filter_condition = snakemake.params[0]

    # open log file
    with open(snakemake.log[0], "w") as f:
        sys.stderr = sys.stdout = f

        rgi_not_loose_ids, rgi_not_loose_dict = parse_rgi_output(in_rgi, filter_condition)
        
        gbk_records = list(SeqIO.parse(in_gbk, "genbank"))
        rgi_gbk = process_gbk_records(gbk_records, rgi_not_loose_ids, rgi_not_loose_dict)

        # write to a file
        SeqIO.write(rgi_gbk, out_gbk, "genbank")

if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake()
    else:
        # For importing in tests
        pass

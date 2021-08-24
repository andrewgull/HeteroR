#!/usr/bin/env python3

import sys
import subprocess
import glob
import pandas as pd
from tqdm import tqdm
import argparse


help_message = "Script to summarize coverage of Nanopore reads across samples. Requires seqkit"
genome_length = 5131220

try:
    strains_file = sys.argv[1]
except IndexError:
    print(help_message)
    sys.exit()

with open(strains_file, 'r') as f:
    strains = [line.rstrip() for line in f.readlines()]

# collect Nanopore all
nanopore_stats = list()
for strain in tqdm(strains):
    try:
        file = glob.glob("data_raw/%s/Nanopore/%s_all.fastq.gz" % (strain, strain))[0]
        proc = subprocess.Popen("seqkit stats %s -T" % file, shell=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
        stats_string = out.decode("utf-8").split("\n")[1].split("\t")
    except IndexError:  # when there is no Nanopore files
        stats_string = ["%s" %strain, "NaN", "NaN", "NaN", "NaN", "NaN", "NaN", "NaN"]

    nanopore_stats.append(stats_string)

nanopore_df = pd.DataFrame.from_records(nanopore_stats, columns=["file", "format", "type", "num_seqs", "sum_len",
                                                                 "min_len", "avg_len", "max_len"])
nanopore_df = nanopore_df.astype({'sum_len': 'float'})
nanopore_df["coverage"] = nanopore_df["sum_len"] / genome_length
nanopore_df.to_csv(path_or_buf="nanopore_stats.tsv", sep="\t", index=False)
print("Coverage ~25x or less is sparse. Use Unicycler.")

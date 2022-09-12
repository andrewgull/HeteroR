#!/usr/bin/env python3
"""Calculate coverage of Nanopore reads"""
import subprocess
import glob
import pandas as pd
from tqdm import tqdm
import argparse
import sys


def get_args():
    """
    Get command line arguments
    """

    parser = argparse.ArgumentParser(
        description="Script for calculating and summarizing coverage of concatenated Nanopore reads across samples.\n"
                    "Requires seqkit - https://bioinf.shenwei.me/seqkit/.\n"
                    "Run this script from the project directory, not from ./scripts.\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("strains_file", metavar="<strains file>",
                        help="file with strain names to copy, one name per line")
    parser.add_argument("genome_length", metavar="<genome length>", help="Approx. genome length, bp (E.coli=5131220)")
    parser.add_argument("output", metavar="<output file>", help="output filename, tab-separated")
    parser.add_argument('raw', help='raw or filtered', nargs='?', choices=('raw', 'filtered'))
    return parser.parse_args()


def main(strains_file, genome_length, raw=True):
    """Calculate and summarize coverage"""
    with open(strains_file, 'r') as f:
        strain_names = [line.rstrip() for line in f.readlines()]

    # collect Nanopore all
    nanopore_stats = list()
    for strain in tqdm(strain_names):
        try:
            if raw:
                file = glob.glob("resources/data_raw/%s/Nanopore/%s_all.fastq.gz" % (strain, strain))[0]
            else:
                file = glob.glob("results/data_filtered/%s/Nanopore/%s_all.fastq.gz" % (strain, strain))[0]
            proc = subprocess.Popen("seqkit stats %s -T" % file, shell=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            stats_string = out.decode("utf-8").split("\n")[1].split("\t")
        except IndexError:  # when there is no Nanopore files
            stats_string = ["%s" % strain, "NaN", "NaN", "NaN", "NaN", "NaN", "NaN", "NaN"]

        nanopore_stats.append(stats_string)

    nanopore_df = pd.DataFrame.from_records(nanopore_stats, columns=["file", "format", "type", "num_seqs", "sum_len",
                                                                     "min_len", "avg_len", "max_len"])
    nanopore_df = nanopore_df.astype({"sum_len": "float"})
    nanopore_df["coverage"] = nanopore_df["sum_len"] / int(genome_length)
    return nanopore_df


if __name__ == '__main__':
    args = get_args()
    strains = args.strains_file
    length = args.genome_length
    output = args.output
    if args.raw == "raw":
        use_raw_files = True
    elif args.raw == "filtered":
        use_raw_files = False
    else:
        print("Error! Neither raw nor filtered files to use!")
        sys.exit(1)

    coverage_stats = main(genome_length=length, strains_file=strains, raw=use_raw_files)
    min_cov, max_cov, avg_cov = coverage_stats["coverage"].min(), coverage_stats["coverage"].max(), \
                                coverage_stats["coverage"].mean()
    coverage_stats.to_csv(path_or_buf=output, sep="\t", index=False)
    print("Batch coverage:\nmin = %f\navg = %f\nmax = %f\nCoverage ~25x or less is sparse, good for Unicycler.\n"
          "Now you can create config and run the pipeline" % (min_cov, avg_cov, max_cov))

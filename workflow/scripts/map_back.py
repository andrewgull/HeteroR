#!/usr/bin/env python3

import argparse
import glob
import subprocess
from tqdm import tqdm


def get_args():
    """
    Get command line arguments
    """

    parser = argparse.ArgumentParser(
        description="This script maps short reads back onto assembly\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("assembly", metavar="assembly/{strain}/assembly.fasta",
                        help="path to hybrid assembly file")
    parser.add_argument("read1", metavar="data_filtered/{strain}/Illumina/{strain}_1.fq.gz",
                        help="forward reads file")
    parser.add_argument("read2", metavar="data_filtered/{strain}/Illumina/{strain}_2.fq.gz",
                        help="reverse reads files")
    parser.add_argument("output", metavar="assemblies/{strain}/assembly.sam")
    parser.add_argument("threads", metavar="<no. threads>", help="number of threads to use")
    return parser.parse_args()


def main(ass, r1, r2, out_sam, t):
    """main function for mapping
    param: ass - assembly
    param: r1 - forward reads
    param: r2 - reverse reads
    param: out_sam - output sam
    param: t - threads
    """
    # INDEX THE GENOME, MAP READS
    subprocess.call(["bwa", "index", ass])  # creates files .sa, .amb, ann, .pac, .bwt, can not be specified as -o
    subprocess.call(["bwa", "mem", "-t", t, ass, r1, r2, ">", out_sam])
    # SAMTOOLS: VIEW, SORT, INDEX
    subprocess.call(["samtools", "view", "-b", "$DIR/assembly/hybrid_uni/assembly.sam > $DIR/assembly/hybrid_uni/assembly.bam"])

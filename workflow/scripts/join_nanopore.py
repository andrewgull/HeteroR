#!/usr/bin/env python
"""Concatenate nanopore reads. This functionality is included in prepare_files.py"""
import subprocess
import argparse


def get_args():
    """
    Get command line arguments
    """

    parser = argparse.ArgumentParser(
        description="Script for concatenating Nanopore reads into one file with prefix '_all'.\n"
                    "Its functionality included in './scripts/prepare_files.py'",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("file_in", metavar="<>",
                        help="Input file?")
    parser.add_argument("file_out", metavar="<>", help="Output file?")
    parser.add_argument("threads", metavar="<no. threads>", help="number of threads to use in pigz")
    return parser.parse_args()


def main():
    args = get_args()

    # the following line is required because pigz doesn't need output name,
    # but the rule join_nanopore needs it, and it must be *_fastq.gz
    true_output = args.file_out[:-3]

    # the following line is required because the rule join_nanopore requires input,
    # but number of input files is always unknown
    true_input = args.file_in[-9:]  # this should be ".fastq.gz"
    threads = args.threads
    # "zcat {input} > {output} && pigz -p {threads} {output}"
    subprocess.call(["zcat", "*%s" % true_input, ">", "%s" % true_output], shell=True)
    subprocess.call(["pigz", "-p", "%s" % threads, "%s" % true_output], shell=True)


if __name__ == '__main__':
    main()

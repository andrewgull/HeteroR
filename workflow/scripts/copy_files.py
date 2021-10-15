#!/usr/bin/env python3
"""Copy read files from ARGOS RAW to ./data_row"""
import subprocess
from tqdm import tqdm
import argparse


def get_args():
    """
    Get command line arguments
    """

    parser = argparse.ArgumentParser(
        description="Script for copying read files from ARGOS directory to data_raw.\n"
                    "Run this script from the project directory, not from ./scripts.\n"
                    "ARGOS DIR='/home/andrei/Data/Argos/imb_sal_raw/500 Sepsis Eco/Sequencing/Strains'\n"
                    "ARGOS directory MUST be mounted!",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("strains_file", metavar="<strains file>",
                        help="file with strain names to copy, one name per line")

    return parser.parse_args()


def main(argos_path="/home/andrei/Data/Argos/imb_sal_raw/500\ Sepsis\ Eco/Sequencing/Strains"):
    """Main function to copy files"""
    args = get_args()
    strain_file = args.strains_file

    with open(strain_file, 'r') as f:
        strains = [line.rstrip() for line in f.readlines()]

    for strain in tqdm(strains):
        subprocess.call("cp -r %s/%s ./data_raw" % (argos_path, strain), shell=True)


if __name__ == '__main__':
    main()
    print("Done! Now you can run 'prepare_files.py'")

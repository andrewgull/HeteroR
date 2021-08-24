#!/usr/bin/env python3
import subprocess
from tqdm import tqdm
import argparse


def get_args():
    """
    Get command line arguments
    """

    parser = argparse.ArgumentParser(
        description="Script for copying read files from mounted ARGOS directory to data_raw.\n"
                    "ARGOS DIR='/mnt/imb_sal_raw/500 Sepsis Eco/Sequencing/Strains'\n"
                    "ARGOS directory MUST be mounted!",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("strains_file", metavar="filename", help="file with strain names to copy, one name per line")

    return parser.parse_args()


def main(argos_path="/mnt/imb_sal_raw/500\ Sepsis\ Eco/Sequencing/Strains"):
    args = get_args()

    strain_file = args.strains_file

    # argos directory must be mounted
    # argos_path = "/mnt/imb_sal_raw/500\ Sepsis\ Eco/Sequencing/Strains"

    with open(strain_file, 'r') as f:
        strains = [line.rstrip() for line in f.readlines()]

    # go to data_raw, you're in ./scripts
    # os.chdir("./data_raw")

    for strain in tqdm(strains):
        subprocess.call("cp -r %s/%s ./data_raw" % (argos_path, strain), shell=True)


if __name__ == '__main__':
    main()

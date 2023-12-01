# copy mutants to the correct dirs and rename them accordingly
# the correct placement is: "results/data_filtered/{parent}/Illumina/mutants/{parent}m_1.fq.gz" 
# where parent is name of a parental strain

# TOC
# read strains table
# turn it into a dict(parent: mutant files path)
# using this dict copy mutants' reads to the new location

import argparse
import shutil
import os
import glob
from tqdm import tqdm
import pandas as pd


def get_args():
    """
    print and parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Script for copying mutants files (short read sequencing) from Argos to the project's directory.\n"
                    "Run this script from the project directory, not from ./scripts.\n"
                    "All paths are hard-coded in the script.\n",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("parent_mutant_table", metavar="<file path>",
                        help="CSV file with parent-mutant correspondence")
    return parser.parse_args()


def table2dict(file):
    """
    param file: path to the file with parents and mutants names
    return: a dictionary {parent: mutant}
    """
    df = pd.read_csv(file, dtype={'strain': 'string',
                                  'mutant': 'string',
                                  'contamination': 'string'})
    df = df.dropna(subset=['mutant'])
    df = df.query('contamination != "cont"')
    df_dict = df.set_index('strain')['mutant'].to_dict()

    return df_dict


def create_paths(df_dict, strain, path):
    """
    param df_dict: dict from a data frame, parent: mutant
    param strain: strain name to recover form the dict
    param path: path to mutants reads files
    return: a list of paths to read1 and read2
    """
    # mut_path = path + mutant
    mut_path = path + df_dict[strain]
    reads = glob.glob(mut_path + "/*.fq.gz")
    return reads


def file_path2dir_path(path):
    # make a pth to directory where the mutant reads file shall be
    # make it from the full path to the file itself
    new_path = "/".join(path.split("/")[:-1])
    return new_path


if __name__ == "__main__":
    # arguments and constants
    args = get_args()
    # "/home/andrei/Data/HeteroR/resources/strain_lists/mutants_strains.csv"
    mut_strains = args.parent_mutant_table
    MUT_SOURCE = "/mnt/data/andrei/Data/Argos/imb_sal_raw/BGI_seqs/F22FTSEUHT1288-02_BACdvnwR_25Apr2023/soapnuke/clean/"

    # get {parent: mutant} dictionary
    strains_dict = table2dict(mut_strains)
    # create paths for each mutant of each parent
    for parent in tqdm(strains_dict.keys()):
        dest_path = f"/home/andrei/Data/HeteroR/results/data_filtered/{parent}/Illumina/mutants/{parent}"
        dest_path_fq1 = dest_path + "m_1.fq.gz"
        dest_path_fq2 = dest_path + "m_2.fq.gz"
        mutants_paths = create_paths(strains_dict, parent, MUT_SOURCE)

        if not os.path.exists(dest_path_fq1):
            os.makedirs(file_path2dir_path(dest_path_fq1), exist_ok=True)
            shutil.copyfile(mutants_paths[0], dest_path_fq1)
        if not os.path.exists(dest_path_fq2):
            os.makedirs(file_path2dir_path(dest_path_fq2), exist_ok=True)
            shutil.copyfile(mutants_paths[1], dest_path_fq2)

    print("Script finished. No errors.")

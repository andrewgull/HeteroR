#!/usr/bin/env python3
import os
import sys
import subprocess
from tqdm import tqdm


help_message = 'Script for linking read files, renaming them properly and joining Nanopore reads for each strain.\n' \
               'Usage: prepare_files.py <strains_file> <threads>'

try:
    # file with strains should be specified here
    strain_file = sys.argv[1]
except IndexError:
    print(help_message)
    sys.exit()

# argos directory must be mounted
argos_path = "/mnt/imb_sal_raw/500\ Sepsis\ Eco/Sequencing/Strains"

with open(strain_file, 'r') as f:
    strains = [line.rstrip() for line in f.readlines()]

# go to data_raw, you're in ./scripts
# os.chdir("./data_raw")

for strain in tqdm(strains):
    subprocess.call("cp -r %s/%s ./data_raw" % (argos_path, strain), shell=True)

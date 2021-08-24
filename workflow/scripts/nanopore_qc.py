#!/usr/bin/env python3
# special script for Nanopore QC run
# designed specifically for use with snakemake

import sys
import glob
import subprocess


dir = sys.argv[1]
cores = sys.argv[2]
data_path = "/home/andrei/Data/HRData"

# collect files with reads
read_files = glob.glob("%s/%s/Nanopore/*.gz" %(data_path, dir))

# run fastqc
for file in read_files:
    subprocess.run(["fastqc", "-t", "%s" %cores, "-o", "%s/qualcheck/Nanopore" %data_path, "%s" %file])

# fastqc -t 10 -q -o qualcheck/Nanopore DA*/Nanopore/*gz

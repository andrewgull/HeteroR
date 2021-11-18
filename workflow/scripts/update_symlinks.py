#!/usr/bin/env python3

"""
to fix broken links pointing to Illumina files
"""

import os
import glob
from tqdm import tqdm


def new_name(filepath, strain_name):

    if "_1.fq.gz" in filepath:
        suffix = "_1.fq.gz"
    else:
        suffix = "_2.fq.gz"

    link_name = os.path.join(os.getcwd(), "/".join(filepath.split("/")[:-1]) + "/renamed/%s" % strain_name + suffix)

    return link_name


strains = glob.glob("resources/data_raw/DA*")

for strain in tqdm(strains):
    read_files = glob.glob("%s/**/*.gz" % strain, recursive=True)
    broken_links = [file for file in read_files if os.path.islink(file) and not os.path.exists(file)]
    illumina_reads = [file for file in read_files if "Illumina" in file and "renamed" not in file]
    illumina_reads.sort()
    # remove broken links
    for link in broken_links:
        os.remove(link)
    # create new links: illumina -> renamed
    for file in illumina_reads:
        destination = new_name(file, os.path.basename(strain))
        source = os.path.join(os.getcwd(), file)
        try:
            os.symlink(source, destination)
        except FileExistsError:
            print("Illumina links exit in %s" % os.path.basename(strain))

print("Done")

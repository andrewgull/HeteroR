#!/usr/bin/env python3
"""
This script creates soft links to read files
Give them proper names!
and joins a strain's nanopore reads into a single file (former functionality of the rule join_nanopore)
naming scheme is as follows
Illumina: STRAIN_[1, 2].fastq.gz
Nanopore: STRAIN_all.fastq.gz
It is assumed that you've copied files from ARGOS to ../data_raw
"""
import os
import glob
import subprocess
from tqdm import tqdm
import argparse


def get_args():
    """
    Get command line arguments
    """

    parser = argparse.ArgumentParser(
        description="This script creates properly named soft links in ./data_raw\n"
                    "And concatenates Nanopore reads into one file prefixed with '_all'\n"
                    "Naming scheme is as follows:\n'STRAIN_[1,2].fq.gz' for Illumina reads\n"
                    "'STRAIN_all.fastq.gz' for Nanopore reads\n"
                    "It is assumed that all read files have been copied to ./data_raw",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("strains", metavar="<strains file>",
                        help="File with strain names, one name per line")
    parser.add_argument("threads", metavar="<no. threads>", help="number of threads to use in pigz")
    return parser.parse_args()


def new_name(path):
    """
    converts path to a filename
    param path: unix path-like string: "data_raw/strain/machine/barcode/filename.fq.gz"
    return: new filename like "strain_14.fq.gz" or "strain_1.fq.gz"
    """
    dirs = path.split("/")
    strain_name, old_name = dirs[1], dirs[-1]
    # suffix is universal for both Illumina and Nanopore files
    # it contains read number and extensions
    suffix = old_name.split("_")[-1]
    file_name = "_".join([strain_name, suffix])
    return file_name


def main():
    # PARAMETERS
    args = get_args()
    strains_file = args.strains
    threads = args.threads

    # 1. COMPRESSING FILES
    # check uncompressed files and compress them
    # to find fastq and fq extensions
    print("1. Checking uncompressed files...\n")

    # this line gives you all .*q files regardless of number of directories in DA*
    uncomp_files = glob.glob("data_raw/DA*/**/*.*q", recursive=True)
    if len(uncomp_files) > 0:
        for line in tqdm(uncomp_files):
            subprocess.run(["pigz", "-p", "%s" % threads, line])
    else:
        print("No uncompressed files found")

    print("2. Collecting compressed files and creating renamed symlinks...\n")

    # to get a list of strains use provided file
    # strains = glob.glob("data_raw/DA*")
    with open(strains_file, 'r') as f:
        strains = [line.rstrip() for line in f.readlines()]

    # 2. CREATING SYMLINKS WITH THE RIGHT NAMES
    for strain in strains:
        # to get a full path to each of GZ files in a given strain:
        read_files = glob.glob("data_raw/%s/**/*.gz" % strain, recursive=True)
        # beginning of each line is a strain name, ending - current filename
        # in the filename - there's a barcode and replicate/read number
        # I need only 2nd barcode from strains with two barcode directories
        # and the only barcode from strains with one barcode directories
        illumina_files = [line for line in read_files if "Illumina" in line]
        nanopore_files = [line for line in read_files if "Nanopore" in line]
        # there can be a case when there's no nanopore or illumina files - we report these cases
        if len(illumina_files) == 0:
            print("No Illumina reads found for %s" % strain)
        if len(nanopore_files) == 0:
            print("No Nanopore reads found for %s" % strain)
            nanopore_clean = list()
        else:
            # filter out nanopore files you don't need
            two_barcodes = [line for line in nanopore_files if line.count("barcode") == 3]
            if len(two_barcodes) == 0:
                nanopore_clean = [line for line in nanopore_files if line.count("barcode") == 2]
                # use these files for symlinks
            else:
                # use two_barcodes for symlinks
                nanopore_clean = two_barcodes

        # join illumina and nanopore lists
        read_files = illumina_files + nanopore_clean
        # iterate and create symlinks
        for line in read_files:
            line_new_name = new_name(line)
            # now, when you have a new name, rename the file (create a soft link)
            # put links into Nanopore or Illumina directories, respectively
            # for symlinks you need absolute paths
            cwd = os.getcwd()
            source = "/".join([cwd, line])
            # create destination - absolute path
            if "Nanopore" in line:
                destination = "/".join([cwd, "data_raw", strain, "Nanopore", line_new_name])
            else:
                # create 'renamed' directory
                os.mkdir("/".join([cwd, "data_raw", strain, "Illumina", "renamed"]))
                destination = "/".join([cwd, "data_raw", strain, "Illumina", "renamed", line_new_name])
            # create link
            try:
                os.symlink(source, destination)
            except FileExistsError:
                print("File exists in destination: %s " % destination)

    # 3. JOINING NANOPORE READS
    # former rule 'join_nanopore': "zcat {input} | pigz -c -p {threads} > {output}"
    print("3. Joining Nanopore reads...")

    for strain in tqdm(strains):
        # strain looks like 'DA62920'
        path = "data_raw/" + strain + "/Nanopore"
        # strain_name = strain.split("/")[-1]
        if os.path.isfile(path + "/" + "%s_all.fastq.gz" % strain):
            print("Joined Nanopore file for strain %s exists" % strain)
        else:
            proc = subprocess.Popen("zcat %s/*.fastq.gz | pigz -c -p %s > %s/%s_all.fastq.gz" %
                                    (path, threads, path, strain),
                                    shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
        # now you can safely remove links, leave only DA[0-9]_all.fastq.gz
        # this is not optimal way but it's easy to embed into this script
        gz_files = glob.glob("%s/*.fastq.gz" % path)
        gz_files = [f for f in gz_files if "_all.fastq.gz" not in f]
        for file in gz_files:
            os.remove(file)
    print("Done!")


if __name__ == '__main__':
    main()
    # TODO: test new functionality
    # cause now fastqc works on all gz files, so it works twice on Nanopore and Illumina data

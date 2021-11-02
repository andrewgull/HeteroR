#!/usr/bin/env python3

"""
this script joins the functionality of three scripts that I use to process raw reads files:
copy_files.py
prepare_files.py
coverage.py
create_config.py
"""

import argparse
import os
import glob
import subprocess
from tqdm import tqdm
import pandas as pd
import yaml


def get_args():
    """
    Get command line arguments
    """

    parser = argparse.ArgumentParser(
        description="This scripts prepares sequencing reads to be fed to the pipeline:\n"
                    "1. Transfers files from Argos to ./data_raw using rsync\n"
                    "2. Renames them like 'STRAIN_[1,2].fq.gz' (Illumina) and 'STRAIN_all.fastq.gz' (Nanopore)\n"
                    "3. Calculates coverage of joined Nanopore reads\n"
                    "4. Creates config file in YAML format\n"
                    "Requirements:\n1. ARGOS must be mounted\n2. seqkit must be installed\n"
                    "--mode flag is not fully functional yet! i.e. only '--mode all' is available.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument("-s", "--strains", type=str, metavar="<file>",
                        help="file with strain names to process, one name per line")
    parser.add_argument("-o", "--output", type=str, metavar="<output file>",
                        help="output file with coverage, tab-separated")
    parser.add_argument("-c", "--config",  type=str, metavar="<config.yaml>",
                        help="config file")
    parser.add_argument("-t", "--threads", type=int, metavar="<number of threads>",
                        help="number of threads to use with pigz", default=4)
    parser.add_argument("-m", "--mode", type=str, metavar="<mode>",
                        help="mode: all, copy, prep, coverage, config (optional)", default="all")
    parser.add_argument("-a", "--argos", type=str, metavar="<path>", help="path to raw data on ARGOS (optional)",
                        default="/home/andrei/Data/Argos/imb_sal_raw/500\ Sepsis\ Eco/Sequencing/Strains")
    parser.add_argument("-l", "--genome_length", type=int, metavar="<genome length>",
                        help="Approx. genome length, bp (optional)", default=5131220)
    parser.add_argument("-v", "--version", action='version', version='%(prog)s 1.0')

    return parser.parse_args()


def copy_files(strain_file, argos_path, destination="data_raw"):
    """A function to copy files from ARGOS except fast5 files"""
    # collect strain names
    with open(strain_file, 'r') as f:
        strains = [line.rstrip() for line in f.readlines()]

    results = list()
    for strain in tqdm(strains):
        # use rsync to copy directories excluding fast5
        # no trailing /
        # no whitespace mirroring in path string
        source = os.path.join(argos_path, strain)
        # if it shows any stdout replace it with subprocess.Popen
        out = subprocess.call("rsync -avrq --exclude='*.fast5' %s %s" % (source, destination), shell=True)
        results.append(out)
        # proc = subprocess.Popen("rsync -avr --exclude='*.fast5' %s %s" % (os.path.join(argos_path, strain),
        # destination), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # out, err = proc.communicate()
    return sum(results)


def new_name(path_string):
    """
    converts path to a filename
    param path: unix path-like string: "data_raw/strain/machine/barcode/filename.fq.gz"
    return: new filename like "strain_14.fq.gz" or "strain_1.fq.gz"
    """
    dirs = path_string.split("/")
    strain_name, old_name = dirs[1], dirs[-1]
    # suffix is universal for both Illumina and Nanopore files
    # it contains read number and extensions
    suffix = old_name.split("_")[-1]
    file_name = "_".join([strain_name, suffix])
    return file_name


def prepare_files(strain_file, threads):
    # 1. COMPRESSING FILES
    # check uncompressed files and compress them
    # to find fastq and fq extensions
    print("2.1. Checking uncompressed files...\n")

    # this line gives you all .*q files regardless of number of directories in DA*
    uncomp_files = glob.glob("data_raw/DA*/**/*.*q", recursive=True)
    if len(uncomp_files) > 0:
        print("Found %i uncompressed fastq files. From\n%s\nto\n%s"
              % (len(uncomp_files), uncomp_files[0], uncomp_files[-1]))
        for line in tqdm(uncomp_files):
            # force overwrite
            subprocess.run(["pigz", "-f", "-p", "%i" % threads, line])
    else:
        print("No uncompressed files found")

    print("2.2. Collecting compressed files and creating renamed symlinks...\n")

    # to get a list of strains use provided file
    # strains = glob.glob("data_raw/DA*")
    with open(strain_file, 'r') as f:
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
                if len(nanopore_clean) == 0:
                    # if dir is called 'Fastq_pass' and has "pass_barcode" in it
                    nanopore_clean = [line for line in nanopore_files if "Fastq_pass" in line]
                # use these files for symlinks
            else:
                # use two_barcodes for symlinks
                nanopore_clean = two_barcodes

        # join illumina and nanopore lists
        read_files = illumina_files + nanopore_clean
        # iterate and create symlinks
        for line in tqdm(read_files):
            line_new_name = new_name(line)
            # now, when you have a new name, rename the file (create a soft link)
            # put links into Nanopore or Illumina directories, respectively
            # for symlinks you need absolute paths
            cwd = os.getcwd()
            source = "/".join([cwd, line])
            # create destination - absolute path
            if "Nanopore" in line:
                destination = os.path.join(cwd, "data_raw", strain, "Nanopore", line_new_name)
            else:
                # create 'renamed' directory
                renamed = os.path.join(cwd, "data_raw", strain, "Illumina", "renamed")
                if not os.path.exists(renamed):
                    os.mkdir(renamed)
                destination = os.path.join(renamed, line_new_name)
            # create link
            try:
                os.symlink(source, destination)
            except FileExistsError:
                print("File exists in destination: %s " % destination)

    # 3. JOINING NANOPORE READS
    # former rule 'join_nanopore': "zcat {input} | pigz -c -p {threads} > {output}"
    print("2.3. Joining Nanopore reads...")

    for strain in tqdm(strains):
        # strain looks like 'DA62920'
        path = "data_raw/" + strain + "/Nanopore"
        # strain_name = strain.split("/")[-1]
        if os.path.isfile(path + "/" + "%s_all.fastq.gz" % strain):
            print("Joined Nanopore file for strain %s exists" % strain)
        else:
            proc = subprocess.Popen("zcat %s/*.fastq.gz | pigz -c -p %i > %s/%s_all.fastq.gz" %
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

    # 4. COMPRESS FAST5 FILES
    print("2.4. Compressing FAST5 files...")
    # collect all fast5 files
    fast5 = glob.glob("data_raw/DA*/**/*.fast5", recursive=True)
    if len(fast5) > 0:
        print("Found %i fastq5 files" % len(fast5))
        for line in tqdm(fast5):
            subprocess.run(["pigz", "-p", "%i" % threads, line])
    else:
        print("No uncompressed fast5 files found!")

    return None


def coverage(strain_file, genome_length):
    """Calculate and summarize coverage"""
    with open(strain_file, 'r') as f:
        strain_names = [line.rstrip() for line in f.readlines()]

    # collect Nanopore all
    nanopore_stats = list()
    for strain in tqdm(strain_names):
        try:
            file = glob.glob("data_raw/%s/Nanopore/%s_all.fastq.gz" % (strain, strain))[0]
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
    nanopore_df["coverage"] = nanopore_df["sum_len"] / genome_length
    return nanopore_df


def create_config(strain_file):
    """create config dictionary"""
    # read file
    with open(strain_file, 'r') as f:
        strains = [line.rstrip() for line in f.readlines()]

    # create dict
    config = dict({"strains": {}})

    # fill in the dict
    for strain in strains:
        config["strains"][strain] = strain

    return config


if __name__ == '__main__':
    args = get_args()

    # 1. Copy files
    print("1. Transfer files from ARGOS...")
    copy_results = copy_files(strain_file=args.strains, argos_path=args.argos)
    if copy_results > 0:
        print("WARNING! Something went wrong during raw files transfer: at least one process finished with exit code 1")

    # 2. Prepare files
    print("\n2. Preparing files...")
    prepare_files(strain_file=args.strains, threads=args.threads)

    # 3. Calculate coverage
    print("\n3. Calculating coverage...")
    coverage_stats = coverage(strain_file=args.strains, genome_length=args.genome_length)
    min_cov, max_cov, avg_cov = coverage_stats["coverage"].min(), coverage_stats["coverage"].max(), coverage_stats["coverage"].mean()
    # write it to a file
    coverage_stats.to_csv(path_or_buf=args.output, sep="\t", index=False)
    # print some stats
    print("Quick stats:\nmin = %f\navg = %f\nmax = %f\nCoverage ~25x or less is sparse, good for Unicycler.\n" % (min_cov, avg_cov, max_cov))

    # 4. Create a config file
    print("\n4. Creating a config file...")
    config_dict = create_config(strain_file=args.strains)
    # write as yaml
    with open(args.config, 'w') as outfile:
        yaml.dump(config_dict, outfile, default_flow_style=False)
    print("All done!\nNow you can run the pipeline.")

#!/usr/bin/env python3
"""
decide to trim or not to trim using a logit regression equation
trim Nanopore reads if probability is higher than a threshold
"""

import math
import subprocess
import shutil
import os


def do_not_trim():
    return 0


def do_trim():
    # run filtlong
    return 0


# inputs
input_file = snakemake.input[0]  # "resources/data_raw/{strain}/Nanopore/{strain}_all.fastq.gz"
output_file = snakemake.output[0]  # "results/data_filtered/{strain}/Nanopore/{strain}_all.fastq.gz"
# coverage and probability params
genome_len = snakemake.params[0]  # 5131220
intercept = snakemake.params[1]  # 2.2401998
coefficient = snakemake.params[2]  # -0.2889426
cut_off = snakemake.params[3]  # 0.70
# filtlong params
min_len = snakemake.params[4]
len_weight = snakemake.params[5]
percent = snakemake.parmas[6]
bases = snakemake.params[7]
threads = snakemake.params[8]

# calculate the coverage
proc = subprocess.Popen("seqkit stats %s -T | cut -f5 | tail -n 1" % input_file, shell=True, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
out_cov, err_cov = proc.communicate()
sum_len = int(out_cov.decode("utf-8"))
coverage = sum_len / genome_len
# here comes the smart formula from the logistic regression model
odds = math.exp(intercept + coefficient * coverage)
# probability of being incomplete after trimming
probability = odds / (1 + odds)

# create output directories if they do not exits
strain = output_file.split("/")[2]
output_directory = "/".join(output_file.split("/")[:-1])  # 'results/data_filtered/{strain}/Nanopore'



if probability >= cut_off:
    # very likely the assembly will be incomplete
    do_not_trim()
    # just copy the file to results
    # create a directory first
    if not os.path.exists(output_directory):
        os.mkdir("results/data_filtered/%s/Nanopore" % strain)
    else:

    shutil.copy(input_file, output_file)

else:
    do_trim()
    # here piping may not work well, split in two commands if so
    subprocess.Popen("filtlong --min_length %i --length_weight %i --keep_percent %i --target_bases %i %s "
                     "| pigz -c -p %i > %s" % (min_len, len_weight, percent, bases, input_file, threads, output_file),
                     shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # both go to the snakemake's log
    out_filter, err_filter = proc.communicate()

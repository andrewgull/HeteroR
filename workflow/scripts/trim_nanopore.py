#!/usr/bin/env python3
"""
decide to trim or not to trim using a logit regression equation
trim Nanopore reads if probability is higher than a threshold
"""

import math
import subprocess
import shutil
import os


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
percent = snakemake.params[6]
bases = snakemake.params[7]
threads = snakemake.params[8]

log_messages = list()

# calculate the coverage
proc = subprocess.Popen("seqkit stats %s -T | cut -f5 | tail -n 1" % input_file, shell=True, stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
out_cov, err_cov = proc.communicate()
log_messages.append("retrieved total length:")
log_messages.append(out_cov)
log_messages.append(err_cov)

sum_len = int(out_cov.decode("utf-8"))
coverage = sum_len / genome_len
log_messages.append("calculated coverage:" + str(coverage))
# here comes the smart formula from the logistic regression model
odds = math.exp(intercept + coefficient * coverage)
# probability of being incomplete after trimming
probability = odds / (1 + odds)
log_messages.append("calculated probability of incomplete assembly:" + str(probability))

# create output directories if they do not exits
strain = output_file.split("/")[2]
output_directory = "/".join(output_file.split("/")[:-1])  # 'results/data_filtered/{strain}/Nanopore'
# 'data_filtered' may not exist, then create 'data_filtered/strain' and 'data_filtered/strainNanopore'
# in a sequential way because this is how os.mkdir works
if not os.path.exists("results/data_filtered"):
    os.mkdir("results/data_filtered")
    os.mkdir("results/data_filtered/%s" % strain)
    os.mkdir("results/data_filtered/%s/Nanopore" % strain)
# now the full path exists

if probability >= cut_off:
    # very likely the assembly will be incomplete
    # just copy the file to results
    shutil.copy(input_file, output_file)
else:
    # output file directory should have been created
    # here piping may not work well, split in two commands if so
    proc = subprocess.Popen("filtlong --min_length %i --length_weight %i --keep_percent %i --target_bases %i %s "
                     "| pigz -c -p %i > %s" % (min_len, len_weight, percent, bases, input_file, threads, output_file),
                     shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # both go to the snakemake's log
    out_filter, err_filter = proc.communicate()
    log_messages.append(out_filter)
    log_messages.append(err_filter)

# write log files
with open(snakemake.log[0], "w") as log_out:
    for message in log_messages:
        log_out.write(str(message) + "\n")

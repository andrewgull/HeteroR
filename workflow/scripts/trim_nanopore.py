#!/usr/bin/env python3
"""
decide to trim or not to trim using a logit regression equation
trim Nanopore reads if probability is higher than a threshold
"""

import math
import subprocess


def do_not_trim():
    return 0


def do_trim():
    # run filtlong
    return 0


# inputs and parameters
input_file = snakemake.input[0]
output_file = snakemake.output[0]
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

proc = subprocess.Popen("seqkit stats %s -T | cut -f5 | tail -n 1" % input_file, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = proc.communicate()
sum_len = int(out.decode("utf-8"))
coverage = sum_len / genome_len
# here comes the smart formula from the logistic regression model
odds = math.exp(intercept + coefficient * coverage)
# probability of being incomplete after trimming
probability = odds / (1 + odds)

if probability >= cut_off:
    # very likely the assembly will be incomplete
    do_not_trim()
else:
    do_trim()
    subprocess.Popen("filtlong --min_length {params.min_len} --length_weight {params.len_weight} --keep_percent {params.perc} --target_bases {params.bases} {input} 2> {log} | pigz -c -p {threads} > {output}")


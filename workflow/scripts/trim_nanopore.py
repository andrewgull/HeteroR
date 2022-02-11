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
genome_len = int(snakemake.params[0])  # 5131220
intercept = float(snakemake.params[1])  # 2.2401998
coefficient = float(snakemake.params[2])  # -0.2889426
cut_off = float(snakemake.params[3])  # 0.70

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


#!/usr/bin/env python
import sys
import subprocess

file_in = sys.argv[1]
file_out = sys.argv[2]
threads = sys.argv[3]
# the following line is required because pigz doesn't need output name,
# but the rule join_nanopore needs it, and it must be *_fastq.gz
true_output = file_out[:-3]

# the following line is required because the rule join_nanopore requires input,
# but number of input files is always unknown
true_input = file_in[-9:]  # this should be ".fastq.gz"

# "zcat {input} > {output} && pigz -p {threads} {output}"
subprocess.call(["zcat", "*%s" % true_input, ">", "%s" % true_output], shell=True)
subprocess.call(["pigz", "-p", "%s" % threads, "%s" % true_output], shell=True)

#!/usr/bin/env python
"""
a script for running quast.py, because direct running of QUAST tool makes it impossible
to connect its rule with unicycler rule
directive 'run' is also not viable cause than it's impossible to use conda envs
"""

import subprocess
import snakemake
import os


def main(indir, outdir, cpus):
    """
    function to run quast.py
    """
    # args = get_args()
    # cpus = args.threads
    # output = args.file_out
    # infile = args.file_in
    # INPUT looks like "assemblies/{strain}" but you need "assemblies/{strain}/assembly.fasta"
    infile = os.path.join(indir, "assembly.fasta")
    subprocess.call(["quast.py", "-t", cpus, "-o", outdir, infile], shell=True)


main(indir=snakemake.input[0], outdir=snakemake.output[0], cpus=snakemake.threads)

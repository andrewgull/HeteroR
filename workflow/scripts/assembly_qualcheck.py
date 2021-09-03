#!/usr/bin/env python
"""
a script for running BUSCO, because direct running of BUSCO tool makes it impossible
to connect its rule with unicycler rule
directive 'run' is also not viable cause than it's impossible to use conda envs
"""

import subprocess
# import snakemake
import os


def main(indir, outdir, cpus):
    """
    function to run busco
    """

    # run BUSCO
    # busco -m genome -i {input} -o busco_results --out_path {output} -l gammaproteobacteria_odb10 --cpu {threads}"
    # INPUT looks like "assemblies/{strain}" but you need "assemblies/{strain}/assembly.fasta"
    infile = os.path.join(indir, "assembly.fasta")
    subprocess.call(["busco", "-m", "genome", "-i", infile, "-o", "busco_results", "--out_path", outdir,
                     "-l", "gammaproteobacteria_odb10", "--cpu", cpus])
    # run QUAST
    output_dir = os.path.join(outdir, "quast_results")
    subprocess.call(["quast.py", "-t", cpus, "-o", output_dir, infile])


main(indir=snakemake.input[0], outdir=snakemake.output[0], cpus=str(snakemake.threads))

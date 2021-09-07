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
    # subprocess.call(["busco", "-m", "genome", "-i", infile, "-o", "busco_results", "--out_path", outdir,
    #                 "-l", "gammaproteobacteria_odb10", "--cpu", cpus])
    proc_busco = subprocess.Popen("busco -m genome -i %s -o busco_results --out_path %s -l gammaproteobacteria_odb10 "
                                  "--cpu %s"
                                  % (infile, outdir, cpus), shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc_busco.communicate()
    strain = indir.split("/")[1]
    # write logs
    with open("logs/%s_busco_stderr.log" % strain, "w") as errout:
        errout.write(err.decode("utf-8"))

    with open("logs/%s_busco_stdout.log" % strain, "w") as stout:
        stout.write(out.decode("utf-8"))

    # run QUAST
    output_dir = os.path.join(outdir, "quast_results")
    subprocess.call(["quast.py", "-t", cpus, "-o", output_dir, infile])


main(indir=snakemake.input[0], outdir=snakemake.output[0], cpus=str(snakemake.threads))

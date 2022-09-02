#!/usr/bin/env python
"""
a script for running BUSCO, because direct running of BUSCO tool makes it impossible
to connect its rule with unicycler rule
directive 'run' is also not viable cause than it's impossible to use conda envs
"""

import subprocess
import os


def write_logs(output, tool, strain_name):
    """
    :param output: stdout or stderr output of communicate(), bytes object
    :param tool: busco or quast
    :param strain_name: name of strain
    """
    if len(output) > 0:
        with open("results/logs/%s_%s.log" % (strain_name, tool), "w") as f:
            f.write(output.decode("utf-8"))
    else:
        pass


def main(in_dir, out_dir, cpus, db_path, tax_ds):
    """
    function to run busco and quast with logging
    :param in_dir: input directory
    :param out_dir: output directory
    :param cpus: number of threads
    :param db_path: path to BUSCO database
    :param tax_ds: taxonomic data set for BUSCO, to print the full list run: busco --list-datasets
    """

    # INPUT looks like "results/assemblies/{strain}" but you need "results/assemblies/{strain}/assembly.fasta"
    infile = os.path.join(in_dir, "assembly.fasta")
    # for QUAST we need a special output_dir
    output_dir = os.path.join(out_dir, "quast_results")
    # for logs we need a strain name
    strain = in_dir.split("/")[1]

    # run BUSCO
    # subprocess.call(["busco", "-m", "genome", "-i", infile, "-o", "busco_results", "--out_path", outdir,
    #                 "-l", "gammaproteobacteria_odb10", "--cpu", cpus])
    proc_busco = subprocess.Popen("busco -m genome -i %s -o busco_results --out_path %s -l %s "
                                  "--cpu %s --download_path %s"
                                  % (infile, out_dir, tax_ds, cpus, db_path), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    busco_stdout, busco_stderr = proc_busco.communicate()

    # run QUAST
    # subprocess.call(["quast.py", "-t", cpus, "-o", output_dir, infile])
    proc_quast = subprocess.Popen("quast.py -t %s -o %s %s" % (cpus, output_dir, infile),
                                  shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    quast_stdout, quast_stderr = proc_quast.communicate()

    # write logs
    output_dict = dict({"busco": (busco_stderr, busco_stdout), "quast": (quast_stderr, quast_stdout)})
    for software in output_dict.keys():
        for item in output_dict[software]:
            write_logs(item, software, strain)


main(in_dir=snakemake.input[0], out_dir=snakemake.output[0], cpus=str(snakemake.threads), db_path=snakemake.input[1], tax_ds=snakemake.params[0])

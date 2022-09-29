# script for adaptive hybrid assembling
# runs unicycler and then flye-medaka-polypolish if chromosomal assembly is not complete
import snakemake
import subprocess
import os

# OPEN LOG
log = open(snakemake.log[0], "w")
# log.write("log text")

# DEFINE INPUTS/OUTPUTS
short_reads_1 = snakemake.input[0]
short_reads_2 = snakemake.input[1]
long_reads = snakemake.input[2]
assembly_dir = snakemake.output[0]
draft_dir = snakemake.output[1]
polish_dir = snakemake.output[2]
unicycler_log_path = os.path.join(assembly_dir, "unicycler.log")

# DEFINE PARAMS
threads = snakemake.threads[0]
basecaller = snakemake.params[0]
genome_size = snakemake.params[1]
coverage = snakemake.params[2]

# RUN UNICYCLER
subprocess.run("unicycler -1 %s -2 %s -l %s -t %i -o %s " % (short_reads_1, short_reads_2, long_reads, threads, assembly_dir), shell=True)

# CHECK ASSEMBLY COMPLETENESS
completeness = subprocess.run("sed -n '/^Component/,/^Polishing/{p;/^Polishing/q}' %s | head -n -3 | tr -s ' ' | "
                              "cut -d ' ' -f 8" % unicycler_log_path, shell=True, capture_output=True)
completeness_stdout = completeness.stdout.decode("utf-8").splitlines()

if completeness_stdout[2] == "incomplete":
    # RUN FLYE
    subprocess.run("flye --nano-raw %s --threads %i --out-dir %s -g %s --asm-coverage %i" %
                   (long_reads, threads, assembly_dir, genome_size, coverage), shell=True)
    # RUN MEDAKA
    subprocess.run("medaka_consensus -i %s -d %s/assembly.fasta -o %s -t %i -m %s" %
                   (long_reads, assembly_dir, draft_dir, threads, basecaller), shell=True)
    # MAP ILLUMINA READS ONTO MEDAKA CONSENSUS
    subprocess.run("bwa index %s/consensus.fasta" % draft_dir, shell=True)
    subprocess.run("bwa mem -t %i -a %s/consensus.fasta %s > %s/alignments_1.sam" %
                   (threads, draft_dir, short_reads_1, polish_dir), shell=True)
    subprocess.run("bwa mem -t %i -a %s/consensus.fasta %s > %s/alignments_2.sam" %
                   (threads, draft_dir, short_reads_1, polish_dir), shell=True)
    # POLYPOLISH
    subprocess.run("polypolish_insert_filter.py --in1 %s/alignments_1.sam --in2 %s/alignments_2.sam "
                   "--out1 %s/filtered_1.sam --out2 %s/filtered_2.sam" %
                   (polish_dir, polish_dir, polish_dir, polish_dir), shell=True)
    subprocess.run("polypolish %s/consensus.fasta %s/filtered_1.sam %s/filtered_2.sam > "
                   "%s/assembly_polished.fasta" % (draft_dir, polish_dir, polish_dir, polish_dir), shell=True)
    # RENAME UNICYCLER ASSEMBLY
    os.rename("%s/assembly.fasta" % assembly_dir, "%s/assembly_unicycler.fasta" % assembly_dir)
    # LINK ASSEMBLY TO RESULTS/ASSEMBLIES
    cwd = os.getcwd()
    # your polished assembly
    source = os.path.join(cwd, "%s/assembly_polished.fasta" % polish_dir)
    # same file as for unicycler
    destination = os.path.join(cwd, "%s/assembly.fasta" % assembly_dir)
    os.symlink(source, destination)
    print("flye-medaka-polypolish's been chosen")
else:
    print("unicycler's been chosen")

# CLOSE LOG
log.close()

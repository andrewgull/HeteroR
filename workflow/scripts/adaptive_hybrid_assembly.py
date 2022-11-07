# script for adaptive hybrid assembling
# runs unicycler and then flye-medaka-polypolish if chromosomal assembly is not complete
import sys
import subprocess
import os

# DEFINE INPUTS/OUTPUTS
short_reads_1 = snakemake.input[0]
short_reads_2 = snakemake.input[1]
long_reads = snakemake.input[2]
assembly_dir = snakemake.output[0]
draft_dir = snakemake.output[1]
polish_dir = snakemake.output[2]
unicycler_log_path = os.path.join(assembly_dir, "unicycler.log")

# DEFINE PARAMS
threads = snakemake.threads
basecaller = snakemake.params[0]
genome_size = snakemake.params[1]
coverage = snakemake.params[2]

# LIST TO KEEP STDOUT/STDERR
outs = list()

# OPEN LOG
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    # RUN UNICYCLER
    unicycler_out = subprocess.run("unicycler -1 %s -2 %s -l %s -t %i -o %s "
                                   % (short_reads_1, short_reads_2, long_reads, threads, assembly_dir),
                                   shell=True, capture_output=True, text=True)
    # ADD CAPTURED OUT TO OUTS
    outs.append(unicycler_out)

    # CHECK ASSEMBLY COMPLETENESS
    unicycler_log_path = "results/assemblies/DA63356/unicycler.log"
    completeness = subprocess.run("sed -n '/^Component/,/^Polishing/{p;/^Polishing/q}' %s | head -n -3 | tr -s ' ' | "
                                  "cut -d ' ' -f 8" % unicycler_log_path, shell=True, capture_output=True, text=True)
    completeness_stdout = completeness.stdout.splitlines()
    # if this table with completeness in unicycler.log contains only 1 assembled component
    # there is no 'total' line in this table
    # therefore the length of completeness_stdout equals 2
    # and index of the chromosome completeness status should be 1, not 2
    # on the next line you should check this
    if len(completeness_stdout) == 2:
        chrom_status_index = 1
    else:
        chrom_status_index = 2
    # now use this index to determine chromosome completeness status
    if completeness_stdout[chrom_status_index] == "incomplete":
        print("Unicycler assembly is not complete.\nFlye-Medaka-Polypolish have been chosen to produce new assembly.")
        os.rename("%s/assembly.fasta" % assembly_dir, "%s/assembly_unicycler.fasta" % assembly_dir)
        # RUN FLYE
        flye_out = subprocess.run("flye --nano-raw %s --threads %i --out-dir %s -g %s --asm-coverage %i" %
                                  (long_reads, threads, assembly_dir, genome_size, coverage),
                                  shell=True, capture_output=True, text=True)
        outs.append(flye_out)

        # RUN MEDAKA
        medaka_out = subprocess.run("medaka_consensus -i %s -d %s/assembly.fasta -o %s -t %i -m %s" %
                                    (long_reads, assembly_dir, draft_dir, threads, basecaller),
                                    shell=True, capture_output=True, text=True)
        outs.append(medaka_out)

        # MAP ILLUMINA READS ONTO MEDAKA CONSENSUS
        os.mkdir(polish_dir)  # you have to create the output dir here
        # 1. INDEX ASSEMBLY
        index_out = subprocess.run("bwa index %s/consensus.fasta" % draft_dir,
                                   shell=True, capture_output=True, text=True)
        outs.append(index_out)
        # 2. MAP READ SET 1
        mem1_out = subprocess.run("bwa mem -t %i -a %s/consensus.fasta %s > %s/alignments_1.sam" %
                                  (threads, draft_dir, short_reads_1, polish_dir),
                                  shell=True, capture_output=True, text=True)
        outs.append(mem1_out)
        # 3. MAP READ SET 2
        mem2_out = subprocess.run("bwa mem -t %i -a %s/consensus.fasta %s > %s/alignments_2.sam" %
                                  (threads, draft_dir, short_reads_2, polish_dir),
                                  shell=True, capture_output=True, text=True)
        outs.append(mem2_out)

        # POLYPOLISH
        # 1. INSERT FILTER
        plp_insert_out = subprocess.run("polypolish_insert_filter.py --in1 %s/alignments_1.sam "
                                        "--in2 %s/alignments_2.sam --out1 %s/filtered_1.sam --out2 %s/filtered_2.sam" %
                                        (polish_dir, polish_dir, polish_dir, polish_dir),
                                        shell=True, capture_output=True, text=True)
        outs.append(plp_insert_out)
        # 2. POLISH
        plp_polish_out = subprocess.run("polypolish %s/consensus.fasta %s/filtered_1.sam %s/filtered_2.sam > "
                                        "%s/assembly_flye_polished.fasta" %
                                        (draft_dir, polish_dir, polish_dir, polish_dir),
                                        shell=True, capture_output=True, text=True)
        outs.append(plp_polish_out)

        # REMOVE SAM FILES
        sam_files = ["alignments_1.sam", "alignments_2.sam", "filtered_1.sam", "filtered_2.sam"]
        for sam in sam_files:
            os.remove(os.path.join(polish_dir, sam))
        
        print("SAM files have been removed")

        # RENAME FLYE ASSEMBLY
        os.rename("%s/assembly.fasta" % assembly_dir, "%s/assembly_flye_raw.fasta" % assembly_dir)
        # LINK ASSEMBLY TO RESULTS/ASSEMBLIES
        cwd = os.getcwd()
        # your polished assembly
        source = os.path.join(cwd, "%s/assembly_flye_polished.fasta" % polish_dir)
        # same file as for unicycler
        destination = os.path.join(cwd, "%s/assembly.fasta" % assembly_dir)
        os.symlink(source, destination)
    else:
        # create dirs: polished, drafts
        os.mkdir(draft_dir)
        os.mkdir(polish_dir)
        print("Unicycler assembly is complete. Empty draft and polish dirs were created.")

    # PRINT STDOUT/STDERR TO LOG
    for out in outs:
        # every out is subrocess.run output not string! 
        print(out.stdout, out.stderr)

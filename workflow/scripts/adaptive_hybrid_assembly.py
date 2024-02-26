# adaptive hybrid assembling
# this script runs unicycler if covergae is shallow (< 20x) or flye-medaka-polypolish
# the latter can also be run if unicycler crashes
# all constants except for UNICYLCER_STATUS are defined in the Snakemake file

import sys
import subprocess
import os

# DEFINE INPUTS/OUTPUTS
SHORT_READS_1 = snakemake.input[0]
SHORT_READS_2 = snakemake.input[1]
LONG_READS = snakemake.input[2]
ASSEMBLY_DIR = snakemake.output[0]
DRAFT_DIR = snakemake.output[1]
POLISH_DIR = snakemake.output[2]

# DEFINE PARAMS
THREADS = snakemake.threads
BASECALLER = snakemake.params[0]
GENOME_SIZE = snakemake.params[1]
COVERAGE = snakemake.params[2]   # coverage for Flye
GENOME_LENGTH = snakemake.params[3]  # expected genome length
COV_THRESHOLD = snakemake.params[4]  # threshold value to choose assembler, < 20x is good for Uni

# LIST TO KEEP STDOUT/STDERR
outs = []

# OPEN LOG
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    
    # CALCULATE COVERAGE in order to decide which assembler to use
    seqkit_out = subprocess.run(f"seqkit stats {LONG_READS} -T", shell=True, capture_output=True, text=True, check=True)
    gen_coverage = int(seqkit_out.stdout.split("\n")[1].split("\t")[4])/GENOME_LENGTH

    # SET UNICYCLER STATUS, it is needed to decide whether to ditch Unicycler if it crashes
    UNICYCLER_STATUS = "ok"

    # RUN UNICYCLER IF COVERAGE IS SHALLOW
    if gen_coverage <= COV_THRESHOLD:
        # RUN UNICYCLER
        unicycler_out = subprocess.run(f"unicycler -1 {SHORT_READS_1} -2 {SHORT_READS_2} -l {LONG_READS} -t {THREADS} -o {ASSEMBLY_DIR}",
                                   shell=True, capture_output=True, text=True, check=True)
        # ADD CAPTURED OUTPUT TO OUTS
        outs.append(unicycler_out)

        # CHECK FOR ERRORS
        if len(unicycler_out.stderr) > 0:
            UNICYCLER_STATUS = "error"
            print("The genome coverage is sparse, but Unicycler has crashed.")
        else:
            # all good, no need to run FMP
            # create dirs: polished, drafts
            os.mkdir(DRAFT_DIR)
            os.mkdir(POLISH_DIR)
            print("The genome coverage is sparse, Unicycler has been chosen.\nEmpty draft and polish dirs have been created.")
    
    # RUN FMP IF UNICYCLER BROKE OR COVERAGE IS SPARSE
    elif gen_coverage > COV_THRESHOLD or UNICYCLER_STATUS == "error":
        # MK ASSEMBLY DIR but it can exist after previous runs of the pipeline
        if not os.path.exists(ASSEMBLY_DIR):
            os.mkdir(ASSEMBLY_DIR)
        
        # RUN FLYE
        print("Flye-Medaka-Polypolish will be chosen to assemble the genome")
        flye_out = subprocess.run(f"flye --nano-raw {LONG_READS} --threads {THREADS} --out-dir {ASSEMBLY_DIR} -g {GENOME_SIZE} --asm-coverage {gen_coverage}",
                                  shell=True, capture_output=True, text=True, check=True)
        outs.append(flye_out)

        # RUN MEDAKA
        medaka_out = subprocess.run(f"medaka_consensus -i {LONG_READS} -d {ASSEMBLY_DIR}/assembly.fasta -o {DRAFT_DIR} -t {THREADS} -m {BASECALLER}",
                                    shell=True, capture_output=True, text=True, check=True)
        outs.append(medaka_out)

        # MAP ILLUMINA READS ONTO MEDAKA CONSENSUS
        os.mkdir(POLISH_DIR)  # you have to create the output dir here
        # 1. INDEX ASSEMBLY
        index_out = subprocess.run(f"bwa index {DRAFT_DIR}/consensus.fasta",
                                   shell=True, capture_output=True, text=True, check=True)
        outs.append(index_out)
        # 2. MAP READ SET 1
        mem1_out = subprocess.run(f"bwa mem -t {THREADS} -a {DRAFT_DIR}/consensus.fasta {SHORT_READS_1} > {POLISH_DIR}/alignments_1.sam",
                                  shell=True, capture_output=True, text=True, check=True)
        outs.append(mem1_out)
        # 3. MAP READ SET 2
        mem2_out = subprocess.run(f"bwa mem -t {THREADS} -a {DRAFT_DIR}/consensus.fasta {SHORT_READS_2} > {POLISH_DIR}/alignments_2.sam",
                                  shell=True, capture_output=True, text=True, check=True)
        outs.append(mem2_out)

        # POLYPOLISH
        # 1. INSERT FILTER
        plp_insert_out = subprocess.run(f"polypolish_insert_filter.py --in1 {POLISH_DIR}/alignments_1.sam "
                                        f"--in2 {POLISH_DIR}/alignments_2.sam --out1 {POLISH_DIR}/filtered_1.sam --out2 {POLISH_DIR}/filtered_2.sam",
                                        shell=True, capture_output=True, text=True, check=True)
        outs.append(plp_insert_out)
        # 2. POLISH
        plp_polish_out = subprocess.run(f"polypolish {DRAFT_DIR}/consensus.fasta {POLISH_DIR}/filtered_1.sam {POLISH_DIR}/filtered_2.sam > "
                                        f"{POLISH_DIR}/assembly_flye_polished.fasta",
                                        shell=True, capture_output=True, text=True, check=True)
        outs.append(plp_polish_out)

        # REMOVE SAM FILES
        sam_files = ["alignments_1.sam", "alignments_2.sam", "filtered_1.sam", "filtered_2.sam"]
        for sam in sam_files:
            os.remove(os.path.join(POLISH_DIR, sam))
        print("SAM files have been removed")

        # RENAME FLYE ASSEMBLY
        os.rename(f"{ASSEMBLY_DIR}/assembly.fasta", f"{ASSEMBLY_DIR}/assembly_flye_raw.fasta")
        # LINK ASSEMBLY TO RESULTS/ASSEMBLIES
        cwd = os.getcwd()
        # your polished assembly
        source = os.path.join(cwd, f"{POLISH_DIR}/assembly_flye_polished.fasta")
        # same file as for unicycler
        destination = os.path.join(cwd, f"{ASSEMBLY_DIR}/assembly.fasta")
        os.symlink(source, destination)

    # PRINT STDOUT/STDERR TO LOG
    for out in outs:
        # every out is subrocess.run output not string!
        print(out.stdout, out.stderr)

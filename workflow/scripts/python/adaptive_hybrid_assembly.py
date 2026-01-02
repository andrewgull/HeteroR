# adaptive hybrid assembling
# this script runs unicycler if covergae is shallow (< 20x) or flye-medaka-polypolish
# the latter can also be run if unicycler crashes
# all constants except for UNICYLCER_STATUS are defined in the Snakemake file

import sys
import subprocess
import os

def calculate_coverage(long_reads, genome_length):
    seqkit_out = subprocess.run(f"seqkit stats {long_reads} -T", shell=True, capture_output=True, text=True, check=True)
    gen_coverage = int(int(seqkit_out.stdout.split("\n")[1].split("\t")[4])/genome_length)
    return gen_coverage, seqkit_out

def run_unicycler(sr1, sr2, lr, threads, assembly_dir, draft_dir, polish_dir):
    # RUN UNICYCLER
    unicycler_out = subprocess.run(f"unicycler -1 {sr1} -2 {sr2} -l {lr} -t {threads} -o {assembly_dir}",
                                   shell=True, capture_output=True, text=True, check=True)
    
    # CHECK FOR ERRORS
    if len(unicycler_out.stderr) > 0:
        return unicycler_out, "error"
    else:
        # all good, no need to run FMP
        # create dirs: polished, drafts
        if not os.path.exists(draft_dir):
            os.mkdir(draft_dir)
        if not os.path.exists(polish_dir):
            os.mkdir(polish_dir)
        return unicycler_out, "ok"

def run_fmp(sr1, sr2, lr, threads, assembly_dir, draft_dir, polish_dir, basecaller, genome_size, gen_coverage):
    outs = []
    # MK ASSEMBLY DIR but it can exist after previous runs of the pipeline
    if not os.path.exists(assembly_dir):
        os.mkdir(assembly_dir)
    
    # RUN FLYE
    print("Flye-Medaka-Polypolish will be chosen to assemble the genome")
    flye_out = subprocess.run(f"flye --nano-raw {lr} --threads {threads} --out-dir {assembly_dir} -g {genome_size} --asm-coverage {gen_coverage}",
                              shell=True, capture_output=True, text=True, check=True)
    outs.append(flye_out)

    # RUN MEDAKA
    medaka_out = subprocess.run(f"medaka_consensus -i {lr} -d {assembly_dir}/assembly.fasta -o {draft_dir} -t {threads} -m {basecaller}",
                                shell=True, capture_output=True, text=True, check=True)
    outs.append(medaka_out)

    # MAP ILLUMINA READS ONTO MEDAKA CONSENSUS
    if not os.path.exists(polish_dir):
        os.mkdir(polish_dir)
    
    # 1. INDEX ASSEMBLY
    index_out = subprocess.run(f"bwa index {draft_dir}/consensus.fasta",
                               shell=True, capture_output=True, text=True, check=True)
    outs.append(index_out)
    # 2. MAP READ SET 1
    mem1_out = subprocess.run(f"bwa mem -t {threads} -a {draft_dir}/consensus.fasta {sr1} > {polish_dir}/alignments_1.sam",
                              shell=True, capture_output=True, text=True, check=True)
    outs.append(mem1_out)
    # 3. MAP READ SET 2
    mem2_out = subprocess.run(f"bwa mem -t {threads} -a {draft_dir}/consensus.fasta {sr2} > {polish_dir}/alignments_2.sam",
                              shell=True, capture_output=True, text=True, check=True)
    outs.append(mem2_out)

    # POLYPOLISH
    # 1. INSERT FILTER
    plp_insert_out = subprocess.run(f"polypolish_insert_filter.py --in1 {polish_dir}/alignments_1.sam "
                                    f"--in2 {polish_dir}/alignments_2.sam --out1 {polish_dir}/filtered_1.sam --out2 {polish_dir}/filtered_2.sam",
                                    shell=True, capture_output=True, text=True, check=True)
    outs.append(plp_insert_out)
    # 2. POLISH
    plp_polish_out = subprocess.run(f"polypolish {draft_dir}/consensus.fasta {polish_dir}/filtered_1.sam {polish_dir}/filtered_2.sam > "
                                    f"{polish_dir}/assembly_flye_polished.fasta",
                                    shell=True, capture_output=True, text=True, check=True)
    outs.append(plp_polish_out)

    # REMOVE SAM FILES
    sam_files = ["alignments_1.sam", "alignments_2.sam", "filtered_1.sam", "filtered_2.sam"]
    for sam in sam_files:
        p = os.path.join(polish_dir, sam)
        if os.path.exists(p):
            os.remove(p)
    print("SAM files have been removed")

    # RENAME FLYE ASSEMBLY
    if os.path.exists(f"{assembly_dir}/assembly.fasta"):
        os.rename(f"{assembly_dir}/assembly.fasta", f"{assembly_dir}/assembly_flye_raw.fasta")
    
    # LINK ASSEMBLY TO RESULTS/ASSEMBLIES
    cwd = os.getcwd()
    source = os.path.join(cwd, f"{polish_dir}/assembly_flye_polished.fasta")
    destination = os.path.join(cwd, f"{assembly_dir}/assembly.fasta")
    if os.path.exists(destination):
        os.remove(destination)
    os.symlink(source, destination)
    
    return outs

def main_logic(sr1, sr2, lr, 
               assembly_dir, draft_dir, polish_dir, 
               threads, basecaller, genome_size, genome_length, cov_threshold):
    outs = []
    
    # CALCULATE COVERAGE
    gen_coverage, seqkit_out = calculate_coverage(lr, genome_length)
    outs.append(seqkit_out)

    # SET UNICYCLER STATUS
    unicycler_status = "ok"

    # RUN UNICYCLER IF COVERAGE IS SHALLOW
    if gen_coverage <= cov_threshold:
        unicycler_out, unicycler_status = run_unicycler(sr1, sr2, lr, threads, assembly_dir, draft_dir, polish_dir)
        outs.append(unicycler_out)
        
        if unicycler_status == "error":
            print("The genome coverage is sparse, but Unicycler has crashed.")
        else:
            print("The genome coverage is sparse, Unicycler has been chosen.\nEmpty draft and polish dirs have been created.")
    
    # RUN FMP IF UNICYCLER BROKE OR COVERAGE IS SPARSE
    if gen_coverage > cov_threshold or unicycler_status == "error":
        fmp_outs = run_fmp(sr1, sr2, lr, threads, assembly_dir, draft_dir, polish_dir, basecaller, genome_size, gen_coverage)
        outs.extend(fmp_outs)
    
    return outs

def run_snakemake():
    # DEFINE INPUTS/OUTPUTS
    SHORT_READS_1 = snakemake.input["short_reads_1"]
    SHORT_READS_2 = snakemake.input["short_reads_2"]
    LONG_READS = snakemake.input["long_reads"]
    ASSEMBLY_DIR = snakemake.output["assembly_dir"]
    DRAFT_DIR = snakemake.output["draft_dir"]
    POLISH_DIR = snakemake.output["polish_dir"]

    # DEFINE PARAMS
    THREADS = snakemake.threads
    BASECALLER = snakemake.params["basecaller"]
    GENOME_SIZE = snakemake.params["genome_size"]
    GENOME_LENGTH = snakemake.params["genome_length"]  # expected genome length
    COV_THRESHOLD = snakemake.params["cov_threshold"]  # threshold value to choose assembler, < 20x is good for Uni

    # OPEN LOG
    with open(snakemake.log[0], "w") as f:
        sys.stderr = sys.stdout = f
        
        outs = main_logic(SHORT_READS_1, SHORT_READS_2, LONG_READS,
                          ASSEMBLY_DIR, DRAFT_DIR, POLISH_DIR,
                          THREADS, BASECALLER, GENOME_SIZE, GENOME_LENGTH, COV_THRESHOLD)

        # PRINT STDOUT/STDERR TO LOG
        for out in outs:
            print(out.stdout, out.stderr)

if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake()
    else:
        # For importing in tests
        pass

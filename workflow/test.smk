from snakemake.io import directory, expand

rule all:
    input:
        expand("tests/final/{strain}_all.done", strain=config['strains'])


rule adaptive_hybrid_assembly:
    input:
        short_reads_1 = "tests/data_filtered/{strain}/Illumina/{strain}_1.fq.gz",
        short_reads_2 = "tests/data_filtered/{strain}/Illumina/{strain}_2.fq.gz",
        long_reads = "tests/data_filtered/{strain}/Nanopore/{strain}_all.fastq.gz"
    output:
        assembly_dir = directory("tests/assemblies/{strain}"),
        draft_dir = directory("tests/drafts/{strain}"),
        polish_dir = directory("tests/polished/{strain}")
    threads: 18
    message:
        "executing assembly script with {threads} threads on {wildcards.strain} reads"
    log:
        "tests/logs/{strain}_assembly.log"
    conda: "envs/hybrid_assembly.yaml"
    params: basecaller="r941_min_fast_g507", genome_size="5m", coverage=1
    script:
        "scripts/adaptive_hybrid_assembly.py"

rule final:
    input:
        ass_dir="tests/assemblies/{strain}",
        draft_dir="tests/drafts/{strain}",
        polish_dir="tests/polished/{strain}"
    output: touch("tests/final/{strain}_all.done")
    shell: "echo 'DONE'"


onsuccess:
    print("Workflow finished, no errors")
# Flye-Medaka-Plypolish
from snakemake.io import touch, directory, temp, expand


# Rule to join together all inputs and outputs
rule all:
    input:
        expand("results/final/{strain}_all.done", strain=config['strains'])

rule filter_nanopore:
    input: "resources/data_raw/{strain}/Nanopore/{strain}_all.fastq.gz"
    output: "results/data_filtered/{strain}/Nanopore/{strain}_all.fastq.gz"
    message: "executing filtlong with {threads} threads on {wildcards.strain} long reads"
    log: "results/logs/{strain}_filtlong.log"
    conda: "envs/filtlong.yaml"
    threads: 18
    params: min_len=3000
    shell: "filtlong --min_length {params.min_len} {input} 2> {log} | pigz -c -p {threads} > {output}"

rule flye:
    input: "resources/data_filtered/{strain}/Nanopore/{strain}_all.fastq.gz"
    output: directory("results/flye/{strain}")
    message: "executing Flye with {threads} threads on {wildcards.strain} filtered long reads"
    log: "results/logs/{strain}_flye.log"
    threads: 14
    conda: "envs/flye.yaml"
    params: genome_size="5m", coverage=50
    shell: "flye --nano-raw {input}  --treads {threads} --out-dir {output} -g {params.genome_size} --asm-coverage {params.coverage}"

rule medaka:
    input: 
        reads="resources/data_filtered/{strain}/Nanopore/{strain}_all.fastq.gz", 
        assembly="results/flye/{strain}"
    output: directory("results/medaka/{strain}")
    message: "executing Medaka with {threads} threads on {wildcards.strain} filtered long reads and Flye assembly"
    log: "results/logs/{strain}_medaka.log"
    threads: 14
    conda: "envs/medaka.yaml"
    params: basecaller="r941_min_fast_g507"
    shell: "medaka_consensus -i {input.reads} -d {input.assembly} -o {output} -t {threads} -m {params.basecaller}"

rule map_reads:
    input: 
        assembly="results/medaka/{strain}",
        r1="resources/data_raw/{strain}/Illumina/renamed/{strain}_1.fq.gz",
        r2="resources/data_raw/{strain}/Illumina/renamed/{strain}_2.fq.gz"
    output: sam1=temp("results/polypolish/alignments_1.sam"), 
            sam2=temp("results/polypolish/alignments_2.sam")
    shell: "bwa index {input.assembly}/consensus.fasta &> {log.index} && "
           "bwa mem -t {threads} -a {input.assembly}/consensus.fasta {input.r1} -o {output.sam1} &> {log.sam1} && "
           "bwa mem -t {threads} -a {input.assembly}/consensus.fasta {input.r2} -o {output.sam2} &> {log.sam12}"
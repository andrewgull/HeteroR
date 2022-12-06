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
    threads: 14
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
    shell: "medaka_consensus -i {input.reads} -d {input.assembly}/assembly.fasta -o {output} -t {threads} -m {params.basecaller}"

rule map_reads:
    input: 
        assembly="results/medaka/{strain}",
        r1="resources/data_raw/{strain}/Illumina/renamed/{strain}_1.fq.gz",
        r2="resources/data_raw/{strain}/Illumina/renamed/{strain}_2.fq.gz"
    output: sam1=temp("results/polypolish/{strain}/alignments_1.sam"), 
            sam2=temp("results/polypolish/{strain}/alignments_2.sam")
    message: "Executing BWA to map short reads onto {wildcards.strain} assembly"
    log: index="results/logs/{strain}_fmp_index.log",
         sam1="results/logs/{strain}_fmp_sam1.log",
         sam2="results/logs/{strain}_fmp_sam2.log"
    threads: 14
    conda: "envs/bwa.yaml"
    shell: "bwa index {input.assembly}/consensus.fasta &> {log.index} && "
           "bwa mem -t {threads} -a {input.assembly}/consensus.fasta {input.r1} -o {output.sam1} &> {log.sam1} && "
           "bwa mem -t {threads} -a {input.assembly}/consensus.fasta {input.r2} -o {output.sam2} &> {log.sam2}"

rule insert_filter:
    input: sam1="results/polypolish/{strain}/alignments_1.sam", sam2="results/polypolish/{strain}/alignments_2.sam"
    output: filt1=temp("results/polypolish/{strain}/filtered_1.sam"), filt2=temp("results/polypolish/{strain}/filtered_2.sam")
    message: "Executing Polypolish inser filter script on {wildcards.strain} SAM files"
    log: "results/logs/{strain}_insert_filter.log"
    threads: 14
    conda: "envs/polypolish.yaml"
    shell: "polypolish_insert_filter.py --in1 {input.sam1} --in2 {input.sam2} --out1 {output.filt1} --out2 {output.filt2}"

rule polypolish:
    input: assembly="results/medaka/{strain}",
           filt1="results/polypolish/{strain}/filtered_1.sam", 
           filt2="results/polypolish/{strain}/filtered_2.sam"
    output: "results/polypolish/{strain}/assembly_fmp.fasta"
    message: "Executing polypolish on {wildcards.strain} draft assembly"
    log: "results/logs/{strain}_polypolish.log"
    threads: 14
    conda: "envs/polypolish.yaml"
    shell: "polypolish {input.assembly} {input.filt1} {input.filt2} > {output}"

rule final:
    input: "results/polypolish/{strain}/assembly_fmp.fasta"
    output: touch("results/final/FMP/{strain}_all.done")
    shell: "echo 'DONE'"

onsuccess:
    print("Workflow finished, no errors")
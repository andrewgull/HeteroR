from snakemake.io import touch, directory, temp, expand


# Rule to join together all inputs and outputs
rule all:
    input:
        expand("results/mutants/final/{parent}_all.done", parent=config['parents'])
        #mutants = expand("results/mutants/final/{mutant}_all.done", mutant=config['mutants'])

rule snippy:
    input:
        r1 = "results/data_filtered/{parent}/Illumina/mutants/{parent}m_1.fq.gz",
        r2 = "results/data_filtered/{parent}/Illumina/mutants/{parent}m_2.fq.gz",
        ref = "results/annotations/{parent}/prokka/{parent}_genomic.gbk"
    output:
        directory("results/snippy/{parent}")
    threads: 10
    message: "executing SNIPPY with {threads} threads on {wildcards.parent} assembly"
    log: "results/logs/{parent}_snippy.log"
    conda: "snippy-env"
    shell:
        "snippy --cpus {threads} --outdir {output} --ref {input.ref} --R1 {input.r1} --R2 {input.r2} > {log}"

rule final:
    input:
        snippy="results/snippy/{parent}"
    output: 
        touch("results/mutants/final/{parent}_all.done")
        #touch("results/mutants/final/{mutant}_all.done")
    shell: "echo 'DONE'"

onsuccess:
    print("Workflow finished, no errors")
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


rule make_reference:
    input: "results/annotations/{parent}/prokka/{parent}_genomic.gff"
    output: "results/variants/{parent}/reference.fasta"
    shell: "sed -n '/##FASTA/,${p}' {input} | sed '1d' > {output}"

rule mapping:
    input: r1 = "results/data_filtered/{parent}/Illumina/mutants/{parent}m_1.fq.gz",
           r2 = "results/data_filtered/{parent}/Illumina/mutants/{parent}m_2.fq.gz",
           ref = "results/variants/{parent}/reference.fasta"
    output: temp("results/variants/{parent}/mutant_mapped.sam")
    threads: 10
    message: ""
    log: "results/logs/variants_mapping.log"
    conda: "varcalling-env"
    shell: "bowtie2-build {input.ref} index && bowtie2 -p {threads} -x index -1 {input.r1} -2 {input.r2} -S {output} &> {log}"

rule sorting:
    input: "results/variants/{parent}/mutant_mapped.sam"
    output: "results/variants/{parent}/mutant_mapped.bam"
    threads: 10
    message: ""
    log: "results/logs/variants_sorting.log"
    conda: "varcalling-env"
    shell: "samtools sort -@ {threads} {input} > {output} 2> {log}"

rule variant_calling:
    input: 
        ref = "results/variants/{parent}/reference.fasta",
        bam = "results/variants/{parent}/mutant_mapped.bam"
    output: "results/variants/{parent}/variants.bcf"
    threads: 10
    message: ""
    log: ""
    conda: "varcalling-env"
    params: max_depth=800, ploidy=1, prior=5.0e-10
    shell: "bcftools mpileup --threads {threads} --max-depth {params.max_depth} -f {input.ref} -Ou {input.bam} | bcftools call --threads {threads} --ploidy {params.ploidy} -mv -Ob --prior {params.prior} -o {output}"

rule variant_filtering:
    input: "results/variants/{parent}/variants.bcf"
    output: "results/variants/{parent}/variants_filtered.bcf"
    threads: 10
    message:
    log: ""
    conda: "varcalling-env"
    params: dist = 3, qual = 30, depth = 20
    shell: "bcftools filter -g{params.dist} -i 'QUAL>{params.qual} && DP>{params.depth}' -Ob {input} > {output}"

rule variant_annotation:
    input: gff = "results/annotations/{parent}/prokka/{parent}_genomic.gff",
           bcf = "results/variants/{parent}/variants_filtered.bcf"
    output: gff_clean = "results/variants/{parent}/{parent}_genomic_clean.gff",
            vcf = "results/variants/{parent}/variants_filtered.vcf",
            gff_annotated = ""results/variants/{parent}/genes_with_mutations.tsv""
    shell: "bcftools view {input.bcf} > {output.vcf} && 
           "sed '/##FASTA/,$d' {input.gff} > {output.gff_clean} && " 
           "bedtools annotate -i {output.gff_clean} -files {output.vcf} > {output.gff_annotated}"

rule final:
    input:
        snippy="results/snippy/{parent}"
    output: 
        touch("results/mutants/final/{parent}_all.done")
    shell: "echo 'DONE'"

onsuccess:
    print("Workflow finished, no errors")
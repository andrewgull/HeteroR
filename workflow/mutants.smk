from snakemake.io import touch, directory, temp, expand


# Rule to join together all inputs and outputs
rule all:
    input:
        expand("results/mutants/final/{parent}_all.done", parent=config['parents'])
        #mutants = expand("results/mutants/final/{mutant}_all.done", mutant=config['mutants'])

rule make_reference:
    input: "results/annotations/{parent}/prokka/{parent}_genomic.gff"
    output: "results/variants/{parent}/reference.fasta"
    shell: "sed -n '/##FASTA/,${{p}}' {input} | sed '1d' > {output}"

rule mapping:
    input: r1 = "results/data_filtered/{parent}/Illumina/mutants/{parent}m_1.fq.gz",
           r2 = "results/data_filtered/{parent}/Illumina/mutants/{parent}m_2.fq.gz",
           ref = "results/variants/{parent}/reference.fasta"
    output: sam = temp("results/variants/{parent}/mutant_mapped.sam")
    threads: 10
    message: "Mapping mutant reads onto {wildcards.parent} genome"
    log: mapping = "results/logs/{parent}_variants_mapping.log",
         index = "results/logs/{parent}_reference_index.log"
    conda: "varcalling-env"
    shell: "bowtie2-build {input.ref} index &> {log.index} && "
           "bowtie2 -p {threads} -x index -1 {input.r1} -2 {input.r2} -S {output.sam} &> {log.mapping}"

rule sorting:
    input: "results/variants/{parent}/mutant_mapped.sam"
    output: "results/variants/{parent}/mutant_mapped.bam"
    threads: 10
    message: "Sorting mapped reads of {wildcards.parent} mutant"
    log: "results/logs/{parent}_variants_sorting.log"
    conda: "varcalling-env"
    shell: "samtools sort -@ {threads} {input} > {output} 2> {log}"

rule variant_calling:
    input: 
        ref = "results/variants/{parent}/reference.fasta",
        bam = "results/variants/{parent}/mutant_mapped.bam"
    output: "results/variants/{parent}/variants.bcf"
    threads: 10
    message: "Calling varinats in {wildcards.parent} mutant"
    log: call = "results/logs/{parent}_variant_calling.log",
         mpileup = "results/logs/{parent}_mpileup.log"
    conda: "varcalling-env"
    params: max_depth=config["max_depth"], ploidy=config["ploidy"], prior=config["prior"]
    shell: "bcftools mpileup --threads {threads} -d {params.max_depth} -f {input.ref} -Ou {input.bam} 2> {log.mpileup} | "
           "bcftools call --threads {threads} --ploidy {params.ploidy} -mv -Ob -P {params.prior} -o {output} 2> {log.call}"

rule variant_filtering:
    input: "results/variants/{parent}/variants.bcf"
    output: "results/variants/{parent}/variants_filtered.bcf"
    threads: 10
    message: "Filtering variants in {wildcards.parent} mutant"
    log: "results/logs/{parent}_variant_filtering.log"
    conda: "varcalling-env"
    params: dist=config["indel_dist"], qual=config["quality"], depth=config["depth"]
    shell: "bcftools filter -g{params.dist} -i 'QUAL>{params.qual} && DP>{params.depth}' -Ob {input} > {output} 2> {log}"

rule variant_annotation:
    input: gff = "results/annotations/{parent}/prokka/{parent}_genomic.gff",
           bcf = "results/variants/{parent}/variants_filtered.bcf"
    output: gff_clean = "results/variants/{parent}/{parent}_genomic_clean.gff",
            vcf = "results/variants/{parent}/variants_filtered.vcf",
            gff_annotated = "results/variants/{parent}/genes_with_mutations.tsv"
    threads: 10
    message: "Annotating variants in {wildcards.parent} mutant"
    log: view = "results/logs/{parent}_bcftools_view.log",
         sed = "results/logs/{parent}_sed_clean_gff.log",
         annotate = "results/logs/{parent}_variant_annotation.log"
    conda: "varcalling-env"
    shell: "bcftools view {input.bcf} > {output.vcf}  2> {log.view} && "
           "sed '/##FASTA/,$d' {input.gff} > {output.gff_clean} 2> {log.sed} && " 
           "bedtools annotate -i {output.gff_clean} -files {output.vcf} > {output.gff_annotated} 2> {log.annotate}"

rule final:
    input:
        gff_annotated= "results/variants/{parent}/genes_with_mutations.tsv"
    output: 
        touch("results/mutants/final/{parent}_all.done")
    shell: "echo 'DONE'"

onsuccess:
    print("Workflow finished, no errors")
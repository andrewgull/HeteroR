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
    conda: "envs/var_calling.yaml"
    shell: "bowtie2-build {input.ref} index &> {log.index} && "
           "bowtie2 -p {threads} -x index -1 {input.r1} -2 {input.r2} -S {output.sam} &> {log.mapping}"

rule sorting:
    input: "results/variants/{parent}/mutant_mapped.sam"
    output: "results/variants/{parent}/mutant_mapped.bam"
    threads: 10
    message: "Sorting mapped reads of {wildcards.parent} mutant"
    log: "results/logs/{parent}_variants_sorting.log"
    conda: "envs/var_calling.yaml"
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
    conda: "envs/var_calling.yaml"
    params: max_depth=config["max_depth"], ploidy=config["ploidy"], prior=config["prior"]
    shell: "bcftools mpileup --threads {threads} -d {params.max_depth} -f {input.ref} -Ou {input.bam} 2> {log.mpileup} | "
           "bcftools call --threads {threads} --ploidy {params.ploidy} -mv -Ob -P {params.prior} -o {output} 2> {log.call}"

rule variant_filtering:
    input: "results/variants/{parent}/variants.bcf"
    output: "results/variants/{parent}/variants_filtered.bcf"
    threads: 10
    message: "Filtering variants in {wildcards.parent} mutant"
    log: "results/logs/{parent}_variant_filtering.log"
    conda: "envs/var_calling.yaml"
    params: dist=config["indel_dist"], qual=config["quality"], depth=config["depth"]
    shell: "bcftools filter -g{params.dist} -i 'QUAL>{params.qual} && DP>{params.depth}' -Ob {input} > {output} 2> {log}"

rule variant_annotation:
    input: gff = "results/annotations/{parent}/prokka/{parent}_genomic.gff",
           bcf = "results/variants/{parent}/variants_filtered.bcf"
    output: gff_clean = "results/variants/{parent}/{parent}_genomic_clean.gff",
            vcf = "results/variants/{parent}/variants_filtered.vcf",
            gff_annotated = "results/variants/{parent}/{parent}_annotated_variants.gff"
    threads: 10
    message: "Annotating variants in {wildcards.parent} mutant"
    log: view = "results/logs/{parent}_bcftools_view.log",
         sed = "results/logs/{parent}_sed_clean_gff.log",
         annotate = "results/logs/{parent}_variant_annotation.log"
    conda: "envs/var_calling.yaml"
    shell: "bcftools view {input.bcf} > {output.vcf}  2> {log.view} && "
           "sed '/##FASTA/,$d' {input.gff} > {output.gff_clean} 2> {log.sed} && " 
           "bedtools annotate -i {output.gff_clean} -files {output.vcf} > {output.gff_annotated} 2> {log.annotate}"

rule filter_annotation:
    input: script = "workflow/scripts/filter_variant_annotations.R",
           gff = "results/variants/{parent}/{parent}_annotated_variants.gff"
    output: "results/variants/{parent}/{parent}_genes_with_variants.tsv"
    message: "Filtering annotated GFF/VCF in {wildcards.parent} mutant"
    log: "results/logs/{parent}_filter_variant_annotation.log"
    conda: "envs/rscripts.yaml"
    shell: "Rscript {input.script} -i {input.gff} -o {output} &> {log}"

rule depth:
    input: "results/variants/{parent}/mutant_mapped.bam"
    output: "results/amplifications/{parent}/depth.tsv.gz"
    threads: 10
    message: "Calculating depth in {wildcards.parent} BAM file"
    log: "results/logs/{parent}_bam_depth.log"
    conda: "envs/var_calling.yaml"
    shell: "samtools depth -@ {threads} {input} | gzip -c > {output} 2> {log}"

rule amplified_regions:
    input: script = "workflow/scripts/find_amplifications.R",
           depth = "results/amplifications/{parent}/depth.tsv.gz"
    output: bed = "results/amplifications/{parent}/amplifications_windows.bed",
            plot = "results/amplifications/{parent}/genome_coverage.png"
    message: "Looking for over-covered windows in {wildcards.parent} mutant genome"
    log: "results/logs/{parent}_amplifications.log"
    conda: "envs/rscripts.yaml"
    params: z = config["z_threshold"], w = config["window_size"]
    shell: "Rscript {input.script} -i {input.depth} -b {output.bed} -l {output.plot} -z {params.z} -w {params.w} &> {log}"

rule final:
    input:
        bed = "results/amplifications/{parent}/amplifications_windows.bed",
        genes_w_snps = "results/variants/{parent}/{parent}_genes_with_variants.tsv"
    output: 
        touch("results/mutants/final/{parent}_all.done")
    shell: "echo 'DONE'"

onsuccess:
    print("Workflow finished, no errors")
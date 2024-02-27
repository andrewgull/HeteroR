from snakemake.io import touch, directory, temp, expand


# Rule to join together all inputs and outputs
rule all:
    input:
        expand("results/mutants/final/{parent}_all.done", parent=config['parents'])
        #mutants = expand("results/mutants/final/{mutant}_all.done", mutant=config['mutants'])

rule trim_reads:
    input: 
        r1 = "resources/data_raw/{parent}/Illumina/mutants/{parent}m_1.fq.gz",
        r2 = "resources/data_raw/{parent}/Illumina/mutants/{parent}m_2.fq.gz"
    output:
        r1 = "results/data_filtered/{parent}/Illumina/mutants/{parent}m_1.fq.gz",
        r2 = "results/data_filtered/{parent}/Illumina/mutants/{parent}m_2.fq.gz"
    threads: 10
    message: "trimming front ends of the {wildcards.parent} reads"
    log: "results/logs/{parent}_mutants_trimming.log"
    conda: "envs/fastp.yaml"
    params: f = 10
    shell: "fastp --in1 {input.r1} --in2 {input.r2} --out1 {output.r1} --out2 {output.r2} --thread {threads} --trim_front1 {params.f} --trim_front2 {params.f} &> {log}" 

rule make_reference:
    input: "results/annotations/{parent}/prokka/{parent}_genomic.gff"
    output: "results/mutants/variants/{parent}/reference.fasta"
    shell: "sed -n '/##FASTA/,${{p}}' {input} | sed '1d' > {output}"

rule mapping_mutant:
    input: r1 = "results/data_filtered/{parent}/Illumina/mutants/{parent}m_1.fq.gz",
           r2 = "results/data_filtered/{parent}/Illumina/mutants/{parent}m_2.fq.gz",
           ref = "results/mutants/variants/{parent}/reference.fasta"
    output: sam = temp("results/mutants/variants/{parent}/mutant_mapped.sam")
    threads: 10
    message: "Mapping mutant reads onto {wildcards.parent} genome"
    log: mapping = "results/logs/{parent}_variants_mapping.log",
         index = "results/logs/{parent}_reference_index.log"
    conda: "envs/var_calling.yaml"
    shell: "bowtie2-build {input.ref} index &> {log.index} && "
           "bowtie2 -p {threads} -x index -1 {input.r1} -2 {input.r2} -S {output.sam} &> {log.mapping}"

rule sorting_mutant:
    input: "results/mutants/variants/{parent}/mutant_mapped.sam"
    output: "results/mutants/variants/{parent}/mutant_mapped.bam"
    threads: 10
    message: "Sorting mapped reads of {wildcards.parent} mutant"
    log: "results/logs/{parent}_mutant_reads_sorting.log"
    conda: "envs/var_calling.yaml"
    shell: "samtools sort -@ {threads} {input} > {output} 2> {log}"

rule mapping_parent:
    # use the same reference file as for the mutant
    input: r1="results/data_filtered/{parent}/Illumina/{parent}_1.fq.gz",
        r2="results/data_filtered/{parent}/Illumina/{parent}_2.fq.gz",
        ref="results/mutants/variants/{parent}/reference.fasta"
    output: sam=temp("results/mutants/copy_number/{parent}/parent_mapped.sam")
    threads: 10
    message: "Mapping parent reads onto {wildcards.parent} reference from GFF file"
    log: mapping="results/logs/{parent}_parent_mapping.log",
        index="results/logs/{parent}_parent_reference_index.log"
    conda: "envs/var_calling.yaml"
    shell: "bowtie2-build {input.ref} index &> {log.index} && "
           "bowtie2 -p {threads} -x index -1 {input.r1} -2 {input.r2} -S {output.sam} &> {log.mapping}"

rule sorting_parent:
    input: "results/mutants/copy_number/{parent}/parent_mapped.sam"
    output: temp("results/mutants/copy_number/{parent}/parent_mapped.bam")
    threads: 10
    message: "Sorting mapped reads of {wildcards.parent} parent"
    log: "results/logs/{parent}_parent_reads_sorting.log"
    conda: "envs/var_calling.yaml"
    shell: "samtools sort -@ {threads} {input} > {output} 2> {log}"

rule depth_parent:
    input: "results/mutants/copy_number/{parent}/parent_mapped.bam"
    output: "results/mutants/copy_number/{parent}/parent_depth.tsv.gz"
    threads: 10
    message: "Calculating depth in {wildcards.parent} parental BAM file"
    log: "results/logs/{parent}_bam_parent_depth.log"
    conda: "envs/var_calling.yaml"
    shell: "samtools depth -@ {threads} {input} | gzip -c > {output} 2> {log}"

rule variant_calling:
    input: 
        ref = "results/mutants/variants/{parent}/reference.fasta",
        bam = "results/mutants/variants/{parent}/mutant_mapped.bam"
    output: "results/mutants/variants/{parent}/variants.bcf"
    threads: 10
    message: "Calling varinats in {wildcards.parent} mutant"
    log: call = "results/logs/{parent}_variant_calling.log",
         mpileup = "results/logs/{parent}_mpileup.log"
    conda: "envs/var_calling.yaml"
    params: max_depth=config["max_depth"], ploidy=config["ploidy"], prior=config["prior"]
    shell: "bcftools mpileup --threads {threads} -d {params.max_depth} -f {input.ref} -Ou {input.bam} 2> {log.mpileup} | "
           "bcftools call --threads {threads} --ploidy {params.ploidy} -mv -Ob -P {params.prior} -o {output} 2> {log.call}"

rule variant_filtering:
    input: "results/mutants/variants/{parent}/variants.bcf"
    output: "results/mutants/variants/{parent}/variants_filtered.bcf"
    threads: 10
    message: "Filtering variants in {wildcards.parent} mutant"
    log: "results/logs/{parent}_variant_filtering.log"
    conda: "envs/var_calling.yaml"
    params: dist=config["indel_dist"], qual=config["quality"], depth=config["depth"]
    shell: "bcftools filter -g{params.dist} -i 'QUAL>{params.qual} && DP>{params.depth}' -Ob {input} > {output} 2> {log}"

rule variant_annotation:
    input: gff = "results/annotations/{parent}/prokka/{parent}_genomic.gff",
           bcf = "results/mutants/variants/{parent}/variants_filtered.bcf"
    output: gff_clean = "results/mutants/variants/{parent}/genomic_clean.gff",
            vcf = "results/mutants/variants/{parent}/variants_filtered.vcf",
            gff_annotated = "results/mutants/variants/{parent}/annotated_variants.gff"
    threads: 10
    message: "Annotating variants in {wildcards.parent} mutant"
    log: view = "results/logs/{parent}_bcftools_view.log",
         sed = "results/logs/{parent}_sed_clean_gff.log",
         annotate = "results/logs/{parent}_variant_annotation.log"
    conda: "envs/var_calling.yaml"
    shell: "bcftools view {input.bcf} > {output.vcf}  2> {log.view} && "
           "sed '/##FASTA/,$d' {input.gff} > {output.gff_clean} 2> {log.sed} && " 
           "bedtools annotate -i {output.gff_clean} -files {output.vcf} > {output.gff_annotated} 2> {log.annotate}"

rule filter_variant_annotation:
    input: script = "workflow/scripts/filter_gff_annotations.R",
           gff = "results/mutants/variants/{parent}/annotated_variants.gff"
    output: "results/mutants/variants/{parent}/genes_with_variants.tsv"
    message: "Filtering annotated GFF/VCF in {wildcards.parent} mutant"
    log: "results/logs/{parent}_filter_variant_annotation.log"
    conda: "envs/rscripts.yaml"
    shell: "Rscript {input.script} -i {input.gff} -o {output} &> {log}"

rule depth_mutant:
    input: "results/mutants/variants/{parent}/mutant_mapped.bam"
    output: "results/mutants/amplifications/{parent}/mutant_depth.tsv.gz"
    threads: 10
    message: "Calculating depth in {wildcards.parent} BAM file"
    log: "results/logs/{parent}_bam_depth.log"
    conda: "envs/var_calling.yaml"
    shell: "samtools depth -@ {threads} {input} | gzip -c > {output} 2> {log}"

rule find_amplified_regions:
    input: script = "workflow/scripts/find_amplifications.R",
           depth = "results/mutants/amplifications/{parent}/mutant_depth.tsv.gz"
    output: bed = "results/mutants/amplifications/{parent}/amplifications_windows.bed",
            plot = "results/mutants/amplifications/{parent}/genome_coverage.png"
    message: "Looking for over-covered windows in {wildcards.parent} mutant genome"
    log: "results/logs/{parent}_amplifications.log"
    conda: "envs/rscripts.yaml"
    params: z = config["z_threshold"], w = config["window_size"]
    shell: "Rscript {input.script} -i {input.depth} -b {output.bed} -l {output.plot} -z {params.z} -w {params.w} &> {log}"

rule merge_amplified_regions:
    input: "results/mutants/amplifications/{parent}/amplifications_windows.bed"
    output: "results/mutants/amplifications/{parent}/amplifications_merged.bed"
    message: "Merging overlapping windows in {wildcards.parent} mutant"
    log: "results/logs/{parent}_merging.log"
    conda: "envs/var_calling.yaml"
    shell: "bedtools merge -i {input} > {output} 2> {log}"

rule annotate_amplified_regions:
    input: bed = "results/mutants/amplifications/{parent}/amplifications_merged.bed",
           gff = "results/mutants/variants/{parent}/genomic_clean.gff"
    output: "results/mutants/amplifications/{parent}/amplifications_annotated.gff"
    message: "Annotating merged amplified regions in {wildcards.parent} mutants"
    log: "results/logs/{parent}_annotate_amplifications.log"
    conda: "envs/var_calling.yaml"
    shell: "touch {output}; bedtools annotate -i {input.gff} -files {input.bed} | grep -v '0.000000' 1>> {output} 2> {log}"

rule filter_annotated_amplified_regions:
    input: script = "workflow/scripts/filter_gff_annotations.R",
           gff = "results/mutants/amplifications/{parent}/amplifications_annotated.gff"
    output: "results/mutants/amplifications/{parent}/amplifications_annotated_filtered.tsv"
    message: "Filtering annotated amplifications in {wildcards.parent} mutant"
    log: "results/logs/{parent}_filter_amplification_annotation.log"
    conda: "envs/rscripts.yaml"
    shell: "Rscript {input.script} -i {input.gff} -o {output} &> {log}"

rule relative_coverage_mutant:
    input: depth = "results/mutants/amplifications/{parent}/mutant_depth.tsv.gz", # this one is from mapping of mutant reads
           script = "workflow/scripts/relative_coverage.R"
    output: "results/mutants/copy_number/{parent}/relative_coverage_mutant.tsv"
    message: "Calculating relative coverage on {wildcards.parent} mutants"
    log: "results/logs/{parent}/relative_coverage_mutants.log"
    conda: "envs/rscripts.yaml"
    params: min_len = config["min_contig_len"]
    shell: "Rscript {input.script} -i {input.depth} -o {output} -m {params.min_len} -l mutant -s {wildcards.parent} &> {log}"

rule relative_coverage_parent:
    input: depth = "results/mutants/copy_number/{parent}/parent_depth.tsv.gz", # this one is from mapping of parental reads
           script = "workflow/scripts/relative_coverage.R"
    output: "results/mutants/copy_number/{parent}/relative_coverage_parent.tsv"
    message: "Calculating relative coverage on {wildcards.parent} parent"
    log: "results/logs/{parent}/relative_coverage_parent.log"
    conda: "envs/rscripts.yaml"
    params: min_len = config["min_contig_len"]
    shell: "Rscript {input.script} -i {input.depth} -o {output} -m {params.min_len} -l parent -s {wildcards.parent} &> {log}"

rule final:
    input:
        bed = "results/mutants/amplifications/{parent}/amplifications_windows.bed",
        genes_w_snps = "results/mutants/variants/{parent}/genes_with_variants.tsv",
        amplifications = "results/mutants/amplifications/{parent}/amplifications_annotated_filtered.tsv",
        rel_cov_mut = "results/mutants/copy_number/{parent}/relative_coverage_mutant.tsv",
        rel_cov_parent = "results/mutants/copy_number/{parent}/relative_coverage_parent.tsv"
    output: 
        touch("results/mutants/final/{parent}_all.done")
    shell: "echo 'DONE'"

onsuccess:
    print("Workflow finished, no errors")
configfile: "config/mutants.yaml"


parents = pd.read_csv(config["parents"], dtype={"parents": str})

mutants_path = config["mutants_path"]
parents_path = config["parents_path"]


rule trim_mutants:
    input:
        r1=lambda wildcards: f"{mutants_path}/{wildcards.parent}_1.fq.gz",
        r2=lambda wildcards: f"{mutants_path}/{wildcards.parent}_2.fq.gz",
    output:
        r1="results/data_filtered/{parent}/short/mutants/{parent}m_1.fq.gz",
        r2="results/data_filtered/{parent}/short/mutants/{parent}m_2.fq.gz",
    threads: 10
    message:
        "trimming front ends of the {wildcards.parent} reads"
    log:
        "results/logs/{parent}_mutants_trimming.log",
    conda:
        "../envs/fastp.yaml"
    container:
        config.get("default_container", None)
    params:
        f=config.get("trim_front"),
        adapter1=config.get("adapter1"),
        adapter2=config.get("adapter2"),
    shell:
        "fastp --in1 {input.r1} --in2 {input.r2} --out1 {output.r1} --out2 {output.r2} --thread {threads} "
        "--trim_front1 {params.f} --trim_front2 {params.f} --adapter_sequence {params.adapter1} "
        "--adapter_sequence_r2 {params.adapter2} &> {log}"


if config.get("trim_parents", False):

    rule trim_parents:
        input:
            r1=lambda wildcards: f"{parents_path}/{wildcards.parent}_1.fq.gz",
            r2=lambda wildcards: f"{parents_path}/{wildcards.parent}_2.fq.gz",
        output:
            r1="results/data_filtered/{parent}/short/{parent}_1.fq.gz",
            r2="results/data_filtered/{parent}/short/{parent}_2.fq.gz",
        threads: 10
        message:
            "trimming front ends of the {wildcards.parent} reads"
        log:
            "results/logs/{parent}_parents_trimming.log",
        conda:
            "../envs/fastp.yaml"
        params:
            f=config.get("trim_front"),
            adapter1=config.get("adapter1"),
            adapter2=config.get("adapter2"),
        shell:
            "fastp --in1 {input.r1} --in2 {input.r2} --out1 {output.r1} --out2 {output.r2} --thread {threads} "
            "--trim_front1 {params.f} --trim_front2 {params.f} --adapter_sequence {params.adapter1} "
            "--adapter_sequence_r2 {params.adapter2} &> {log}"


rule create_links:
    # required for ISmapper - it recognizes files only if 'fastq.gz' is in the name
    input:
        r1="results/data_filtered/{parent}/short/mutants/{parent}m_1.fq.gz",
        r2="results/data_filtered/{parent}/short/mutants/{parent}m_2.fq.gz",
    output:
        r1="results/data_filtered/{parent}/short/mutants/{parent}_1.fastq.gz",
        r2="results/data_filtered/{parent}/short/mutants/{parent}_2.fastq.gz",
    log:
        "../results/logs/{parent}_links.log",
    shell:
        "ln -s {input.r1} {output.r2} && ln -s {input.r2} {output.r2} 2> {log}"


rule make_reference:
    input:
        "results/annotations/{parent}/prokka/{parent}_genomic.gff",
    output:
        "results/mutants/variants/{parent}/reference.fasta",
    log:
        "results/logs/{parent}_ref.log",
    shell:
        "sed -n '/##FASTA/,${{p}}' {input} | sed '1d' > {output} 2> {log}"


rule mapping_mutant:
    input:
        r1="results/data_filtered/{parent}/short/mutants/{parent}m_1.fq.gz",
        r2="results/data_filtered/{parent}/short/mutants/{parent}m_2.fq.gz",
        ref="results/mutants/variants/{parent}/reference.fasta",
    output:
        sam=temp("results/mutants/variants/{parent}/mutant_mapped.sam"),
    threads: 10
    message:
        "Mapping mutant reads onto {wildcards.parent} genome"
    log:
        mapping="results/logs/{parent}_variants_mapping.log",
        index="results/logs/{parent}_reference_index.log",
    conda:
        "../envs/varcalling.yaml"
    container:
        config.get("default_container", None)
    shell:
        "bowtie2-build {input.ref} index &> {log.index} && "
        "bowtie2 -p {threads} -x index -1 {input.r1} -2 {input.r2} -S {output.sam} &> {log.mapping}"


rule sorting_mutant:
    input:
        "results/mutants/variants/{parent}/mutant_mapped.sam",
    output:
        "results/mutants/variants/{parent}/mutant_mapped.bam",
    threads: 10
    message:
        "Sorting mapped reads of {wildcards.parent} mutant"
    log:
        "results/logs/{parent}_mutant_reads_sorting.log",
    conda:
        "../envs/varcalling.yaml"
    container:
        config.get("default_container", None)
    shell:
        "samtools sort -@ {threads} {input} > {output} 2> {log}"


rule mapping_parent:
    # use the same reference file as for the mutant
    input:
        r1="results/data_filtered/{parent}/short/{parent}_1.fq.gz",
        r2="results/data_filtered/{parent}/short/{parent}_2.fq.gz",
        ref="results/mutants/variants/{parent}/reference.fasta",
    output:
        sam=temp("results/mutants/copy_number/{parent}/parent_mapped.sam"),
    threads: 10
    message:
        "Mapping parent reads onto {wildcards.parent} reference from GFF file"
    log:
        mapping="results/logs/{parent}_parent_mapping.log",
        index="results/logs/{parent}_parent_reference_index.log",
    conda:
        "../envs/varcalling.yaml"
    container:
        config.get("default_container", None)
    shell:
        "bowtie2-build {input.ref} index &> {log.index} && "
        "bowtie2 -p {threads} -x index -1 {input.r1} -2 {input.r2} -S {output.sam} &> {log.mapping}"


rule sorting_parent:
    input:
        "results/mutants/copy_number/{parent}/parent_mapped.sam",
    output:
        temp("results/mutants/copy_number/{parent}/parent_mapped.bam"),
    threads: 10
    message:
        "Sorting mapped reads of {wildcards.parent} parent"
    log:
        "results/logs/{parent}_parent_reads_sorting.log",
    conda:
        "../envs/varcalling.yaml"
    container:
        config.get("default_container", None)
    shell:
        "samtools sort -@ {threads} {input} > {output} 2> {log}"


rule depth_parent:
    input:
        "results/mutants/copy_number/{parent}/parent_mapped.bam",
    output:
        "results/mutants/copy_number/{parent}/parent_depth.tsv.gz",
    threads: 10
    message:
        "Calculating depth in {wildcards.parent} parental BAM file"
    log:
        "results/logs/{parent}_bam_parent_depth.log",
    conda:
        "../envs/varcalling.yaml"
    container:
        config.get("default_container", None)
    shell:
        "samtools depth -@ {threads} {input} | gzip -c > {output} 2> {log}"


rule variant_calling:
    input:
        ref="results/mutants/variants/{parent}/reference.fasta",
        bam="results/mutants/variants/{parent}/mutant_mapped.bam",
    output:
        "results/mutants/variants/{parent}/variants.bcf",
    threads: 10
    message:
        "Calling varinats in {wildcards.parent} mutant"
    log:
        call="results/logs/{parent}_variant_calling.log",
        mpileup="results/logs/{parent}_mpileup.log",
    conda:
        "../envs/varcalling.yaml"
    container:
        config.get("default_container", None)
    params:
        max_depth=config.get("max_depth"),
        ploidy=config.get("ploidy"),
        prior=config.get("prior"),
    shell:
        "bcftools mpileup --threads {threads} -d {params.max_depth} -f {input.ref} -Ou {input.bam} 2> {log.mpileup} | "
        "bcftools call --threads {threads} --ploidy {params.ploidy} -mv -Ob -P {params.prior} -o {output} 2> {log.call}"


rule variant_filtering:
    input:
        "results/mutants/variants/{parent}/variants.bcf",
    output:
        "results/mutants/variants/{parent}/variants_filtered.bcf",
    threads: 10
    message:
        "Filtering variants in {wildcards.parent} mutant"
    log:
        "results/logs/{parent}_variant_filtering.log",
    conda:
        "../envs/varcalling.yaml"
    container:
        config.get("default_container", None)
    params:
        dist=config.get("indel_dist"),
        qual=config.get("quality"),
        depth=config.get("depth"),
    shell:
        "bcftools filter -g{params.dist} -i 'QUAL>{params.qual} && DP>{params.depth}' -Ob {input} > {output} 2> {log}"


rule variant_annotation:
    input:
        gff="results/annotations/{parent}/prokka/{parent}_genomic.gff",
        bcf="results/mutants/variants/{parent}/variants_filtered.bcf",
    output:
        gff_clean="results/mutants/variants/{parent}/genomic_clean.gff",
        vcf="results/mutants/variants/{parent}/variants_filtered.vcf",
        gff_annotated="results/mutants/variants/{parent}/annotated_variants.gff",
    threads: 10
    message:
        "Annotating variants in {wildcards.parent} mutant"
    log:
        view="results/logs/{parent}_bcftools_view.log",
        sed="results/logs/{parent}_sed_clean_gff.log",
        annotate="results/logs/{parent}_variant_annotation.log",
    conda:
        "../envs/varcalling.yaml"
    container:
        config.get("default_container", None)
    shell:
        "bcftools view {input.bcf} > {output.vcf}  2> {log.view} && "
        "sed '/##FASTA/,$d' {input.gff} > {output.gff_clean} 2> {log.sed} && "
        "bedtools annotate -i {output.gff_clean} -files {output.vcf} > {output.gff_annotated} 2> {log.annotate}"


rule filter_variant_annotation:
    input:
        "results/mutants/variants/{parent}/annotated_variants.gff",
    output:
        "results/mutants/variants/{parent}/genes_with_variants.tsv",
    message:
        "Filtering annotated GFF/VCF in {wildcards.parent} mutant"
    log:
        "results/logs/{parent}_filter_gff_annotations.log",
    conda:
        "../envs/rscripts.yaml"
    container:
        config.get("rscripts_container", None)
    script:
        "../scripts/filter_gff_annotations.R"


rule depth_mutant:
    input:
        "results/mutants/variants/{parent}/mutant_mapped.bam",
    output:
        "results/mutants/amplifications/{parent}/mutant_depth.tsv.gz",
    threads: 10
    message:
        "Calculating depth in {wildcards.parent} BAM file"
    log:
        "results/logs/{parent}_bam_depth.log",
    conda:
        "../envs/varcalling.yaml"
    container:
        config.get("default_container", None)
    shell:
        "samtools depth -@ {threads} {input} | gzip -c > {output} 2> {log}"


rule find_amplified_regions:
    input:
        "results/mutants/amplifications/{parent}/mutant_depth.tsv.gz",
    output:
        bed="results/mutants/amplifications/{parent}/amplifications_windows.bed",
        plot="results/mutants/amplifications/{parent}/genome_coverage.png",
    message:
        "Looking for windows with increased coverage in {wildcards.parent} mutant genome"
    log:
        "results/logs/{parent}_amplifications.log",
    conda:
        "../envs/rscripts.yaml"
    container:
        config.get("rscripts_container", None)
    params:
        z=config.get("z_threshold"),
        w=config.get("window_size"),
    script:
        "../scripts/find_amplifications.R"


rule merge_amplified_regions:
    input:
        "results/mutants/amplifications/{parent}/amplifications_windows.bed",
    output:
        "results/mutants/amplifications/{parent}/amplifications_merged.bed",
    message:
        "Merging overlapping windows in {wildcards.parent} mutant"
    log:
        "results/logs/{parent}_merging.log",
    conda:
        "../envs/varcalling.yaml"
    container:
        config.get("default_container", None)
    shell:
        "bedtools merge -i {input} > {output} 2> {log}"


rule annotate_amplified_regions:
    input:
        bed="results/mutants/amplifications/{parent}/amplifications_merged.bed",
        gff="results/mutants/variants/{parent}/genomic_clean.gff",
    output:
        "results/mutants/amplifications/{parent}/amplifications_annotated.gff",
    message:
        "Annotating merged amplified regions in {wildcards.parent} mutants"
    log:
        "results/logs/{parent}_annotate_amplifications.log",
    conda:
        "../envs/varcalling.yaml"
    container:
        config.get("default_container", None)
    shell:
        "touch {output}; bedtools annotate -i {input.gff} -files {input.bed} | grep -v '0.000000' 1>> {output} 2> {log}"


rule filter_annotated_amplified_regions:
    input:
        "results/mutants/amplifications/{parent}/amplifications_annotated.gff",
    output:
        "results/mutants/amplifications/{parent}/amplifications_annotated_filtered.tsv",
    message:
        "Filtering annotated amplifications in {wildcards.parent} mutant"
    log:
        "results/logs/{parent}_filter_amplification_annotation.log",
    conda:
        "../envs/rscripts.yaml"
    container:
        config.get("rscripts_container", None)
    script:
        "../scripts/filter_gff_annotations.R"


rule relative_coverage_mutant:
    input:
        depth="results/mutants/amplifications/{parent}/mutant_depth.tsv.gz",  # this one is from mapping of mutant reads
        ref="results/mutants/variants/{parent}/reference.fasta",
    output:
        rel_cov="results/mutants/copy_number/{parent}/relative_coverage_mutant.tsv",
    message:
        "Calculating relative coverage on {wildcards.parent} mutants"
    log:
        "results/logs/{parent}/relative_coverage_mutants.log",
    conda:
        "../envs/biostrings.yaml"
    container:
        config.get("biostrings_container", None)
    params:
        min_len=config.get("min_contig_len"),
        label="mutant",
    script:
        "../scripts/relative_coverage.R"


rule relative_coverage_parent:
    input:
        depth="results/mutants/copy_number/{parent}/parent_depth.tsv.gz",  # this one is from mapping of parental reads
        ref="results/mutants/variants/{parent}/reference.fasta",
    output:
        rel_cov="results/mutants/copy_number/{parent}/relative_coverage_parent.tsv",
    message:
        "Calculating relative coverage on {wildcards.parent} parent"
    log:
        "results/logs/{parent}/relative_coverage_parent.log",
    conda:
        "../envs/biostrings.yaml"
    container:
        config.get("biostrings_container", None)
    params:
        min_len=config.get("min_contig_len"),
        label="parent",
    script:
        "../scripts/relative_coverage.R"


rule collect_all_IS:
    input:
        expand(
            "results/isescan/{parent}/regions/regions_joined_final.fasta.is.fna",
            parent=parents["parents"],
        ),
    output:
        "results/mutants/ismapper/query_collection.fasta",
    log:
        "results/logs/collect_all_is.log",
    shell:
        "cat {input} > {output} 2> {log}"


rule find_best_IS:
    input:
        # TODO: must be file
        "results/isescan/",  # trailing / is important!
    output:
        "results/mutants/ismapper/best_IS_representatives.tsv",
    message:
        "Looking for best representatives of each IS family"
    log:
        "results/logs/mutants_best_IS.log",
    conda:
        "../envs/rscripts.yaml"
    container:
        config.get("rscripts_container", None)
    script:
        "../scripts/find_best_IS_examples.R"


rule extract_IS_headers:
    input:
        table="results/mutants/ismapper/best_IS_representatives.tsv",
        collection="results/mutants/ismapper/query_collection.fasta",
    output:
        "results/mutants/ismapper/best_IS_representatives_IDs.txt",
    log:
        "results/logs/is_headers.log",
    shell:
        "while IFS=$'\t' read -r col1 col2 _; do grep $col1 {input.collection} | grep $col2 >> {output}; done < {input.table} && sed -i 's/>//g' {output} 2> {log}"


rule extract_best_IS:
    input:
        id_file="results/mutants/ismapper/best_IS_representatives_IDs.txt",
        collection="results/mutants/ismapper/query_collection.fasta",
    output:
        "results/mutants/ismapper/best_IS_from_each_family.fasta",
    log:
        "results/logs/best_is.log",
    container:
        config.get("default_container", None)
    shell:
        "seqkit grep -n -f {input.id_file} {input.collection} -o {output} &> {log}"


rule map_new_insertions:
    input:
        is_queries="results/mutants/ismapper/best_IS_from_each_family.fasta",
        parent_ref="results/mutants/variants/{parent}/reference.fasta",
        mut_reads_1="results/data_filtered/{parent}/short/mutants/{parent}_1.fastq.gz",
        mut_reads_2="results/data_filtered/{parent}/short/mutants/{parent}_2.fastq.gz",
    output:
        directory("results/mutants/ismapper/new_insertions/{parent}"),
    threads: 10
    message:
        "Looking for new IS insertions in {wildcards.parent} mutant"
    log:
        "results/logs/{parent}_ISmapper.log",
    conda:
        "../envs/ismapper.yaml"
    container:
        config.get("default_container", None)
    shell:
        "ismap --queries {input.is_queries} --reads {input.mut_reads_1} {input.mut_reads_2} "
        "--reference {input.parent_ref} --t {threads} --output_dir {output} &> {log}"


rule final:
    input:
        bed="results/mutants/amplifications/{parent}/amplifications_windows.bed",
        genes_w_snps="results/mutants/variants/{parent}/genes_with_variants.tsv",
        amplifications="results/mutants/amplifications/{parent}/amplifications_annotated_filtered.tsv",
        rel_cov_mut="results/mutants/copy_number/{parent}/relative_coverage_mutant.tsv",
        rel_cov_parent="results/mutants/copy_number/{parent}/relative_coverage_parent.tsv",
        new_insert="results/mutants/ismapper/new_insertions/{parent}",
    output:
        touch("results/final/{parent}_mutants_all.done"),
    shell:
        "echo 'DONE'"


onsuccess:
    print("Mutants workflow finished, no errors")

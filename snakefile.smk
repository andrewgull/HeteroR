configfile: "config.yaml"

rule all:
    input:
        expand("busco_results/{strain}", strain=config['strains'])
# all quality_check_* rules produce all the expected files but job crashes with the following message (for ex):
# Job Missing files after 5 seconds:
# qualcheck/DA62942/Nanopore/DA62942_all.fastqc.zip
rule quality_check_illumina:
    input:
        "data_raw/{strain}/Illumina/{strain}_1.fq.gz"
    output:
        "qualcheck/{strain}/Illumina/{strain}_1.fastqc.zip"
    log: "logs/{strain}_illumina_qc.log"
    conda: "envs/rscripts.yaml"
    script:
        "scripts/run_qualcheck.R"

rule quality_check_trimmed:
    input:
        "data_filtered/{strain}/Illumina/{strain}_1.fq.gz"
    output:
        "qualcheck/{strain}/Illumina_trimmed/{strain}_1.fastqc.zip"
    log: "logs/{strain}_trimmed_qc.log"
    conda: "envs/rscripts.yaml"
    script:
        "scripts/run_qualcheck.R"

rule quality_check_nanopore:  # also builds length distribution plots
    input:
        "data_raw/{strain}/Nanopore/{strain}_all.fastq.gz"
    output:
        "qualcheck/{strain}/Nanopore/{strain}_all.fastqc.zip"
    log: "logs/{strain}_nanopore_qc.log"
    conda: "envs/rscripts.yaml"
    script:
        "scripts/run_qualcheck.R"

rule quality_check_filtered:
    input:
        "data_filtered/{strain}/Nanopore/{strain}_all.fastq.gz"
    output:
        "qualcheck/{strain}/Nanopore_filtered/{strain}_all.fastqc.zip"
    log: "logs/{strain}_filtered_qc.log"
    conda: "envs/rscripts.yaml"
    script:
        "scripts/run_qualcheck.R"

rule auto_trim_illumina:
    input:
        short_read_1 = "data_raw/{strain}/Illumina/{strain}_1.fq.gz",
        short_read_2 = "data_raw/{strain}/Illumina/{strain}_2.fq.gz"
    output:
        short_read_1 = "data_filtered/{strain}/Illumina/{strain}_1.fq.gz",
        short_read_2 = "data_filtered/{strain}/Illumina/{strain}_2.fq.gz",
        report_html = "qualcheck/fastp_reports/{strain}_report.html",
        report_json = "qualcheck/fastp_reports/{strain}_report.json"
    threads: 14
    message: "executing fastp with {threads} threads on {wildcards.strain} short reads"
    log: "logs/{strain}_fastp.log"
    conda: "envs/fastp.yaml"
    shell:
        "fastp --in1 {input.short_read_1} --in2 {input.short_read_2} --out1 {output.short_read_1} --out2 {output.short_read_2} "
        "--thread {threads} --qualified_quality_phred 20 --cut_window_size 4 "
        "--cut_right 20 --length_required 50 --html {output.report_html} --json {output.report_json}"

rule filtlong:
    input:
        "data_raw/{strain}/Nanopore/{strain}_all.fastq.gz"
    output:
        "data_filtered/{strain}/Nanopore/{strain}_all.fastq.gz"
    message:
        "executing filtlong on {wildcards.strain} long reads"
    log: "logs/{strain}_filtlong.log"
    conda: "envs/filtlong.yaml"
    threads: 14
    shell:
        "filtlong --min_length 1000 --keep_percent 90 {input} | pigz -c -p {threads} > {output}"

rule unicycler:
    input:
        short_read_1 = "data_filtered/{strain}/Illumina/{strain}_1.fq.gz",
        short_read_2 = "data_filtered/{strain}/Illumina/{strain}_2.fq.gz",
        long_read = "data_filtered/{strain}/Nanopore/{strain}_all.fastq.gz"
    output:
        directory("assemblies/{strain}")
    threads: 14
    message:
        "executing Unicycler with {threads} threads on {wildcards.strain} reads"
    log:
        "logs/{strain}_unicycler.log"
    conda: "envs/unicycler.yaml"
    shell:
        "unicycler -1 {input.short_read_1} -2 {input.short_read_2} -l {input.long_read} -t {threads} -o {output}"

rule quast:
    input:
        "assemblies/{strain}/assembly.fasta"
    output:
        directory("quast_results/{strain}")
    threads: 14
    message:
        "executing QUAST with {threads} threads on {wildcards.strain} assembly"
    log:
        "logs/{strain}_quast.log"
    conda: "envs/quast.yaml"
    shell:
        "quast.py -t {threads} -o {output} {input}"

rule busco:
    input:
        "assemblies/{strain}/assembly.fasta"
    output:
        directory("busco_results/{strain}")
    threads: 14
    message: "executing BUSCO with {threads} threads on {wildcards.strain} assembly"
    log: "logs/{strain}_busco.log"
    conda: "envs/busco.yaml"
    shell:
        "busco -m genome -i {input} -o {output} -l gammaproteobacteria_odb10 --cpu {threads}"

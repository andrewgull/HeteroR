configfile: "config/assembly.yaml"


strains = pd.read_csv(config["strains"], dtype={"strains": str})

reads_path = config["reads_path"]


rule trim_short:
    input:
        short_read_1=lambda wildcards: f"{reads_path}/{wildcards.strain}_1.fq.gz",
        short_read_2=lambda wildcards: f"{reads_path}/{wildcards.strain}_2.fq.gz",
    output:
        short_read_1="results/data_filtered/{strain}/short/{strain}_1.fq.gz",
        short_read_2="results/data_filtered/{strain}/short/{strain}_2.fq.gz",
        report_html="results/qualcheck_reads/fastp_reports/{strain}_report.html",
        report_json="results/qualcheck_reads/fastp_reports/{strain}_report.json",
    wildcard_constraints:
        strain="DA[0-9]*",
    threads: 18
    message:
        "executing fastp with {threads} threads on {wildcards.strain} short reads"
    log:
        "results/logs/{strain}_fastp.log",
    conda:
        "envs/fastp.yaml"
    container:
        config.get("default_container", "")
    params:
        q=config.get("quality", ""),
        W=config.get("window_size", ""),
        r=config.get("cut_right", ""),
        l=config.get("length_required", ""),
        f=config.get("trim_front", ""),
    shell:
        "fastp --in1 {input.short_read_1} --in2 {input.short_read_2} --out1 {output.short_read_1} "
        "--out2 {output.short_read_2} "
        "--thread {threads} --qualified_quality_phred {params.q} --cut_window_size {params.W} "
        "--cut_right {params.r} --length_required {params.l} --trim_front1 {params.f} --low_complexity_filter "
        "--html {output.report_html} --json {output.report_json} &> {log}"


# simple trimming of long reads
rule filter_long:
    input:
        lambda wildcards: f"{reads_path}/{wildcards.strain}.fastq.gz",
    output:
        "results/data_filtered/{strain}/long/{strain}.fastq.gz",
    message:
        "executing filtlong on {wildcards.strain} long reads"
    log:
        "results/logs/{strain}_filtlong.log",
    conda:
        "envs/filtlong.yaml"
    container:
        config.get("default_container", "")
    threads: 18
    params:
        min_len=config.get("min_nanopore_length", ""),
    shell:
        "filtlong --min_length {params.min_len} {input} 2> {log} | pigz -c -p {threads} > {output}"


# Make an assembly with Unicycler or FLye-Medaka-Polypolish depending on long coverage
rule adaptive_hybrid_assembly:
    input:
        short_reads_1="results/data_filtered/{strain}/short/{strain}_1.fq.gz",
        short_reads_2="results/data_filtered/{strain}/short/{strain}_2.fq.gz",
        long_reads="results/data_filtered/{strain}/long/{strain}.fastq.gz",
    output:
        assembly_dir=directory("results/assemblies/{strain}"),
        draft_dir=directory("results/drafts/{strain}"),
        polish_dir=directory("results/polished/{strain}"),
    threads: 18
    message:
        "executing assembly script with {threads} threads on {wildcards.strain} reads"
    log:
        "results/logs/{strain}_assembly.log",
    conda:
        "envs/hybrid_assembly.yaml"
    container:
        config.get("assembly_container", "")
    params:
        basecaller=config.get("basecaller", ""),
        genome_size=config.get("genome_size", ""),
        coverage=config.get("coverage", ""),
        genome_length=config.get("genome_length", ""),
        cov_threshold=config.get("cov_threshold", ""),
    script:
        "scripts/adaptive_hybrid_assembly.py"


# mapping of short reads on the assembly
rule map_back:
    input:
        assembly_dir="results/assemblies/{strain}",
        short_read_1="results/data_filtered/{strain}/short/{strain}_1.fq.gz",
        short_read_2="results/data_filtered/{strain}/short/{strain}_2.fq.gz",
    output:
        # an output file marked as temp is deleted after all rules that use it as an input are completed
        temp("results/mapping/{strain}/assembly.sam"),
    threads: 18
    message:
        "executing BWA with {threads} threads on {wildcards.strain} assembly"
    log:
        mem="results/logs/{strain}_bwa_mem.log",
        index="results/logs/{strain}_bwa_index.log",
    conda:
        "envs/bwa.yaml"
    container:
        config.get("default_container", "")
    shell:
        "bwa index {input.assembly_dir}/assembly.fasta &> {log.index} && "
        "bwa mem -t {threads} {input.assembly_dir}/assembly.fasta {input.short_read_1} {input.short_read_2} -o {output} &> {log.mem}"


# sorting and converting of the mapped short reads
rule sort_mapping:
    input:
        "results/mapping/{strain}/assembly.sam",
    output:
        "results/mapping/{strain}/assembly.bam",
    threads: 18
    message:
        "executing SAMTOOLS: VIEW-SORT with {threads} threads on {wildcards.strain} mapping file"
    log:
        "results/logs/{strain}_samtools.log",
    conda:
        "envs/samtools.yaml"
    container:
        config.get("default_container", "")
    shell:
        "samtools view -b {input} | samtools sort -o {output} -O BAM -@ {threads} &> {log} && samtools index {output}"


# collecting unmapped reads
rule collect_unmapped:
    input:
        "results/mapping/{strain}/assembly.bam",
    output:
        r1="results/mapping/{strain}/unmapped_1.fastq",
        r2="results/mapping/{strain}/unmapped_2.fastq",
    threads: 18
    message:
        "executing SAMTOOLS: VIEW-FASTQ with {threads} threads on {wildcards.strain} BAM file"
    log:
        "results/logs/{strain}_unmapped.log",
    conda:
        "envs/samtools.yaml"
    container:
        config.get("default_container", "")
    shell:
        "samtools view -@ {threads} -u -f 12 -F 256 {input} | samtools fastq -1 {output.r1} -2 {output.r2} -@ {threads} &> {log}"


# assembling plasmids from unmapped reads
rule additional_plasmid_assembly:
    input:
        r1="results/mapping/{strain}/unmapped_1.fastq",
        r2="results/mapping/{strain}/unmapped_2.fastq",
    output:
        directory("results/plasmids/{strain}"),
    threads: 18
    message:
        "executing SPAdes in plasmid mode with {threads} threads on unmapped reads of {wildcards.strain}"
    log:
        "results/logs/{strain}_spades.log",
    conda:
        "envs/spades.yaml"
    container:
        config.get("default_container", "")
    shell:
        # || true prevents the rule from failing when spades throws an error; this happens when unmapped files are too small
        "spades.py --plasmid -1 {input.r1} -2 {input.r2} -t {threads} -o {output} &> {log} || true"


# summarize assembly statistics based on Unicycler's log and SPAdes output
rule assembly_summary:
    input:
        "results/assemblies/{strain}",
        "results/plasmids/{strain}",
    output:
        "results/assemblies_joined/{strain}/summary.tsv",
    threads: 1
    params:
        position=config["position"],
    message:
        "summarizing unicycler and SPAdes assemblies of strain {wildcards.strain}"
    log:
        "results/logs/{strain}_assembly_summary.log",
    conda:
        "envs/biopython.yaml"
    container:
        config.get("biopython_container", "")
    script:
        "scripts/assembly_summary.py"


# merging Unicycler and SPAdes assemblies into single file
rule merge_assemblies:
    input:
        "results/assemblies/{strain}",
        "results/plasmids/{strain}",
    output:
        "results/assemblies_joined/{strain}/assembly.fasta",
    threads: 1
    message:
        "joining Unicycler assembly and SPAdes plasmid assembly together, strain {wildcards.strain}"
    log:
        "results/logs/{strain}_joiner.log",
    conda:
        "envs/biopython.yaml"
    container:
        config.get("biopython_container", "")
    script:
        "scripts/join_two_fastas.py"


# join outputs together
rule final:
    input:
        assembly_joined="results/assemblies_joined/{strain}/assembly.fasta",
        summary="results/assemblies_joined/{strain}/summary.tsv",
    output:
        touch("results/final/{strain}_assembly_all.done"),
    shell:
        "echo 'DONE'"


onsuccess:
    print("Assembly workflow finished, no errors")

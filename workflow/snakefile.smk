##############################################################################################
# the main pipeline for the Hetero-resistance project
# source: github.com/andrewgull/HeteroR
# author: Andrei Guliaev
# this pipeline assembles genomes, annotates resistance genes, direct repeats and IS elements
##############################################################################################

from snakemake.io import touch, directory, temp, expand


# Rule to join together all inputs and outputs
rule all:
    input:
        expand("results/final/{strain}_all.done", strain=config['strains'])

# Automated short read trimming
# fastp parameters:
        # -q, --qualified_quality_phred
        # the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
        # -W, --cut_window_size
        # the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4 (int [=4])
        # -r, --cut_right
        # move a sliding window from front to tail, if meet one window with mean quality < threshold,
        # drop the bases in the window and the right part, and then stop.
        # -l, --length_required
        # reads shorter than length_required will be discarded, default is 15. (int [=15])
        # -y, --low_complexity_filter
        # enable low complexity filter. The complexity is defined as the percentage of base that is
        # different from its next base (base[i] != base[i+1]).
        # -f, --trim_front1
        # trimming how many bases in front of read1, default is 0 (int [=0])
rule trim_short:
    input:
        short_read_1 = "resources/data_raw/{strain}/short/renamed/{strain}_1.fq.gz",
        short_read_2 = "resources/data_raw/{strain}/short/renamed/{strain}_2.fq.gz"
    output:
        short_read_1 = "results/data_filtered/{strain}/short/{strain}_1.fq.gz",
        short_read_2 = "results/data_filtered/{strain}/short/{strain}_2.fq.gz",
        report_html = "results/qualcheck_reads/fastp_reports/{strain}_report.html",
        report_json = "results/qualcheck_reads/fastp_reports/{strain}_report.json"
    wildcard_constraints: strain="DA[0-9]*"
    threads: 18
    message: "executing fastp with {threads} threads on {wildcards.strain} short reads"
    log: "results/logs/{strain}_fastp.log"
    conda: "envs/fastp.yaml"
    container: "containers/fastp.sif"
    params: q=config["quality"], W=config["window_size"], r=config["cut_right"], l=config["length_required"], f=config["trim_front"]
    shell:
        "fastp --in1 {input.short_read_1} --in2 {input.short_read_2} --out1 {output.short_read_1} "
        "--out2 {output.short_read_2} "
        "--thread {threads} --qualified_quality_phred {params.q} --cut_window_size {params.W} "
        "--cut_right {params.r} --length_required {params.l} --trim_front1 {params.f} --low_complexity_filter "
        "--html {output.report_html} --json {output.report_json} &> {log}"

# simple trimming of long reads
rule filter_long:
    input: "resources/data_raw/{strain}/long/{strain}_all.fastq.gz"
    output: "results/data_filtered/{strain}/long/{strain}_all.fastq.gz"
    message: "executing filtlong on {wildcards.strain} long reads"
    log: "results/logs/{strain}_filtlong.log"
    conda: "envs/filtlong.yaml"
    container: "containers/filtlong.sif"
    threads: 18
    params: min_len=config["min_long_length"]
    shell: "filtlong --min_length {params.min_len} {input} 2> {log} | pigz -c -p {threads} > {output}"

# Make an assembly with Unicycler or FLye-Medaka-Polypolish depending on long coverage
rule adaptive_hybrid_assembly:
    input:
        short_reads_1 = "results/data_filtered/{strain}/short/{strain}_1.fq.gz",
        short_reads_2 = "results/data_filtered/{strain}/short/{strain}_2.fq.gz",
        long_reads = "results/data_filtered/{strain}/long/{strain}_all.fastq.gz"
    output:
        assembly_dir = directory("results/assemblies/{strain}"),
        draft_dir = directory("results/drafts/{strain}"),
        polish_dir = directory("results/polished/{strain}")
    threads: 18
    message: "executing assembly script with {threads} threads on {wildcards.strain} reads"
    log: "results/logs/{strain}_assembly.log"
    conda: "envs/hybrid_assembly.yaml"
    container: "containers/hybrid_assembly.sif"
    params: basecaller=config["basecaller"], genome_size=config["genome_size"], coverage=config["coverage"], genome_length=config["genome_length"], cov_threshold=config["cov_threshold"]
    script: "scripts/adaptive_hybrid_assembly.py"

# assembly quality control
rule qc_assembly:
    input: "results/assemblies/{strain}", "resources/busco_downloads"
    output: directory("results/qualcheck_assembly/{strain}")
    threads: 18
    message: "executing BUSCO and QUAST with {threads} threads on {wildcards.strain} assembly"
    log: "results/logs/{strain}_assembly_qc.log"
    conda: "envs/busco_quast.yaml"
    container: "containers/assembly_qc.sif"
    params: tax_dataset=config["tax_dataset"]
    script: "scripts/QC_assembly.py"

# mapping of short reads on the assembly
rule map_back:
    input:
        assembly_dir = "results/assemblies/{strain}",
        short_read_1= "results/data_filtered/{strain}/short/{strain}_1.fq.gz",
        short_read_2 = "results/data_filtered/{strain}/short/{strain}_2.fq.gz"
    output:
        # an output file marked as temp is deleted after all rules that use it as an input are completed
        temp("results/mapping/{strain}/assembly.sam")
    threads: 18
    message: "executing BWA with {threads} threads on {wildcards.strain} assembly"
    log: mem = "results/logs/{strain}_bwa_mem.log",
         index = "results/logs/{strain}_bwa_index.log"
    conda: "envs/bwa.yaml"
    container: "containers/bwa.sif"
    shell:
        "bwa index {input.assembly_dir}/assembly.fasta &> {log.index} && "
        "bwa mem -t {threads} {input.assembly_dir}/assembly.fasta {input.short_read_1} {input.short_read_2} -o {output} &> {log.mem}"

# sorting and converting of the mapped short reads
rule sort_mapping:
    input:
        "results/mapping/{strain}/assembly.sam"
    output:
        "results/mapping/{strain}/assembly.bam"
    threads: 18
    message: "executing SAMTOOLS: VIEW-SORT with {threads} threads on {wildcards.strain} mapping file"
    log: "results/logs/{strain}_samtools.log"
    conda: "envs/samtools.yaml"
    container: "containers/samtools.sif"
    shell:
        "samtools view -b {input} | samtools sort -o {output} -O BAM -@ {threads} &> {log} && samtools index {output}"

# collecting unmapped reads
rule collect_unmapped:
    input:
        "results/mapping/{strain}/assembly.bam"
    output:
        r1 = "results/mapping/{strain}/unmapped_1.fastq",
        r2 = "results/mapping/{strain}/unmapped_2.fastq"
    threads: 18
    message: "executing SAMTOOLS: VIEW-FASTQ with {threads} threads on {wildcards.strain} BAM file"
    log: "results/logs/{strain}_unmapped.log"
    conda: "envs/samtools.yaml"
    container: "containers/samtools.sif"
    shell:
        "samtools view -@ {threads} -u -f 12 -F 256 {input} | samtools fastq -1 {output.r1} -2 {output.r2} -@ {threads} &> {log}"

# assembling plasmids from unmapped reads
rule additional_plasmid_assembly:
    input:
        r1 = "results/mapping/{strain}/unmapped_1.fastq",
        r2 = "results/mapping/{strain}/unmapped_2.fastq"
    output: directory("results/plasmids/{strain}")
    threads: 18
    message: "executing SPAdes in plasmid mode with {threads} threads on unmapped reads of {wildcards.strain}"
    log: "results/logs/{strain}_spades.log"
    conda: "envs/spades.yaml"
    container: "containers/spades.sif"
    shell:
        # || true prevents the rule from failing when spades throws an error; this happens when unmapped files are too small
        "spades.py --plasmid -1 {input.r1} -2 {input.r2} -t {threads} -o {output} &> {log} || true"

# summarize assembly statistics based on Unicycler's log and SPAdes output
rule assembly_summary:
    input:
        "results/assemblies/{strain}",
        "results/plasmids/{strain}"
    output: "results/assemblies_joined/{strain}/summary.tsv"
    threads: 1
    params: position=config["position"]
    message: "summarizing unicycler and SPAdes assemblies of strain {wildcards.strain}"
    log: "results/logs/{strain}_assembly_summary.log"
    conda: "envs/biopython.yaml"
    container: "containers/biopython.sif"
    script: "scripts/assembly_summary.py"

# merging Unicycler and SPAdes assemblies into single file
rule merge_assemblies:
    input:
        "results/assemblies/{strain}",
        "results/plasmids/{strain}"
    output:
        "results/assemblies_joined/{strain}/assembly.fasta"
    threads: 1
    message: "joining Unicycler assembly and SPAdes plasmid assembly together, strain {wildcards.strain}"
    log: "results/logs/{strain}_joiner.log"
    conda: "envs/biopython.yaml"
    container: "containers/biopython.sif"
    script: "scripts/join_two_fastas.py"

# Annotate merged assembly
rule assembly_annotation:
    # proteins (prodigal) and rRNA (barrnap)
    input:
        "results/assemblies_joined/{strain}/assembly.fasta"
    output:
        directory("results/annotations/{strain}/prokka")
    threads: 18
    message: "executing PROKKA with {threads} threads on full assembly of {wildcards.strain}"
    log: "results/logs/{strain}_prokka.log"
    conda: "envs/prokka.yaml"
    container: "containers/prokka.sif"
    params: centre=config["centre"], minlen=config["minlen"], genus=config["genus"], species=config["species"]
    shell:
        # skip tRNAs search?
        "prokka --addgenes --addmrna --compliant --notrna --outdir {output} --prefix {wildcards.strain}_genomic --centre {params.centre} --genus {params.genus} "
        "--species {params.species} --strain {wildcards.strain} --kingdom Bacteria --cpus {threads} "
        "--mincontiglen {params.minlen} {input} &> {log}"

# Change default names in GBK files to compatible with assembly's sequence IDs
rule rename_annotations:
    input: "results/annotations/{strain}/prokka"
    output: "results/annotations/{strain}/prokka_renamed/{strain}_genomic.gbk"
    params: filename="{strain}_genomic.gbk"
    message: "renaming {wildcards.strain} annotations"
    log: "results/logs/{strain}_rename_annot.log"
    conda: "envs/biopython.yaml"
    container: "containers/biopython.sif"
    script: "scripts/rename_genomic_gbk.py"

# Annotate tRNAs separately
rule trna_annotation:
    # tRNA genes only - tRNAScan-SE
    input:
        "results/assemblies_joined/{strain}/assembly.fasta"
    output:
        general = "results/annotations/{strain}/trna/trna_gen.txt",
        struct = "results/annotations/{strain}/trna/trna_struct.txt",
        iso = "results/annotations/{strain}/trna/trna_iso.txt",
        stats = "results/annotations/{strain}/trna/trna_stat.txt",
        bed = "results/annotations/{strain}/trna/trna_coords.bed",
        gff = "results/annotations/{strain}/trna/trna_feat.gff",
        fasta = "results/annotations/{strain}/trna/trna_seq.fasta"
    threads: 18
    message: "executing tRNAScan-SE with {threads} threads on full assembly of {wildcards.strain} strain"
    log: "results/logs/{strain}_trnascan.log"
    conda: "envs/trnascan.yaml"
    container: "containers/trnascan.sif"
    shell:
        "tRNAscan-SE -B --forceow -o {output.general} -f {output.struct} -s {output.iso} -m {output.stats} -b {output.bed} "
        "-j {output.gff} -a {output.fasta} -l {log} --thread {threads} {input} &> {log}"

# Merge both annotations
rule join_annotations:
    input:
        prokka="results/annotations/{strain}/prokka",
        trnascan="results/annotations/{strain}/trna/trna_seq.fasta"
    output: "results/annotations/{strain}/joined/annotation.fasta"  # it's not supposed to be used as input for RGI tool
    message: "joining {wildcards.strain} annotations"
    log: "results/logs/{strain}_join_annot.log"
    conda: "envs/biopython.yaml"
    container: "containers/biopython.sif"
    script: "scripts/join_two_fastas.py"

# Map short reads onto the joined assemblies to get coverage
rule genome_coverage:
    input:
        assembly="results/assemblies_joined/{strain}/assembly.fasta",
        fastq1="results/data_filtered/{strain}/short/{strain}_1.fq.gz",
        fastq2="results/data_filtered/{strain}/short/{strain}_2.fq.gz"
    output:
        bam=temp("results/genome_coverage/{strain}/assembly.bam"),
        depth="results/genome_coverage/{strain}/depth.txt"
    threads: 18
    message: "mapping reads onto {wildcards.strain} joined assembly"
    log: bwa_index="results/logs/{strain}_genome_coverage_bwa_index.log",
         bwa_mem="results/logs/{strain}_genome_coverage_bwa_mem.log",
         samtools_sort="results/logs/{strain}_genome_coverage_samtools_sort.log",
         samtools_index="results/logs/{strain}_genome_coverage_samtools_index.log",
         samtools_depth="results/logs/{strain}_genome_coverage_samtools_depth.log"
    conda: "envs/hybrid_assembly.yaml"
    container: "containers/hybrid_assembly.sif"
    shell: "bwa index {input.assembly} &> {log.bwa_index} && bwa mem -t {threads} {input.assembly} {input.fastq1} "
           "{input.fastq2} 2> {log.bwa_mem} | samtools sort -@ {threads} -o {output.bam} &> {log.samtools_sort} && "
           "samtools index {output.bam} &> {log.samtools_index} && "
           "samtools depth -a {output.bam} 2> {log.samtools_depth} > {output.depth}"

# Find resistance genes
# it ignores tRNAs
# to get RGI data base, run: rgi load --card_json ./card_database/card.json --local
rule resistance_genes:
    input:
        "results/annotations/{strain}/prokka"
    output:
        "results/resistance_genes/{strain}/rgi_table.txt" # IT'S JUST A PREFIX!
    threads: 18
    message: "executing RGI with {threads} threads on predicted genes/proteins from {wildcards.strain}"
    log: "results/logs/{strain}_rgi.log"
    conda: "envs/rgi.yaml"
    container: "containers/rgi.sif"
    shell:
        "output=$(echo '{output}' | cut -d'.' -f 1) && "
        "rgi main --input_sequence {input}/{wildcards.strain}_genomic.faa --output_file $output  "
        "--input_type protein --local  --num_threads {threads} --include_loose --clean &> {log}"

# Make separate file with annotations of the resistance genes
rule rg_annotation:
    input:
        "results/resistance_genes/{strain}/rgi_table.txt",
        "results/annotations/{strain}/prokka"
    output: "results/annotations/{strain}/resistance_genes/{strain}_resistance_genes.gbk"
    message: "annotating {wildcards.strain} resistance genes"
    log: "results/logs/{strain}_rgi_annot.log"
    conda: "envs/biopython.yaml"
    params: filter_criterion=config["filter_criterion"]
    script: "scripts/rgi2gff.py"

# make bed files with coords of regions around the resistance genes to find repeats in them
rule regions_coords:
    input:
        "results/assemblies_joined/{strain}/assembly.fasta",
        "results/annotations/{strain}/prokka",
        "results/resistance_genes/{strain}/rgi_table.txt"
    output:
       "results/direct_repeats/{strain}/regions/regions_within.bed",
       "results/direct_repeats/{strain}/regions/regions_overlapping_5_end.bed",
       "results/direct_repeats/{strain}/regions/regions_overlapping_3_end.bed"
    message: "creating BED files for RGs flanking regions in {wildcards.strain} assembly"
    log: "results/logs/{strain}_getbed.log"
    conda: "envs/biopython.yaml"
    container: "containers/biopython.sif"
    params: span=config["span"], min_plasmid_size=config["min_plasmid_size"]
    script: "scripts/flanking_regions.py"

# retrieve regions as fasta according to their coordinates
# bed tools returns an empty file if a bed file is empty
rule regions_seqs:
    input:
        assembly="results/assemblies_joined/{strain}/assembly.fasta",
        bed_normal="results/direct_repeats/{strain}/regions/regions_within.bed",
        bed_5_end="results/direct_repeats/{strain}/regions/regions_overlapping_5_end.bed",
        bed_3_end="results/direct_repeats/{strain}/regions/regions_overlapping_3_end.bed"
    output:
        normal="results/direct_repeats/{strain}/regions/regions_within.fasta",
        left="results/direct_repeats/{strain}/regions/regions_overlapping_5_end.fasta",
        right="results/direct_repeats/{strain}/regions/regions_overlapping_3_end.fasta"
    message: "retrieving regions' sequences from {wildcards.strain} assembly"
    log: "results/logs/{strain}_bedtools.log"
    conda: "envs/bedtools.yaml"
    container: "containers/bedtools.sif"
    shell:
        "samtools faidx {input.assembly} &&"
        "bedtools getfasta -fi {input.assembly} -bed {input.bed_normal} -nameOnly -fo {output.normal} &> {log}; "
        "bedtools getfasta -fi {input.assembly} -bed {input.bed_5_end} -nameOnly -fo {output.left} &>> {log}; "
        "bedtools getfasta -fi {input.assembly} -bed {input.bed_3_end} -nameOnly -fo {output.right} &>> {log}"

# join regions' ends overlapping 5' or 3'-ends
rule ends_overlaps:
    input:
        "results/direct_repeats/{strain}/regions/regions_within.fasta",
        "results/direct_repeats/{strain}/regions/regions_overlapping_5_end.fasta",
        "results/direct_repeats/{strain}/regions/regions_overlapping_3_end.fasta"
    output:
        "results/direct_repeats/{strain}/regions/regions_within_joined.fasta",
        "results/direct_repeats/{strain}/regions/regions_joined_5_end.fasta",
        "results/direct_repeats/{strain}/regions/regions_joined_3_end.fasta"
    message: "joining regions overlapping chromosome ends in {wildcards.strain} assembly"
    log: "results/logs/{strain}_join_ends.log"
    conda: "envs/biopython.yaml"
    container: "containers/biopython.sif"
    script: "scripts/join_ends.py"

# merge fasta files with the regions into a single file
rule merge_regions:
    input:
        normal="results/direct_repeats/{strain}/regions/regions_within_joined.fasta",
        left_joined="results/direct_repeats/{strain}/regions/regions_joined_5_end.fasta",
        right_joined="results/direct_repeats/{strain}/regions/regions_joined_3_end.fasta"
    output: "results/direct_repeats/{strain}/regions/regions_joined_final.fasta"
    message: "concatenating regions from {wildcards.strain} assembly"
    log: "results/logs/{strain}_concatenate_regions.txt"
    container: "containers/base.sif"
    shell: "cat {input.normal} {input.left_joined} {input.right_joined} > {output} 2> {log}"

# find direct repeat pairs in the regions
rule direct_repeats:
    input:
        "results/direct_repeats/{strain}/regions/regions_joined_final.fasta"
    output:
        directory("results/direct_repeats/{strain}/repeats_no_mismatch")
    threads: 18
    message: "executing GRF with {threads} threads on {wildcards.strain} assembly"
    log: "results/logs/{strain}_grf_perfect.log"
    conda: "envs/grf.yaml"
    container: "containers/grf.sif"
    params: mode=config["mode"], min_size=config["min_size"], format=config["format"], mism=config["mism"], seed_mism=config["seed_mism"], max_dist=config["max_dist"], min_dist=config["min_dist"]
    shell:
        "grf-main -i {input} -c {params.mode} -o {output} -t {threads} --min_tr {params.min_size} -f {params.format} "
        "--max_mismatch {params.mism} --seed_mismatch {params.seed_mism} --max_space {params.max_dist} --min_space {params.min_dist} &> {log} "

# make gff files for repeats found
# coordinates are LOCAL, i.e. related to a region around RG not the whole chromosome
rule dr_annotation:
    input:
        "results/direct_repeats/{strain}/repeats_no_mismatch",
        "results/annotations/{strain}/prokka"
    output:
        "results/annotations/{strain}/repeats/{strain}_repeats_no_mismatch_perfect.gff",
        "results/annotations/{strain}/repeats/{strain}_repeats_no_mismatch_imperfect.gff"
    message: "executing GFF_parser.py on {wildcards.strain} perfect repeats data"
    log: "results/logs/{strain}_gff_perfect.log"
    conda: "envs/biopython.yaml"
    container: "containers/biopython.sif"
    params: min_len=config["min_repeat_length"]
    script: "scripts/GRF_parser.py"

# make a table with repeats coordinates, remove duplicates
rule dr_table:
    input:
        "results/direct_repeats/{strain}/repeats_no_mismatch"
    output:
        "results/annotations/{strain}/repeats/{strain}_repeats.csv"
    message: "making CSV files for repeat pairs in strain {wildcards.strain}"
    log: "results/logs/{strain}_repeats_csv.log"
    conda: "envs/biopython.yaml"
    container: "containers/biopython.sif"
    script: "scripts/make_repeat_tables.py"

# join repeat tables, label repeat pairs as spanning RG centers or not, calculate AR length
rule dr_summary:
    input:
        expand("results/annotations/{strain}/repeats/{strain}_repeats.csv", strain=config["strains"]),
        expand("results/direct_repeats/{strain}/regions/regions_within.bed", strain=config["strains"])
    output: "results/tables/repeats_summary.csv"
    message: "making summary table with repeat coordinates for all strains"
    log: "results/logs/repeat_summary.log"
    conda: "envs/biopython.yaml"
    container: "containers/biopython.sif"
    script: "scripts/make_repeat_summary_table.py"

rule isescan:
    input: "results/direct_repeats/{strain}/regions/regions_joined_final.fasta"
    output: directory("results/isescan/{strain}")
    threads: 18
    message: "executing ISEScan on {wildcards.strain}"
    conda: "envs/isescan.yaml"
    container: "containers/isescan.sif"
    log: "results/logs/{strain}_isescan.log"
    shell: "isescan.py --seqfile {input} --output {output} --nthread {threads} &> {log}"

# join outputs together
rule final:
    input:
        qc_ass="results/qualcheck_assembly/{strain}",
        trnascan="results/annotations/{strain}/trna/trna_gen.txt",
        assembly="results/assemblies/{strain}",
        assembly_joined="results/assemblies_joined/{strain}/assembly.fasta",
        summary="results/assemblies_joined/{strain}/summary.tsv",
        gff_nomism="results/annotations/{strain}/repeats/{strain}_repeats_no_mismatch_perfect.gff",
        direct_repeats="results/direct_repeats/{strain}/regions/regions_joined_final.fasta",
        rg_gbk="results/annotations/{strain}/resistance_genes/{strain}_resistance_genes.gbk",
        renamed_gbk="results/annotations/{strain}/prokka_renamed/{strain}_genomic.gbk",
        rep_csv="results/annotations/{strain}/repeats/{strain}_repeats.csv",
        rep_sum="results/tables/repeats_summary.csv",
        depth="results/genome_coverage/{strain}/depth.txt",
        isescan="results/isescan/{strain}"
    output: touch("results/final/{strain}_all.done")
    shell: "echo 'DONE'"

onsuccess:
    print("Workflow finished, no errors")

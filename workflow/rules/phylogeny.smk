configfile: "config/phylogeny.yaml"


rule annotate_genes:
    # proteins (prodigal) and rRNA (barrnap)
    input:
        "results/assemblies_joined/{strain}/assembly.fasta",
    output:
        directory("results/annotations/{strain}/prokka"),
    threads: 10
    message:
        "executing PROKKA with {threads} threads on full assembly of {wildcards.strain}"
    log:
        "results/logs/{strain}_prokka.log",
    conda:
        "../envs/prokka.yaml"
    container:
        config.get("prokka_container", None)
    params:
        centre=config.get("centre", "centre_name"),
        minlen=config.get("minlen", 200),
        genus=config.get("genus", "genus_name"),
        species=config.get("species", "species_name"),
    shell:
        "prokka --addgenes --addmrna --compliant --notrna --outdir {output} --prefix {wildcards.strain}_genomic --centre {params.centre} --genus {params.genus} "
        "--species {params.species} --strain {wildcards.strain} --kingdom Bacteria --cpus {threads} "
        "--mincontiglen {params.minlen} {input} &> {log}"


rule core_genome:
    input:
        expand("results/annotations/{strain}/prokka/{strain}_genomic.gff", strain=get_strains(config)),
    output:
        directory("results/phylogeny/core_genome/"),
    threads: 10
    message:
        "executing basic roary"
    log:
        "results/logs/roary.log",
    conda:
        "../envs/roary.yaml"
    container:
        config.get("roary_container", None)
    shell:
        "roary -f {output} -p {threads} -e -n -v {input} &> {log}"


rule phylogenetic_tree:
    input:
        alignment="results/phylogeny/core_genome",
    output:
        directory("results/phylogeny/tree"),
    threads: 10
    message:
        "executing IQTree on core gene alignment"
    log:
        "results/logs/iqtree.log",
    conda:
        "../envs/iqtree.yaml"
    container:
        config.get("iqtree_container", None)
    params:
        bootstrap=config.get("bootstrap", 1000),
        alrt=config.get("alrt", 1000),
        prefix=config.get("prefix", "core_tree"),
    shell:
        # no partitions here, otherwise it will take forever
        "iqtree -s {input.alignment}/core_gene_alignment.aln --seqtype DNA -m GTR+I+G "
        "-T AUTO -ntmax {threads} --prefix {params.prefix} -B {params.bootstrap} -alrt {params.alrt} -bnni --safe &> {log}; "
        "[ -d {output} ] || mkdir {output}; mv {params.prefix}.* {output}"


rule plot_tree:
    input:
        tree_dir="results/phylogeny/tree",
        labels="notebooks/modelling/data/heteroresistance_testing.csv",
    output:
        "results/phylogeny/tree/core_genome_tree.pdf",
    message:
        "plotting the core gneome tree"
    log:
        "results/log/plot_tree.log",
    conda:
        "../envs/plottreer.yaml"
    container:
        config.get("rscripts_container", None)
    params:
        filename=config.get("filename", "core_tree"),
        width=config.get("width", 10),
        height=config.get("height", 10),
        units=config.get("units", "in"),
        outgroup=config.get("outgroup", "outgroup"),
    script:
        "../scripts/R/plot_tree.R"


rule final:
    input:
        annot="results/annotations/{strain}/prokka/{strain}_genomic.gff",
        pan="results/phylogeny/core_genome",
        plot="results/phylogeny/tree/core_genome_tree.pdf",
    output:
        touch("results/final/{strain}_phylogeny_all.done"),
    shell:
        "echo 'DONE'"

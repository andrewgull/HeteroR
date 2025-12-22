configfile: "config/phylogeny.yaml"


strains = pd.read_csv(config["strains"], dtype={"strains": str})


rule annotate_genes:
    # proteins (prodigal) and rRNA (barrnap)
    input:
        "results/assemblies_joined/{strain}/assembly.fasta",
    output:
        directory("results/annotations/{strain}/prokka"),
    threads: 18
    message:
        "executing PROKKA with {threads} threads on full assembly of {wildcards.strain}"
    log:
        "results/logs/{strain}_prokka.log",
    conda:
        "envs/prokka.yaml"
    params:
        centre=config["centre"],
        minlen=config["minlen"],
        genus=config["genus"],
        species=config["species"],
    shell:
        "prokka --addgenes --addmrna --compliant --notrna --outdir {output} --prefix {wildcards.strain}_genomic --centre {params.centre} --genus {params.genus} "
        "--species {params.species} --strain {wildcards.strain} --kingdom Bacteria --cpus {threads} "
        "--mincontiglen {params.minlen} {input} &> {log}"


rule core_genome:
    input:
        expand("results/annotations/{strain}/prokka/{strain}_genomic.gff", strain=strains["strains"]),
    output:
        directory("results/phylogeny/core_genome/"),
    threads: 14
    message:
        "executing basic roary"
    log:
        "results/logs/roary.log",
    conda:
        "envs/roary.yaml"
    shell:
        # check options
        "roary -f {output} -p {threads} -e -n -v {input} &> {log}"


rule phylogenetic_tree:
    input:
        alignment="results/phylogeny/core_genome",
    output:
        directory("results/phylogeny/tree"),
    threads: 14
    message:
        "executing IQTree on core gene alignment"
    log:
        "results/logs/iqtree.log",
    conda:
        "envs/iqtree.yaml"
    params:
        bootstrap=config["bootstrap"],
        alrt=config["alrt"],
        prefix=config["prefix"],
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
        "envs/plottreer.yaml"
    params:
        filename=config["filename"],
        width=config["width"],
        height=["height"],
        units=config["units"],
        outgroup=config["outgroup"],
    script:
        "scripts/plot_tree.R"


rule final:
    input:
        annot="results/annotations/{strain}/prokka/{strain}_genomic.gff",
        pan="results/phylogeny/core_genome",
        plot="results/phylogeny/tree/core_genome_tree.pdf",
    output:
        touch("results/final/{strain}_phylogeny_all.done"),
    shell:
        "echo 'DONE'"

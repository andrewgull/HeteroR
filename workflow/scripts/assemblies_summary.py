import subprocess
import pandas as pd
import os
from Bio import SeqIO


strains = list()
col_names = ["Component", "Segments", "Links", "Length", "N50", "Longest_component", "Status"]
for strain in strains:
    summary = subprocess.run("sed -n '/^Component/,/^Polishing/{p;/^Polishing/q}' "
                             "assemblies/%s/unicycler.log | head -n -3" % strain, shell=True, capture_output=True)
    summary_stdout = summary.stdout.decode("utf-8").splitlines()
    # i can skip first two elements (header and total): summary_lines[2:len(summary_lines)]
    summary_lines = [item.split() for item in summary_stdout][2:len(summary_stdout)]
    summary_df = pd.DataFrame(summary_lines, columns=col_names)
    summary_df["Strain"] = strain
    summary_df['Type'] = ['Chromosome'] + ['Plasmid'] * (len(summary_df) - 1)

    # get stats from plasmid assemblies
    # finished plasmids assembly is called scaffolds.fasta
    plasmid_assembly = "plasmids/%s/scaffolds.fasta" % strain
    if os.path.isfile(plasmid_assembly) and os.path.getsize(plasmid_assembly) > 0:
        plasmid = [len(seq) for seq in SeqIO.parse(plasmid_assembly, 'fasta')]



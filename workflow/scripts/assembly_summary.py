import subprocess
import pandas as pd
import os
from Bio import SeqIO


input_assembly = snakemake.input[0]  # it's a dir path like 'results/assemblies/DA00000'
input_plasmid = snakemake.input[1]  # dir name like 'results/plasmids/DA00000'
output = snakemake.output[0]

# get strain name
strain = input_assembly.split('/')[2]
# names for future table
col_names = ["Component", "Segments", "Links", "Length", "N50", "Longest_component", "Status"]

# parse unicycler's log file
summary = subprocess.run("sed -n '/^Component/,/^Polishing/{p;/^Polishing/q}' "
                         "%s/unicycler.log | head -n -3" % input_assembly, shell=True, capture_output=True)
summary_stdout = summary.stdout.decode("utf-8").splitlines()
# i can skip first two elements (header and total): summary_lines[2:len(summary_lines)]
if len(summary_stdout) > 2:
    # there is 'total' which should be removed ( as well as the original header)
    summary_lines = [item.split() for item in summary_stdout][2:len(summary_stdout)]
else:
    # there is no 'total', just one header and one chromosome
    # remove only original header
    summary_lines = [item.split() for item in summary_stdout[1:]]

summary_df = pd.DataFrame(summary_lines, columns=col_names)
summary_df["Strain"] = strain
# ideally there should be one chromosomal sequence and the rest should be plasmids
summary_df['Type'] = ['Chromosome'] + ['Plasmid'] * (len(summary_df) - 1)

# get stats from plasmid assemblies: count and length
# finished plasmids assembly is called 'scaffolds.fasta'
# create a DataFrame
plasmid_assembly = "%s/scaffolds.fasta" % input_plasmid
if os.path.isfile(plasmid_assembly) and os.path.getsize(plasmid_assembly) > 0:
    plasmid_len = [len(seq) for seq in SeqIO.parse(plasmid_assembly, 'fasta')]
    plasmid_df = pd.DataFrame({'Length': plasmid_len})
    plasmid_df['Component'] = None
    plasmid_df['Segments'] = 1
    plasmid_df['Links'] = None
    plasmid_df['N50'] = None
    plasmid_df['Longest_component'] = plasmid_len
    plasmid_df['Status'] = 'scaffold'
    plasmid_df['Strain'] = strain
    plasmid_df['Type'] = 'Plasmid'
    # concatenate summary_df and plasmid_df
    joined_df = pd.concat([summary_df, plasmid_df])
    # rewrite 'Component' to include new plasmids
    joined_df['Component'] = [i for i in range(1, len(joined_df) + 1)]
else:
    joined_df = summary_df
# write this DataFrame to a file
# output is declared on the top of he script (and comes from snakemake)
joined_df.to_csv(path_or_buf=output, sep="\t", index=False)

import subprocess
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
import sys


input_assembly = snakemake.input[0]  # it's a dir path like 'results/assemblies/DA00000'
input_plasmid = snakemake.input[1]  # dir name like 'results/plasmids/DA00000'
strain_position = snakemake.params[0]  # index of strain name position in input assembly, 2 for normal assembling, 3 for reassembling
output = snakemake.output[0]

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    # GET STRAIN NAME
    strain = input_assembly.split('/')[strain_position]

    # READ POLISHED ASSEMBLY INFO OR PARSE UNICYCLER LOG
    polished_info = "%s/assembly_info.txt" % input_assembly
    if os.path.isfile(polished_info):
        summary_df = pd.read_csv(polished_info, delimiter="\t")
        # leave seq_name as component, length as length, complete as status, add strain and type
        summary_df.drop(["alt_group", "graph_path"], axis=1, inplace=True)
        summary_df.columns = ["Component", "Length", "Coverage", "Status", "Repeat", "Mult"]
        summary_df["Type"] = ["Chromosome"] + ["Plasmid"]*(len(summary_df)-1)
        summary_df["Strain"] = strain
        summary_df["Segments"] = np.nan
        summary_df["Links"] = np.nan
        summary_df["N50"] = np.nan
        summary_df["Longest_component"] = np.nan
        # replace N/C in Status with incomplete/complete
        summary_df["Status"] = summary_df['Status'].str.replace("N", "incomplete")
        summary_df["Status"] = summary_df['Status'].str.replace("C", "complete")
    else:
        # PARSE UNICYCLER'S LOG
        col_names = ["Component", "Segments", "Links", "Length", "N50", "Longest_component", "Status"]
        # parse unicycler's log file
        summary = subprocess.run("sed -n '/^Component/,/^Polishing/{p;/^Polishing/q}' "
                                "%s/unicycler.log | head -n -3" % input_assembly,
                                shell=True, capture_output=True, text=True)
        summary_stdout = summary.stdout.splitlines()
        # I can skip first two elements (header and total): summary_lines[2:len(summary_lines)]
        if len(summary_stdout) > 2:
            # there is 'total' which should be removed ( as well as the original header)
            summary_lines = [item.split() for item in summary_stdout][2:len(summary_stdout)]
        else:
            # there is no 'total', just one header and one chromosome
            # remove only original header
            summary_lines = [item.split() for item in summary_stdout[1:]]

        summary_df = pd.DataFrame(summary_lines, columns=col_names)
        summary_df["Strain"] = strain
        summary_df["Repeat"] = np.nan
        summary_df["Mult"] = np.nan
        summary_df["Coverage"] = np.nan
        # ideally there should be one chromosomal sequence and the rest should be plasmids
        summary_df['Type'] = ['Chromosome'] + ['Plasmid'] * (len(summary_df) - 1)

    # summary_df always has the same set of columns: when it comes from unicycler
    # or when it comes from flye-medaka-polyploish

    # GET STATS FROM PLASMID ASSEMBLY
    # finished plasmids assembly is called 'scaffolds.fasta'
    # create a DataFrame with the same set of columns
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
        plasmid_df["Repeat"] = np.nan
        plasmid_df["Mult"] = np.nan
        plasmid_df["Coverage"] = np.nan
        # concatenate summary_df and plasmid_df
        joined_df = pd.concat([summary_df, plasmid_df])
        # rewrite 'Component' to include new plasmids
        joined_df['Component'] = [i for i in range(1, len(joined_df) + 1)]
    else:
        joined_df = summary_df
    # write this DataFrame to a file
    # output is declared on the top of the script (and comes from snakemake)
    joined_df.to_csv(path_or_buf=output, sep="\t", index=False)

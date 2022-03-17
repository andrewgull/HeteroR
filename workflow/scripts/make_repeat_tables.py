# script to make a csv
import sys
sys.path.append("/home/andrei/GitProjects/HeteroR/workflow/scripts")
from GRF_parser import parse_grf_output_no_ranges
import pandas as pd


def make_repeats_df(spacer_file, strain_index=0):
    # remove duplicated rows!
    with open(spacer_file) as f:
        output_lines = [line.rstrip() for line in f.readlines()]
    parsed_lines = [parse_grf_output_no_ranges(line) for line in output_lines]
    spacer_df = pd.DataFrame(columns=["record_id", "start_1", "end_1", "start_2", "end_2", "length"], data=parsed_lines)
    spacer_df["strain"] = spacer_file.split("/")[strain_index]
    return spacer_df


# snakemake.input[0] looks like "results/direct_repeats/{strain}/repeats_no_mismatch/perfect.spacer.id"
repeat_df = make_repeats_df(snakemake.input[0], strain_index=2)
repeat_df.drop_duplicates(inplace=True)
repeat_df.to_csv(snakemake.output[0])

# script to make a csv
import sys
sys.path.append("/home/andrei/GitProjects/HeteroR/workflow/scripts")
from GRF_parser import parse_grf_output_no_ranges
import pandas as pd
import os


def make_repeats_df(spacer_file, strain_index=0):
    # remove duplicated rows!
    with open(spacer_file) as f:
        output_lines = [line.rstrip() for line in f.readlines()]
    parsed_lines = [parse_grf_output_no_ranges(line) for line in output_lines]
    spacer_df = pd.DataFrame(columns=["record_id", "start_1", "end_1", "start_2", "end_2", "length"], data=parsed_lines)
    spacer_df["strain"] = spacer_file.split("/")[strain_index]
    return spacer_df


if __name__ == '__main__':
    with open(snakemake.log[0], "w") as f:
        sys.stderr = sys.stdout = f
        in_file = os.path.join(snakemake.input[0], "perfect.spacer.id")
        repeat_df = make_repeats_df(in_file, strain_index=2)
        repeat_df.drop_duplicates(inplace=True)
        repeat_df.to_csv(snakemake.output[0])

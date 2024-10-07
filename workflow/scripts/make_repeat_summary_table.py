#!/usr/bin/env python3
# script to join repeat summary tables and calculate AR length and span center
import pandas as pd
import numpy as np
import sys

# open log
with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    # inputs
    last_file = int(len(snakemake.input)/2)  # index of the first file in the bed input
    repeat_files_lst = [snakemake.input[i] for i in range(last_file)]  # doesn't include last_file
    bed_files_lst = [snakemake.input[i] for i in range(last_file, len(snakemake.input))]  # does include last file

    # read csv files with repeat pairs and join them into a DataFrame
    repeat_df_list = [pd.read_csv(f) for f in repeat_files_lst]
    repeat_df = pd.concat(repeat_df_list)

    # calculate Amplifiable Region length
    repeat_df["AR_length"] = repeat_df.start_2 - repeat_df.end_1 + 1

    # read bed files and join them into a DataFrame
    bed_df_list = [pd.read_csv(f, delimiter="\t", header=None) for f in bed_files_lst]
    bed_df = pd.concat(bed_df_list)

    # rename V4 (4th col) to 'record_id'
    bed_df = bed_df.rename(columns={3: "record_id"})
    # calculate gene center location
    bed_df["gene_center"] = (bed_df[2] - bed_df[1] + 1)/2
    # select only required columns
    bed_df = bed_df[["record_id", "gene_center"]]

    # join repeats and bed data frames
    # ow=outer - use union of keys from both data frames, analogous to inner join that I used in dplyr
    repeat_df_merged = pd.merge(repeat_df, bed_df, on="record_id", how="outer")

    # make spans_center column
    repeat_df_merged["spans_center"] = np.where((repeat_df_merged['end_1'] <= repeat_df_merged["gene_center"]) &
                                                (repeat_df_merged["start_2"] >= repeat_df_merged["gene_center"]), 'yes', 'no')
    # drop the first column
    repeat_df_merged.drop("Unnamed: 0", axis=1, inplace=True)

    # write to output
    repeat_df_merged.to_csv(snakemake.output[0])

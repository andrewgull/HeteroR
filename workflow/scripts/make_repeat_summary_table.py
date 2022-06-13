# script to join repeat summary tables and calculate AR length and span center
import pandas as pd
import glob

repeat_files_lst = ["/home/andrei/Data/HeteroR/results/annotations/DA62886/repeats/DA62886_repeats.csv",
                    "/home/andrei/Data/HeteroR/results/annotations/DA62888/repeats/DA62888_repeats.csv"]

repeat_df_list = [pd.read_csv(f) for f in repeat_files_lst]
repeat_df = pd.concat(repeat_df_list)

repeat_df["AR_length"] = repeat_df.start_2 - repeat_df.end_1 + 1

bed_files_lst = ["/home/andrei/Data/HeteroR/results/direct_repeats/DA62886/regions/regions_within.bed",
                 "/home/andrei/Data/HeteroR/results/direct_repeats/DA62888/regions/regions_within.bed"]

bed_df_list = [pd.read_csv(f, delimiter="\t", header=None) for f in bed_files_lst]
bed_df = pd.concat(bed_df_list)

bed_df.rename(columns={3: "record_id"})

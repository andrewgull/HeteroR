# script for making the feature table - should it be SQLite?
# list of features:
# 1. number of RG
# 2. number of RG on plasmids
# 3. nearest distance to oriC
# 4. median distance to oriC
# 5. total number of repeats
# 6. median number of repeats
# 7. median repeat length
# 8. longest repeat length
# 9. match concetration (=number of mismatches)
# 10. amplifiable region (AR) min length
# 11. AR median length
# 12. BLACK BOX

import pandas as pd


# strains list should come from the config
strains = ["DA63084", "DA63186", "DA63322", "DA63946", "DA64026"]

# paths to tables with repeats and RGs
path_to_repeats_csv = "results/annotations/{strain}/repeats/%s_repeats.csv"
path_to_rgi = "/home/andrei/Data/HeteroR/results/resistance_genes/%s/rgi_table.txt"

# make full tables
repeats_df_lst = [pd.read_csv(path_to_repeats_csv % strain, delimiter="\t") for strain in strains]
rgi_df_list = [pd.read_csv(path_to_rgi % strain, delimiter="\t") for strain in strains]

repeat_df = pd.concat(repeats_df_lst)
rgi_df = pd.concat(rgi_df_list)
rgi_notLoose = rgi_df[rgi_df["Cut_Off"] != "Loose"]

repeat_df.head()

# get list of all available ABs
antibiotics = list(set(rgi_notLoose["Drug Class"]))

# filter macrolides resistant
rgi_macrolides = rgi_notLoose[rgi_notLoose["Drug Class"].str.contains("macrolide")]
rgi_fqn = rgi_notLoose[rgi_notLoose["Drug Class"].str.contains("fluoroquinolone")]


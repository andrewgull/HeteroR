from Bio import SeqIO
from glob import glob
import pandas as pd


# cd '/home/andrei/Data/HeteroR/'

# for each strain
# go through annotation gbk files,
# collect genes only from plasmids
# make a df
# collect these genes from RG tables (left-join)
# write results in one table with one column for strain name

annotation_files = glob("results/annotations/DA*/resistance_genes/DA*_resistance_genes.gbk")
total_df_list = list()
for f in annotation_files:
    strain = f.split("/")[2]
    record = [rec for rec in SeqIO.parse(f, "genbank") if rec.id != "1"]
    if len(record) > 0:
        plasmid_dict = dict()
        final_df_list = []
        for rec in record:
            loci = [feature.qualifiers['locus_tag'][0] for feature in rec.features]
            plasmid_id = [rec.id] * len(loci)
            plasmid_dict["gene_id"] = loci
            plasmid_dict["plasmid_id"] = plasmid_id
            plasmid_df = pd.DataFrame.from_dict(plasmid_dict)
            final_df_list.append(plasmid_df)

        all_plasmid_df = pd.concat(final_df_list)
        all_plasmid_df["strain"] = strain
        # read the corresponding rgi file
        rg_df = pd.read_csv("results/resistance_genes/%s/rgi_table.txt" % strain, delimiter="\t")
        # transform the first column cause it's stupid
        rg_df["gene_id"] = list(map(lambda x: x.split(" ")[0], rg_df["ORF_ID"]))
        # join rg_df and all_plasmid_df
        merged_df = pd.merge(all_plasmid_df, rg_df, on='gene_id', how='left')
        total_df_list.append(merged_df)

output_df = pd.concat(total_df_list)
output_df.to_csv("plasmid_RGs.csv", index=False, sep="\t")

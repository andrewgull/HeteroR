# script to filter out repeat pairs that do not span a region's center
library(dplyr)

filter_center <- function(bed_file, repeat_table){
  # function to filter out repeat pairs not spanning a region's center
  # bed_file: results/direct_repeats/{strain}/regions/regions_within.bed
  # repeat_table: results/annotations/{strain}/repeats/{strain}_repeats.csv
  bed_df <- read.delim(bed_file, header=FALSE)
  bed_df$gene_center <- round((bed_df$V3 - bed_df$V2 + 1)/2)
  bed_df <- rename(bed_df, "record_id"=V4)
  bed_df <- select(bed_df, record_id, gene_center)
  # read {strain}_repeats.csv
  repeat_df <- read.csv(repeat_table, header=TRUE)
  repeat_df <- full_join(bed_df, repeat_df, by="record_id")
  repeat_df$spans_center <- if_else(repeat_df$end_1 <= repeat_df$gene_center & repeat_df$start_2 >= repeat_df$gene_center, "yes", "no")
  repeat_df_center <- filter(repeat_df, spans_center == "yes")
  return(repeat_df_center)
}

centered_repeats_df <- filter_center(snakemake@input[[1]], snakemake@input[[2]])
write.csv(x=centered_repeats_df, file=snakemake@output[[1]], row.names = FALSE)

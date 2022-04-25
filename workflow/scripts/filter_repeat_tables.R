# filter out repeat pairs that do not span a region's center

# libraries
library(dplyr)
library(optparse)

# CLI parsing
option_list = list(
   make_option(c("-b", "--bed"),
               type = "character",
               default = NULL,
               help = "A path to a bed file",
               metavar = "character"),
   make_option(c("-r", "--table"),
               type = "character",
               default = NULL,
               help = "A path to a table with repeat pairs (csv)",
               metavar = "character"),
	make_option(c("-o", "--out"),
                type = "character",
                default = NULL,
                help = "output file name",
                metavar = "character")
);

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if (is.null(opt$bed)){
 print_help(opt_parser)
 stop("A filepath to a bed file must be provided", call. = FALSE)
}
if (is.null(opt$table)){
 print_help(opt_parser)
 stop("A filepath to a repeat table must be provided", call. = FALSE)
}

# get file paths
bed_file_path <- opt$bed
repeat_table_path <- opt$table
output_file <- opt$out

# main function
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

# old way to execute the function via snakemake object
# centered_repeats_df <- filter_center(snakemake@input[[1]], snakemake@input[[2]])
# write.csv(x=centered_repeats_df, file=snakemake@output[[1]], row.names = FALSE)

# newer way - via Rscript (from @lanchlandeer github)
centered_repeats_df <- filter_center(bed_file_path, repeat_table_path)
write.csv(x=centered_repeats_df, file=output_file, row.names = FALSE)
# script to read depth files and calculate relative coverage
library(optparse)

# CLI parsing
option_list <- list(
    make_option(c("-i", "--input"),
                type = "character",
                default = NULL,
                help = "depth file",
                metavar = "path"),
    make_option(c("-s", "--strain"),
                type = "character",
                default = NULL,
                help = "strain name",
                metavar = "char"),
    make_option(c("-l", "--label"),
                type = "character",
                default = NULL,
                help = "mutant or parent",
                metavar = "char"),
    make_option(c("-m", "--min_contig_len"),
                type = "integer",
                default = 1000,
                help = "length in nt",
                metavar = "int"),
    make_option(c("-o", "--output"),
                type = "character",
                default = NULL,
                help = "output file (tsv)",
                metavar = "path")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
 print_help(opt_parser)
 stop("Input file must be provided", call. = FALSE)
}
if (is.null(opt$output)){
 print_help(opt_parser)
 stop("Output file must be provided", call. = FALSE)
}

#### Libraries ####
suppressPackageStartupMessages(library(dplyr))
library(readr)

#### Functions ####
parse_depth <- function(file_path, min_contig_length = 1000, strain, label) {
  # param file_path: depth file path
  # min_contig_length: contigs to filter out
  # param strain: strain name
  # param label: mutant or parent
  depth_df <- read_tsv(file_path, col_names = F, show_col_types = F)
  
  # check correct column names
  stopifnot(exprs = {"X1" %in% names(depth_df) & "X3" %in% names(depth_df)})
  # check if the table is empty
  stopifnot(nrow(depth_df) > 0)
  
  # calculate length and avg. coverage of each contig
  cov_df <-
    depth_df %>%
    group_by(X1) %>%
    summarise(length = n(),
              coverage = sum(X3)) %>%
    filter(length > min_contig_length) %>%
    mutate(avg.cov = round(coverage / length, 3))
  
  # relative coverage
  chrom_cov <- cov_df %>%
    filter(length == max(length)) %>%
    pull(avg.cov)
  
  cov_df$rel.cov <- cov_df$avg.cov / chrom_cov
  cov_df$strain <- strain
  cov_df$label <- label
  
  return(cov_df)
}

#### Main ####
coverage_table <- parse_depth(opt$input, opt$min_contig_len, opt$strain, opt$label) # nolint: line_length_linter.
write_delim(coverage_table, file = opt$output, delim = "\t")
print("Finished, no errors.")

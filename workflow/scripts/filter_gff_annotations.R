# script for reading GFF file and filtering based on the 10th column + retrieving gene names
# this column contains the fraction of GFF file covered by VCF file
# if values in this column greater than 0, then the corresponding gene contains a variant
# outputs a nice TSV table witg gene names and their coordinates
library(optparse)

# CLI parsing
option_list <- list(
    make_option(c("-i", "--input"),
                type = "character",
                default = NULL,
                help = "GFF file to parse",
                metavar = "path"),
    make_option(c("-o", "--output"),
                type = "character",
                default = NULL,
                help = "output file",
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

input_gff <- opt$input
output_gff <- opt$output

# Libraries
library(dplyr)
library(stringr)

# functions
read_gff <- function(gff_filename) {
  # read gff fileand retiurn it as plain data frame
  df <-
    as_tibble(read.delim(gff_filename, header = F, comment.char = "#")) %>%
    select(V1, V3, V4, V5, V7, V9, V10) %>%
    filter(V3 == "gene")
  return(df)
}

# main
gff_table <- read_gff(input_gff)
gff_annotated <- filter(gff_table, V10 > 0) %>%
  mutate(V11 = str_extract(V9, pattern="(?<=gene=).*(?=;locus_tag)"),
         V12 = str_extract(V9, pattern="(?<=ID=).*(?=_gene;)")) %>%
  select(-c(V3, V9, V10))

names(gff_annotated) <- c("contig", "start", "stop", "strand", "gene_name", "locus_tag")

# write results
write.table(gff_annotated, output_gff, sep = "\t", row.names = FALSE)
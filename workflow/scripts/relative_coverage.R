# script to read depth files and calculate relative coverage
library(optparse)

# CLI parsing
option_list <- list(
    make_option(c("-d", "--depth"),
                type = "character",
                default = NULL,
                help = "depth file",
                metavar = "depth.tsv.gz"),
    make_option(c("-r", "--ref"),
                type = "character",
                default = NULL,
                help = "reference genome file",
                metavar = "ref.fasta"),
    make_option(c("-s", "--strain"),
                type = "character",
                default = NULL,
                help = "strain name",
                metavar = "DA00000"),
    make_option(c("-l", "--label"),
                type = "character",
                default = NULL,
                help = "mutant or parent",
                metavar = "<mutant/parent>"),
    make_option(c("-m", "--min_contig_len"),
                type = "integer",
                default = 1000,
                help = "length in nt",
                metavar = "int"),
    make_option(c("-o", "--output"),
                type = "character",
                default = NULL,
                help = "output file (tsv)",
                metavar = "relative_depth.tsv")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$depth)){
 print_help(opt_parser)
 stop("Input depth file must be provided", call. = FALSE)
}
if (is.null(opt$ref)){
 print_help(opt_parser)
 stop("Input reference file must be provided", call. = FALSE)
}
if (is.null(opt$output)){
 print_help(opt_parser)
 stop("Output file must be provided", call. = FALSE)
}

#### Libraries ####
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Biostrings))

#### Functions ####
parse_depth <- function(depth_path, ref_path, min_contig_length=1000, strain="DA00000", label="no_label") {
  # param depth_path: a path to the depth file
  # param ref_path: a path to the reference genome in fasta format
  # min_contig_length: contigs to filter out
  # param strain: strain name
  # param label: mutant or parent
  # read reference
  ref_gen <- readDNAStringSet(ref_path)
  # make a table with contig names and lengths
  contig_len_df <- tibble("contig.name" = names(ref_gen),
                         "contig.len" = width(ref_gen))
  # read depth file
  depth_df <- read.table(depth_path, sep = "\t", header = FALSE)
  # check correct column names
  stopifnot(exprs = {"V1" %in% names(depth_df) & "V3" %in% names(depth_df)})
  # check if the table is empty
  stopifnot(nrow(depth_df) > 0)

  # calculate length and avg. coverage of each contig
  cov_df <- depth_df %>%
    group_by(V1) %>%
    summarise(tot.depth = sum(V3)) %>%
    dplyr::rename("contig.name" = V1) %>%
    full_join(contig_len_df, by = "contig.name") %>%
    mutate(avg.depth = tot.depth / contig.len) %>%
    arrange(-contig.len) %>%
    filter(contig.len >= min_contig_length)

  # relative coverage
  chrom_cov <- cov_df %>%
    filter(contig.len == max(contig.len)) %>%
    pull(avg.depth)

  cov_df$rel.depth <- cov_df$avg.depth / chrom_cov
  cov_df$strain <- strain
  cov_df$label <- label
  # make a column indicating what is chrom and what is plasmid
  chrom_plasm_col <- c("chrom", rep("plasmid", nrow(cov_df) - 1))
  cov_df$genome.part <- chrom_plasm_col
  
  return(cov_df)
}

#### Main ####
coverage_table <- parse_depth(opt$depth, opt$ref, opt$min_contig_len, opt$strain, opt$label) # nolint: line_length_linter.
write.table(x = coverage_table, file = opt$output, sep = "\t", row.names = FALSE)
print("Finished, no errors.")

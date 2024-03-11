# scipt to find and save best IS representatives of each family
library(optparse)

# CLI parsing
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              default = NULL,
              help = "A path to ISEscan results",
              metavar = "results/isescan"),
  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "output file name (tsv)",
              metavar = "IS_examples.tsv")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("An input directory must be provided", call. = FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("An output filename must be provided", call. = FALSE)
}

#### Libraries ####
suppressPackageStartupMessages(library(dplyr))
library(readr)

#### Functions ####
filter_margins <- function(df, filt_term, marg) {
  # filters ISEscan table
  # keeps rows with isLen & orfLen within +/- 10% of the normal IS size
  # for specified IS family
  filter(
    df,
    family == filt_term,
    between(isLen, db.length - marg * db.length, db.length + marg *
              db.length),
    between(orfLen, db.length - marg * db.length, db.length + marg *
              db.length)
  )
}

select_plausible_length <- function(fam_name, big_is_df, mrg=c(0.05, 0.1, 0.5)) {
    # filter those entries in IS table that have both
    # IS and ORF lengths within specified margins
    # try the first margin
    out_df <- filter_margins(big_is_df, fam_name, mrg[1])
    if (nrow(out_df) == 0) {
      # if nothing found
      # try another margin
      message("try margin 2...")
      out_df <- filter_margins(big_is_df, fam_name, mrg[2])
    }
    # if still 0, try another margin
    if (nrow(out_df) == 0) {
      message("try margin 3...")
      out_df <- filter_margins(big_is_df, fam_name, mrg[3])
    }
    # if still 0, issue a warning
    if (nrow(out_df) == 0) {
      message(paste0("No entries found for ", fam_name))
    }
    # select some columns
    out_df %>%
      select(seqID, family, cluster, isLen, orfLen)
}

#### Main ####

# read and join all files with IS characteristics
# collect results
isescan_results <- dir(path = opt$input, pattern = "regions_joined_final.fasta.tsv", recursive = TRUE)

# read them as TSV and bind by rows
is_df <-
  map_dfr(isescan_results, ~ tryCatch({
    read_tsv(paste0(opt$input, .), show_col_types = FALSE)
  }, error = function(e) {
    message(paste("No ISEscan results exist for ", ., sep = " "))
    return(NULL)
  }))

# know your ISfamilies
is_families <- unique(is_df$family)

# apply the functions from above, select the first row of each output table
best_is_representatives <- map_dfr(is_families, ~ select_plausible_length(., is_df) %>%
  filter(orfLen < isLen) %>%
    slice_head(n=1))

# save
write_delim(best_is_representatives, file = opt$output, delim = "\t")
message("Finished. No errors.")
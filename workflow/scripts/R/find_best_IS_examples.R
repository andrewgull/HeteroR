#################################################################
# scipt to find and save best IS representatives of each family
# input: path to ISEscan results
# output: output file name (tsv)
################################################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### LIBRARIES ####
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
library(readr)

#### FUNCTIONS ####
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

select_plausible_length <- function(fam_name,
                                    big_is_df,
                                    mrg = c(0.05, 0.1, 0.5)) {
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

join_results <- function(results_path) {
  # read and join all files with IS characteristics
  # collect results
  isescan_results <- dir(path = results_path,
                         pattern = "regions_joined_final.fasta.tsv",
                         recursive = TRUE)

  # read them as TSV and bind by rows
  is_df <-
    map_dfr(isescan_results, ~ tryCatch({
      read_tsv(paste0(results_path, .), show_col_types = FALSE)
    }, error = function(e) {
      message(paste("No ISEscan results exist for ", ., sep = " "))
      return(NULL)
    }))

  # know your ISfamilies
  is_families <- unique(is_df$family)

  # apply the functions from above, select the first row of each output table
  best_is <-
    map_dfr(is_families, ~ select_plausible_length(., is_df) %>%
              filter(orfLen < isLen) %>%
              slice_head(n = 1))
  return(best_is)
}

#### RUN ####
best_is_representatives <- join_results(snakemake@input[[1]])

# save to file
write_delim(best_is_representatives,
            file = snakemake@output[[1]],
            delim = "\t")
print("Finished. No errors.")

#### CLOSE LOG ####
sink()

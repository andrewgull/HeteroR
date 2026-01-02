#############################################################
# Script to read depth files and calculate relative coverage
#############################################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### LIBRARIES ####
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Biostrings))

#### FUNCTIONS ####
parse_depth <- function(depth_path,
                        ref_path,
                        min_con_len = 1000,
                        strain = "DA00000",
                        label = "no_label") {
  # param depth_path: a path to the depth file
  # param ref_path: a path to the reference genome in fasta format
  # min_contig_length: contigs to filter out
  # param strain: strain name
  # param label: mutant or parent
  # read reference
  ref_gen <- readDNAStringSet(ref_path)
  # make a table with contig names and lengths
  contig_len_df <- tibble(
    "contig.name" = names(ref_gen),
    "contig.len" = width(ref_gen)
  )
  # read depth file
  depth_df <- read.table(depth_path, sep = "\t", header = FALSE)
  # check correct column names
  stopifnot(exprs = {
    "V1" %in% names(depth_df) & "V3" %in% names(depth_df)
  })
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
    filter(contig.len >= min_con_len)

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

# to write the table
main <- function(cov_table, filename) {
  write.table(
    x = cov_table,
    file = filename,
    sep = "\t",
    row.names = FALSE
  )
}

#### RUN ####
coverage_table <- parse_depth(
  depth_path = snakemake@input[["depth"]],
  ref_path = snakemake@input[["ref"]],
  strain = snakemake@wildcards[["parent"]],
  label = snakemake@params[["label"]],
  min_con_len = snakemake@params[["min_len"]]
)

# save it
main(
  cov_table = coverage_table,
  filename = snakemake@output[["rel_cov"]]
)

print("Finished, no errors.")

#### CLOSE LOG ####
sink()

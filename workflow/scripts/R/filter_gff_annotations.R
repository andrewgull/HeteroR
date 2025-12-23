#####################################################################
#
# A script for reading GFF file and filtering
# based on the 10th column + retrieving gene names.
# This column contains the fraction of GFF file covered by VCF file.
# If values in this column greater than 0
# then the corresponding gene contains a variant.
# The script outputs a nice TSV table with gene names and their coordinates.
#
#####################################################################

#### OPEN LOG ####
sink(snakemake@log[[1]])

#### LIBRARIES ####
suppressPackageStartupMessages(library(dplyr))
library(stringr)

#### FUNCTIONS ####
# read gff file and return it as a plain data frame
read_gff <- function(gff_filename) {
  df <-
    as_tibble(read.delim(gff_filename, header = FALSE, comment.char = "#")) %>%
    select(V1, V3, V4, V5, V7, V9, V10) %>%
    filter(V3 == "gene")
  return(df)
}

# main function
gff2tsv <- function(input_gff, output_tsv) {
  # check if annotated GFF is empty (in case when there was nothing to annotate)
  if (file.size(input_gff) == 0L) {
    # make an empty DF
    gff_annotated <- data.frame("contig" = c(NA), "start" = c(NA),
                                "stop" = c(NA), "strand" = c(NA),
                                "gene_name" = c(NA), "locus_tag" = c(NA))
    print("An empty GFF file is detected.")
  } else {
    # parse annotated GFF
    gff_table <- read_gff(input_gff)
    gff_annotated <- filter(gff_table, V10 > 0) %>%
      mutate(V11 = str_extract(V9, pattern = "(?<=gene=).*(?=;locus_tag)"),
             V12 = str_extract(V9, pattern = "(?<=ID=).*(?=_gene;)")) %>%
      select(-c(V3, V9, V10))

    names(gff_annotated) <- c("contig", "start", "stop", 
                              "strand", "gene_name", "locus_tag")
    print("GFF file is annotated.")
  }

  # write results to the output
  write.table(gff_annotated, output_tsv, sep = "\t", row.names = FALSE)
  print("GFF file is converted successfully.")
}

#### RUN ####
gff2tsv(snakemake@input[[1]], snakemake@output[[1]])
print("Script completed successfully.")

#### CLOSE LOG ####
sink()

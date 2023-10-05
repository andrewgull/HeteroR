# The alignment file for the iqtree was initially built over all strains. 
# This code is to remove the strains that were not used for modelling. 

# libs
library(Biostrings)

# constants
annotations_dir <- "~/Data/HeteroR/results/pangenome/annotations_db/"
features_data <- "data/features_strain.csv"
lab_testing_data <- "data/heteroresistance_testing.csv"
file_alignment_raw <- "~/Data/HeteroR/results/pangenome/roary_results_gtr/core_gene_alignment.aln"
file_alignment_clean <- "~/Data/HeteroR/results/pangenome/roary_results_gtr/core_gene_alignment_clean.aln"

# find all the strains that were used for the alignment 
files <- list.files(annotations_dir)
files <- files[grep("^DA", files)] 
files <- gsub(".gff", "", files)

# find list of strains that were used for modelling
data_strain <- readr::read_csv(features_data, na = c("NA", "-Inf"), show_col_types = F) 

# HR labels
hr_testing <- readr::read_csv(lab_testing_data, show_col_types = F)

data_strain <- dplyr::left_join(hr_testing, data_strain, by = "strain")

# find strains to remove 
if (length(setdiff(files, data_strain$strain)) > 0) {
  to_remove <- paste0(files[!files %in% data_strain$strain], "_genomic")
  raw_alignment <-readDNAStringSet(file_alignment_raw)
  
  # remove
  clean_alignment <-  raw_alignment[!names(raw_alignment) %in% to_remove]
  clean_alignment
  
  # save
  writeXStringSet(clean_alignment, file_alignment_clean)
} else {
  print("Nothing to remove from the alignment file")
}

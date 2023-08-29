# The alignment file for the iqtree was initially built over all strains. 
# This code is to remove the strains that were not used for modelling. 

# libs
library(Biostrings)

# find all the strains that were used for the alignment 
files <- list.files("~/HeteroR/results/pangenome/annotations_db/")
files <- files[grep("^DA", files)] 
files <- gsub(".gff", "", files)

# find list of strains that were used for modelling
data_strain <- read_csv("data/features_strain.csv", na = c("NA", "-Inf")) 

# HR labels
hr_testing12 <- read_csv("data/heteroresistance_testing.csv", col_select = c(strain, Gr12)) %>% 
  filter(!is.na(strain)) %>% 
  rename("resistance" = Gr12 )

data_strain <- data_strain %>% 
  left_join(hr_testing12, by = "strain") %>% 
  filter(resistance != "R", !is.na(resistance)) %>% pull(strain)

# find strains to remove 
to_remove <- paste0(files[!files %in% data_strain], "_genomic")
raw_alignment <-readDNAStringSet("~/HeteroR/results/pangenome/roary_results_gtr/core_gene_alignment.aln")

# remove
clean_alignment <-  raw_alignment[!names(raw_alignment) %in% to_remove]
clean_alignment

# save
writeXStringSet(clean_alignment, '~/HeteroR/results/pangenome/roary_results_gtr/core_gene_alignment_clean.aln')

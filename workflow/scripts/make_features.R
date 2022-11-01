# make features table for ML trainig and validation

# libraries
library(optparse)
library(tidyverse)

# CLI options
option_list <- list(
   make_option(c("-b", "--bed"),
               type = "character",
               default = NULL,
               help = "A path to BED files",
               metavar = "character"),
   make_option(c("-r", "--repeat_tables"),
               type = "character",
               default = NULL,
               help = "A path to tables with repeat pairs (csv)",
               metavar = "character"),
   make_option(c("-g", "--resistance_genes"),
               type = "character",
               default = NULL,
               help = "A path to tables with resistance genes (csv)",
               metavar = "character"),
   make_option(c("-t", "--resistance_testing"),
               type = "character",
               default = NULL,
               help = "A path to tables with resistance testing results (csv)",
               metavar = "character"),
   make_option(c("-a", "--assembly_summaries"),
               type = "character",
               default = NULL,
               help = "A path to assemblies summaries tables (tsv)",
               metavar = "character"),
   make_option(c("-o", "--out"),
                type = "character",
                default = NULL,
                help = "output file name",
                metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$bed)) {
 print_help(opt_parser)
 stop("bed input files must be provided", call. = FALSE)
}
if (is.null(opt$table)) {
 print_help(opt_parser)
 stop("repeat table input files must be provided", call. = FALSE)
}

# SHORTCUTS FOR INPUTS/OUTPUTS
bed_file_paths <- opt$bed # character vector from expand()
repeat_csv_paths <- opt$table # character vector from expand()
rgi_csv_paths <- opt$resistance_genes # character vector from expand()
testing_csv <- opt$resistance_testing # one file
assembly_summaries <- opt$assembly_summaries # character vector from expand()
output_file <- opt$out # a tsv file
rgi_identity <- 90
rgi_length_percent <- 90
target_rg <- "beta-lactamase"

###################
### READ INPUTS ###
###################

# READ REPEAT TABLES
repeat_df <- bind_rows(lapply(repeat_csv_paths, function(x) {
    read_csv(x)
}))

# READ BED FILES
# to determine center spanning repeat pairs
bed_df <- bind_rows(lapply(bed_file_paths, function(x) {
  df <- read_delim(x, header = FALSE)
  # parse strain name from input path
  df$strain <- "strain"
  return(df)
}))

# READ RGI DATA
rgi <- bind_rows(lapply(rgi_csv_paths, function(x) {
  df <- read_delim(x, na.strings = "n/a")
  # parse strain name from input path
  df$strain <- "strain"
  return(df)
}))

# GET AB GROUPS
ab_groups <- unique(rgi$Drug.Class)

# READ LAB TESTING RESULTS
resistance_df <- read_csv(testing_csv)
resistance_df$resistance <- as.factor(resistance_df$resistance)

# READ ASSEMBLY SUMMARIES
assembly_summary <- lapply(assembly_summaries, function(x) {
  strain <- "strain function here"
  df <- read_delim(x)
  df$Links <- as.integer(gsub(",", "", df$Links, fixed = TRUE))
  df$Length <- as.integer(gsub(",", "", df$Length, fixed = TRUE))
  df$N50 <- as.integer(gsub(",", "", df$N50, fixed = TRUE))
  df$Longest_component <- as.integer(gsub(",", "", df$Longest_component, fixed = TRUE)) # nolint
  return(df)
})

# READ GFF FILES
# a function to read GFF files using STANDARD read.delim
read_gff <- function(gff_filename){
  # read and process
  df <- as_tibble(read.delim(gff_filename, header=F, comment.char="#")) %>%
  select(V1, V3, V4, V5, V7, V9) %>% 
  filter(V3=="gene")
  return(df)
  }


#############################
### SOME INPUT PROCESSING ###
#############################

# bind assembly summaries to a single file
assembly_summary <- bind_rows(assembly_summary)

# calculate plasmid counts
plasmid_counts <- assembly_summary %>%
  group_by(Strain, Type) %>%
  summarise(n.plasmids = n()) %>%
  filter(Type == "Plasmid") %>%
  select(-Type) %>%
  rename("strain" = Strain)

# calculate AR length
repeat_df$AR_length <- repeat_df$start_2 - repeat_df$end_1 + 1

# find gene center in bed files
bed_df$gene_center <- round((bed_df$V3 - bed_df$V2 + 1) / 2)
bed_df <- rename(bed_df, "record_id" = V4)
bed_df <- select(bed_df, strain, record_id, gene_center)

# bed_df contains slightly more genes than repeat_df (+624)
# BUT repeat_df is important, that's why I use left_join() below

# join BED and Repeats
repeat_df <- left_join(repeat_df, bed_df, by = "record_id")
repeat_df$spans_center <- if_else(
    repeat_df$end_1 <= repeat_df$gene_center &
    repeat_df$start_2 >= repeat_df$gene_center,
    "yes", "no")
repeat_df$record_id <- sub("_gene", "", repeat_df$record_id)
repeat_df <- repeat_df %>%
  rename("repeat_length" = length) %>%
  select(-strain.y) %>%
  rename("strain" = strain.x) %>%
  relocate(strain, .before = record_id)

# RGI filtering
# add record id and put it in front of the table
rgi$record_id <- map_chr(rgi$ORF_ID, function(x) strsplit(x, " ")[[1]][1])
rgi <- relocate(rgi, strain, record_id, .before = ORF_ID)

# filter out Loose hits and some columns
rgi <- rgi %>%
  filter(Cut_Off != "Loose") %>%
  select(-c(Predicted_Protein, CARD_Protein_Sequence,
            Contig, Start, Stop, Orientation, Predicted_DNA, Note))

# filter by Identity
rgi <- rgi %>% filter(Best_Identities >= rgi_identity)

# filter by Percentage Length
rgi <- filter(rgi, Percentage.Length.of.Reference.Sequence >= rgi_length_percent) # nolint

# filter by Resistance Mechanism
stay <- c("antibiotic target replacement", "antibiotic target protection", "antibiotic inactivation") # nolint
rgi <- filter(rgi, Resistance.Mechanism %in% stay)

# filter by SNP
rgi <- filter(rgi, is.na(SNPs_in_Best_Hit_ARO))

##############################
### MAKE BASIC FEATURES DF ###
##############################
features <- repeat_df %>%
  select(-c(gene_center, start_1, end_1, start_2, end_2, X))

# join it with filtered RGI
features <- full_join(features, rgi, by = c("strain", "record_id"))

# remove some columns
features <- features %>%
  select(-c(ARO, Cut_Off, Pass_Bitscore, Best_Hit_Bitscore, Best_Hit_ARO, 
            Best_Identities, Model_ID, Model_type, Nudged, SNPs_in_Best_Hit_ARO,
            Other_SNPs, ID, Percentage.Length.of.Reference.Sequence))

# filter genes without Resistance.Mechanism
features <- filter(features, !is.na(Resistance.Mechanism))

# join features and testing results
features <- left_join(features, resistance_df, by = "strain")

# join features and plasmid counts
features <- left_join(features, plasmid_counts, by = "strain")
# NA in n.plasmids means 0
features$n.plasmids <- ifelse(is.na(features$n.plasmids), 0, features$n.plasmids) # nolint
features <- features %>% relocate(n.plasmids, .before = "ORF_ID")

##########################
### BL LABELS & COUNTS ###
##########################

# label beta-lactamases (or other RGs of interest)
f_bl_filt <- features %>%
  select(strain, record_id, AMR.Gene.Family) %>%
  distinct() %>%
  mutate(is_beta_lac = grepl(target_rg, AMR.Gene.Family))
# if mutate() doesn't work make it on a separate line

# count BL genes (or other type)
n_beta_lac <- f_bl_filt %>%
  group_by(strain) %>%
  summarise(n.beta.lac = sum(is_beta_lac))

# add is.beta.lac to features for future calculations
features <- left_join(features, f_bl_filt %>%
  select(record_id, is_beta_lac), by = "record_id") %>%
    relocate(is_beta_lac, .before = "ORF_ID")

###########################
### FEATURES PER STRAIN ###
###########################
features_strain <- features %>%
  select(strain, resistance, n.plasmids) %>%
  distinct()

# join BL counts with features_strain
features_strain <- left_join(features_strain, n_beta_lac, by = "strain")

###########################
### LOCATION ON PLASMID ###
###########################

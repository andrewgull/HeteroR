# make features table for ML trainig and validation

# libraries
library(optparse)
library(tidyverse)

#### CLI options ####
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
   make_option(c("-f", "--gff_files"),
               type = "character",
               default = NULL,
               help = "GFF files names",
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

#### Assign values from opt ####
bed_file_paths <- opt$bed # character vector from expand()
repeat_csv_paths <- opt$table # character vector from expand()
rgi_csv_paths <- opt$resistance_genes # character vector from expand()
testing_csv <- opt$resistance_testing # one file
assembly_summaries <- opt$assembly_summaries # character vector from expand()
gff_files <- opt$gff_files # character vector from expand()
output_file <- opt$out # a tsv file
rgi_identity <- 90
rgi_length_percent <- 90
target_rg <- "beta-lactamase"

#### Read inputs ####

# READ REPEAT TABLES
repeat_df <- bind_rows(lapply(repeat_csv_paths, function(x) {
    read_csv(x)
}))

# READ BED FILES
# to determine center spanning repeat pairs
bed_df <- bind_rows(lapply(bed_file_paths, function(x) {
  df <- read_delim(x, col_names = FALSE)
  # parse strain name from input path
  df$strain <- "strain"
  return(df)
}))

# READ RGI DATA
rgi <- bind_rows(lapply(rgi_csv_paths, function(x) {
  df <- read_delim(x, na = "n/a")
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
read_gff <- function(gff_filename) {
  # read and process
  # read_delim() didn't work
  df <- read.delim(gff_filename, header = FALSE, comment.char = "#") %>%
    dplyr::select(V1, V3, V4, V5, V7, V9) %>%
    dplyr::filter(V3 == "gene")
  return(tibble::as_tibble(df))
}
# read gff files, it takes time
gff <- map_dfr(gff_files, ~ read_gff(.))


#### Input processing ####

# bind assembly summaries to a single file
assembly_summary <- bind_rows(assembly_summary)

# PLASMID COUNTS
plasmid_counts <- assembly_summary %>%
  group_by(Strain, Type) %>%
  summarise(n.plasmids = n()) %>%
  filter(Type == "Plasmid") %>%
  select(-Type) %>%
  rename("strain" = Strain)

# AR LENGTH
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

# RGI FILTERING
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

#### Basic features DF ####
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

#### Beta-Lac labels & counts ####

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

#### Features per strain ####
features_strain <- features %>%
  select(strain, resistance, n.plasmids) %>%
  distinct()

# join BL counts with features_strain
features_strain <- left_join(features_strain, n_beta_lac, by = "strain")

#### Location on plasmid ####
# process GFFs
# transform the last column to get the same record_id as in features
gff <- separate(gff, V9, sep = ";locus_tag=", into = c("V10", "V11")) %>%
  select(-V10)
# give new names
names(gff) <- c("seq.region", "type", "start", "stop", "strand", "record_id")
# filter out genes that are not in features
gff_in_features <- gff %>%
  filter(record_id %in% unique(features$record_id))
# remove these >2 million rows
rm(gff)
# add 'located on plasmid': if there is '1' in seq.region then its a chromosome
# (the first element of an assembly)
gff_in_features$on.plasmid <- if_else(str_detect(gff_in_features$seq.region, "1", negate = T), 1, 0)
#gff_in_features$on.plasmid <- if_else(gff_in_features$located.on.plasmid, 1, 0)

#### Distance to oriC ####
gff_in_features <- rename(gff_in_features, "dist.to.oriC" = start)

#### On plus strand ####
gff_in_features$pos.strand <- if_else(gff_in_features$strand == "+", 1, 0)

# remove some columns
gff_in_features <- select(gff_in_features, c(record_id, dist.to.oriC, pos.strand, on.plasmid))

# add is beta lac
gff_in_features <- left_join(gff_in_features, select(features, c(record_id, is_beta_lac)) %>% distinct(), by = "record_id")
gff_bl <- gff_in_features %>%  filter(is_beta_lac) %>% select(-is_beta_lac)

#### Add dist to oriC and plus strand to features ####
features <- left_join(features, gff_bl, by="record_id")
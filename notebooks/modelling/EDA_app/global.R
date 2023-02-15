# Shiny web application for EDA

library(shiny)
library(tidyverse)
library(plotly)
library(recipes)
library(bestNormalize)
library(ggforce)

# Read the main data table with features
features_strain <- read_csv("data/features_strain.csv") 

# get vars to use later in UI (selectInput)
vars <- names(select(features_strain, -c(strain, resistance)))
num_vars <- names(select(features_strain, -c(strain, resistance, n.beta.lac.3, n.beta.lac.4, chrom.status)))
char_vars <- c("resistance", "n.beta.lac.3", "n.beta.lac.4", "chrom.status")

# get strain names to use in UI (heat map)
strains <- features_strain$strain

# get ampC and non-ampC counts 
bl_count <- features_strain %>% 
  select(resistance, n.beta.lac, ampC.type.beta.lactamase) %>% 
  mutate(non.ampC = n.beta.lac - ampC.type.beta.lactamase) %>% 
  select(-n.beta.lac)

# make it tidy
bl_count_tidy <- gather(bl_count, key = "gene", value = "n", 2:3) %>% 
  group_by(resistance, gene) %>% 
  summarize(sum = sum(n))

###################
## NORMALIZATION ##
###################

data_rec <- recipe(resistance ~., data = features_strain) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_orderNorm(all_numeric_predictors()) 

data_prep <- prep(data_rec, retain = TRUE)

features_norm <- data_prep$template

vars_norm <- names(features_norm %>% select(-strain, resistance))

###########
# HEATMAP #
###########

# get per gene data
data <- read_csv("../data/features.csv")

# get unique AMR types
amr_types <- unique(data$AMR.Gene.Family)

# get unique gene rows with AMR.Gene.Family
bl_amr_gene_types <- data %>% 
  select(strain, record_id, AMR.Gene.Family) %>% 
  distinct()

# make a DF of each AMR type counts in each gene
bl_blr_gene_presence <- as_tibble(
  as.data.frame(
    map(amr_types, function(x) ifelse(grepl(x, bl_amr_gene_types$AMR.Gene.Family, fixed = T), 1, 0)),
    col.names = amr_types))

# add strain and summarize counts
bl_blr_gene_presence$strain <- bl_amr_gene_types$strain
bl_blr_gene_presence <- relocate(bl_blr_gene_presence, strain, .before = "ampC.type.beta.lactamase")
# here sum is applied to all variables
bl_amr_types_strain <- bl_blr_gene_presence %>% 
  group_by(strain) %>% 
  summarise_each(~sum(.))
# we need shorter names
#names(bl_amr_types_strain) <- c("strain", "ampC", "DFR", "APH6", "APH3.1", "SUL", "TEM", "SAT", "ANT3", "MPH", "APH3.2", "CTX.M", "CAT", "AAC3", "OXA", "AAC6", "ANT2", "FTT", "SHV", "TR.RPP", "APH4", "QNR")

# make tidy for plotting
bl_amr_types_strain_td <- gather(bl_amr_types_strain, 
                                 key = "AMR.type", 
                                 value = "N", 
                                 2:ncol(bl_amr_types_strain))

##########
## UMAP ##

library(embed)

# a special function
plot_validation_results <- function(dat) {
  dat %>%
    select(-strain) %>% 
    # Create the scatterplot matrix
    ggplot(aes(x = .panel_x, y = .panel_y, color = resistance, fill = resistance)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_autodensity(alpha = .3) +
    facet_matrix(vars(-resistance), layer.diag = 2) + 
    scale_color_brewer(palette = "Set1") + 
    scale_fill_brewer(palette = "Set1")
}

# read UMAP 12 data with falsely identified strains
data_umap12_postEDA <- read_csv("data/data_umap12_postEDA.csv")

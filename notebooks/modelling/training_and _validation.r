# a script to run training and validatio of selected ML models
# using workflowsets package

# packages
library(tidymodels)
library(workflowsets)
library(bonsai)
library(themis)
library(colino)
library(bestNormalize)
library(baguette)
library(lightgbm)
library(tidyposterior)
library(ggplot2)
source("functions.R")

# TOC
read_data()
process_data()
make_recipe()
make_spec()
run_workflow()
save_diagnostic_plots()

# this functions reads two files: with features per strain, must contain column 'strain'
# and with laboratory testing results, must contain columns 'strain' and 'resistance'
read_data <- function(file_path_1 = "data/features_strain.csv",
                      file_path_2 = "data/heteroresistance_testing.csv") {
  # check if 'strain' is present in both files
  # check that both files contain more than one column and more than zero rows
  data_strain <- readr::read_csv(file_path_1,
                                 na = c("NA", "-Inf"),
                                 show_col_types = FALSE)
  
  hr_testing <- readr::read_csv(file_path_2,
                                show_col_types = FALSE)
  
  data_strain <- data_strain %>%
    left_join(hr_testing, by = "strain")
  
  return(data_strain)
}
# process the data the same way as in model_analysis.Rmd
data_strain <- readr::read_csv("data/features_strain.csv",
                               na = c("NA", "-Inf"),
                               show_col_types = FALSE)

hr_testing <- readr::read_csv("data/heteroresistance_testing.csv",
                              show_col_types = FALSE)

data_strain <- data_strain %>%
  left_join(hr_testing, by = "strain")

# strain filtering: DA63310 doesn't have ampC and repeats
# the other two clustered with E. fergusoni and E. albertii on the tree
data_strain <- data_strain %>%
  mutate(n.beta.lac.4 = factor(ifelse(n.beta.lac > 4, "yes", "no"))) %>%
  relocate(n.beta.lac.4, .before = "n.plasmids") %>%
  filter(resistance != "R", 
         !(strain %in% c("DA63310", "DA63246", "DA63068"))) %>%
  mutate(
    resistance = factor(resistance, levels = c("HR", "nonHR")),
    chrom.status = factor(chrom.status)
  )

data_strain$N50 <- NULL
data_strain[is.na(data_strain)] <- 0
data_strain <- data_strain %>% 
  select(-contains("oriC")) %>% 
  select(-contains("plus"))

#### SPLIT ####
# same seed number as in modelling.Rmd
set.seed(124)

# splitting proportion should be the same
data_split <- initial_split(data_strain,
                            prop = 0.8,
                            strata = resistance)

df_train <- training(data_split)
df_test <- testing(data_split)

cv_folds <- vfold_cv(df_train,
                     strata = "resistance",
                     v = 10,
                     repeats =10)

# metrics for imbalanced classes
imbalanced_metrics <- metric_set(roc_auc, j_index, mcc, pr_auc)


cores <- 8
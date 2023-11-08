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

#### TOC ####
read_input()
process_input()
make_recipe()
make_spec()
run_workflow()
save_diagnostic_plots()

#### MINOR FUNCTIONS ####
# to check if a provided DF contains column 'strain'
column_exists <- function(df, colname) {
  if (colname %in% names(df)){
    TRUE
  }  else {
    FALSE
  }
}

# to check if a proivided DF has no rows
empty <- function(df) {
  if (nrow(df) == 0) {
    TRUE
  } else {
    FALSE
  }
}


#### MAIN FUNCTIONS ####
# this functions reads two files: with features per strain, must contain column 'strain'
# and with laboratory testing results, must contain columns 'strain' and 'resistance'
read_input <- function(file_path_1 = "data/features_strain.csv",
                      file_path_2 = "data/heteroresistance_testing.csv") {
  # file_path_1: a path to file with features per strain
  # file_path_2: a path to file with testing results (outcome)
  # return: left join between 1 and 2
  
  # read the first file
  data_strain <- readr::read_csv(file_path_1,
                                 na = c("NA", "-Inf"),
                                 show_col_types = FALSE)
  
  # check if column 'strain' is present
  if (!column_exists(data_strain, "strain")){
    message("Input data does not contain column 'strain'")
    stop()
  }
  
  # check if df is empty
  if (empty(data_strain)) {
    message("One of the provided tables is empty!")
    stop()
  }
  
  # read the second file
  hr_testing <- readr::read_csv(file_path_2,
                                show_col_types = FALSE)
  
  # check if column 'strain' is present
  if (!column_exists(hr_testing, "strain")){
    message("Input data does not contain column 'strain'")
    stop()
  }
  
  # check if df is empty
  if (empty(hr_testing)) {
    message("One of the provided tables is empty!")
    stop()
  }
  
  # join the two data frames
  data_strain <- data_strain %>%
    left_join(hr_testing, by = "strain")
  
  return(data_strain)
}


# this function filters data and adds/removes columns
process_input <- function(df){
  # df: a data frame with features, outcome and strain names
  # return: filtered df with new/removed columns
  
  # check if df is empty
  if (empty(df)) {
    message("Provided table with features is empty!")
    stop()
  }
  
  # strain filtering: DA63310 doesn't have ampC and repeats
  # the other two clustered with E. fergusoni and E. albertii on the tree
  df <- df %>%
    mutate(n.beta.lac.4 = factor(ifelse(n.beta.lac > 4, "yes", "no"))) %>%
    relocate(n.beta.lac.4, .before = "n.plasmids") %>%
    filter(resistance != "R", 
           !(strain %in% c("DA63310", "DA63246", "DA63068"))) %>%
    mutate(
      resistance = factor(resistance, levels = c("HR", "nonHR")),
      chrom.status = factor(chrom.status)
    ) %>% 
    select(-contains("oriC"), -contains("plus"), -N50) %>% 
    replace(is.na(.), 0)
  
  return(df)
}



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
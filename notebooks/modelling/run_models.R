# run models and save them
library(tidymodels)
library(themis)
library(probably)
library(vip)
library(skimr)
library(stacks)
library(bestNormalize) # for ord QQ norm
library(embed)

# DATA
data_strain <- readr::read_csv("/home/andrei/GitProjects/HeteroR/notebooks/modelling/data/features_strain.csv", na = c("NA", "-Inf"))

hr_testing <- readr::read_csv("/home/andrei/GitProjects/HeteroR/notebooks/modelling/data/heteroresistance_testing_gr12.csv")

data_strain <- data_strain %>% 
  left_join(hr_testing, by = "strain")

data_strain <- data_strain %>% 
  mutate(n.beta.lac.3 = factor(ifelse(n.beta.lac > 3, "yes", "no"))) %>% 
  mutate(n.beta.lac.4 = factor(ifelse(n.beta.lac > 4, "yes", "no"))) %>% 
  relocate(n.beta.lac.3, n.beta.lac.4, .before = "n.plasmids") %>% 
  filter(resistance != "R") %>% 
  mutate(resistance = factor(resistance, levels = c("HR", "nonHR")),
         chrom.status = factor(chrom.status))

data_strain$N50 <- NULL
data_strain$NA. <- NULL
data_strain[is.na(data_strain)] <- 0

# SPLIT
set.seed(124)

data_split <- initial_split(data_strain, prop = 0.8, strata = resistance)

df_train <- training(data_split)
df_test <- testing(data_split)

# RECIPES
main_recipe <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>% 
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100)


lr_recipe <- main_recipe %>%
  step_corr(threshold = 0.75)

# Stratified, repeated 10-fold cross-validation is used to resample the model:
cv_folds <- vfold_cv(df_train, strata = "resistance", v = 10, repeats = 10) # is better than v=5
# metrics for imbalanced classes
cls_metrics <- metric_set(roc_auc, j_index)


# LR + Bayes
set.seed(333)
# set model type/engine
lr_mod <- 
  logistic_reg(
    penalty = tune(), 
    mixture = 1) %>% 
  set_engine("glmnet")

# define the workflow
lr_workflow <- 
  workflow() %>% 
  add_model(lr_mod) %>% 
  add_recipe(lr_recipe)

# extract settings
lr_set <- extract_parameter_set_dials(lr_workflow)

set.seed(12)

lr_bres <-
  lr_workflow %>% 
  tune_bayes(
    resamples = cv_folds,
    # To use non-default parameter ranges
    param_info = lr_set,
    # Generate five at semi-random to start
    initial = 5,
    iter = 50,
    # How to measure performance?
    metrics = metric_set(roc_auc),
    control = control_bayes(no_improve = 30, verbose = FALSE, save_pred = TRUE, save_workflow = TRUE)
  )

save.image("/home/andrei/Data/HeteroR/notebooks/lrb.RData")

# run models and save them
library(optparse)

# CLI parsing
option_list <- list(
  make_option(c("-m", "--model"),
              type = "character",
              default = NULL,
              help = "A model type (one of the following; lr, mars, svm, rf, knn, bt)",
              metavar = "character"),
  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "output file (.RData)",
              metavar = "character"),
  make_option(c("-t", "--threads"),
              type = "integer",
              default = 2,
              help = "number of threads for Random Forest and Boosted Trees",
              metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

model_type <- opt$model
output_file <- opt$output
threads <- opt$threads

# libs
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


non_corr_recipe <- main_recipe %>%
  step_corr(threshold = 0.75)

# Stratified, repeated 10-fold cross-validation is used to resample the model:
cv_folds <- vfold_cv(df_train, strata = "resistance", v = 10, repeats = 10) # is better than v=5
# metrics for imbalanced classes
cls_metrics <- metric_set(roc_auc, j_index)

###########################################################
# evaluate models with resampling and Bayesian grid search

set.seed(333)
# set model type/engine

set_model <- function(mod, cores) {
  if (mod == "lr") {
    my_mod <- logistic_reg(
        penalty = tune(), 
        mixture = 1) %>% 
      set_engine("glmnet")
  } else if (mod == "mars") {
    my_mod <- mars(
      mode = "classification",
      engine = "earth",
      num_terms = tune(),
      prod_degree = tune(),
      prune_method = "backward") %>% 
      translate()
  } else if(mod == "svm") {
    my_mod <- svm_linear(
        cost = tune()) %>% # margin - for regression only
      set_mode("classification") %>%
      set_engine("kernlab")  # default
  } else if(mod == "rf") {
    my_mod <- rand_forest(
      mtry = tune(), 
      min_n = tune(), 
      trees = 1000) %>% 
      set_engine("ranger", num.threads = cores) %>% 
      set_mode("classification")
  } else if(mod == "bt") {
    my_mod <- boost_tree(
      trees = 50, 
      mtry = tune(), 
      min_n = tune(), 
      tree_depth = tune(), 
      learn_rate = tune(), 
      loss_reduction = tune(), 
      sample_size = tune(), 
      stop_iter = tune()) %>% 
      set_engine("xgboost", num.threads = cores) %>% 
      set_mode("classification")
  } else if(mod == "knn") {
    my_mod <- nearest_neighbor(
      neighbors = tune(),
      weight_func = tune(),
      dist_power = tune()) %>% #
      set_engine("kknn") %>%
      set_mode("classification")
  }
  return(my_mod)
}

# define the model
my_model <- set_model(mod = model_type, cores = threads)

# make a workflow with no_corr recipe
if (model_type == "rf" | model_type == "bt"){
  my_workflow <- workflow() %>% 
    add_model(my_model) %>% 
    add_recipe(main_recipe)
  
param_set <- extract_parameter_set_dials(my_workflow) %>%
  finalize(x = df_train %>% select(-resistance))

} else {
  my_workflow <- workflow() %>% 
    add_model(my_model) %>% 
    add_recipe(non_corr_recipe)
  
  # extract settings for bayesian grid search
  param_set <- extract_parameter_set_dials(my_workflow)
}



set.seed(12)

# resample and evaluate
model_bres <-
  my_workflow %>% 
  tune_bayes(
    resamples = cv_folds,
    # To use non-default parameter ranges
    param_info = param_set,
    # Generate five at semi-random to start
    initial = 5,
    iter = 50,
    # How to measure performance?
    metrics = metric_set(roc_auc),
    control = control_bayes(no_improve = 30, 
                            verbose = FALSE, 
                            save_pred = TRUE, 
                            save_workflow = TRUE)
  )

# save
save.image(output_file)

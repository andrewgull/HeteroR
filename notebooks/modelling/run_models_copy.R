# run models and save them
library(optparse)

#### CLI parsing ####
option_list <- list(
  make_option(c("-m", "--model"),
              type = "character",
              default = NULL,
              help = "A model type (should be one of the following: lr, mars, svm, rf, knn, bt)",
              metavar = "character"),
  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "output file (.rds)",
              metavar = "character"),
  make_option(c("-t", "--threads"),
              type = "integer",
              default = 2,
              help = "number of threads for Random Forest and Boosted Trees",
              metavar = "integer"),
  make_option(c("-r", "--recipe"),
              type = "character",
              default = NULL,
              help = "recipe to use (should be one of the follwing: main, ncorr, pca, umap)",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

#### LIBS ####
suppressPackageStartupMessages(library(tidymodels)) # to keep quiet
library(themis) # for SMOTE
library(bestNormalize) # for ORQ-norm
library(embed) # for UMAP

#### DATA ####
path_data <- "/home/andrei/GitProjects/HeteroR/notebooks/modelling/data/features_strain.csv"
path_labels <- "/home/andrei/GitProjects/HeteroR/notebooks/modelling/data/heteroresistance_testing_gr12.csv"
data_strain <- readr::read_csv(path_data, 
                              na = c("NA", "-Inf"),
                              show_col_types = FALSE)

hr_testing <- readr::read_csv(path_labels, 
                              show_col_types = FALSE)

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

#### SPLIT ####
set.seed(124)

data_split <- initial_split(data_strain, prop = 0.8, strata = resistance)

df_train <- training(data_split)
df_test <- testing(data_split)

#### RECIPES ####
main_recipe <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>% 
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100)

ncorr_recipe <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>% 
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_corr(threshold = 0.75) %>%
  step_smote(resistance, over_ratio = 1, seed = 100)

pca_recipe <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_orderNorm(all_numeric_predictors()) %>% 
  step_normalize(all_predictors()) %>%
  step_pca(all_numeric_predictors(), num_comp = tune()) %>% 
  step_normalize(all_numeric_predictors()) %>% 
  step_smote(resistance, over_ratio = 1, seed = 100)

umap_recipe <- recipe(resistance ~., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_orderNorm(all_numeric_predictors()) %>% 
  step_normalize(all_predictors()) %>% 
  step_umap(all_numeric_predictors(), outcome = "resistance", num_comp = 20) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) 

#### FOLDS & METRICS ####
cv_folds <- vfold_cv(df_train, 
                     strata = "resistance", 
                     v = 10, 
                     repeats = 10) # is better than v=5

cls_metrics <- metric_set(roc_auc, j_index) # metrics for imbalanced classes

#### FUN: MODEL SPECIFICATION ####
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
      dist_power = tune()) %>%
      set_engine("kknn") %>%
      set_mode("classification")
  }
  return(my_mod)
}

#### FUN: WORKFLOW
set_wf <- function(mod, rec, cores){
  # mod: model type, one of: lr, knn, mars, svm, rf, bt
  # rec: recipe object (one of: main, ncorr, pca, umap)
  # rec must be in GlobalEnv
  
  if (rec == "main"){
    rc <- main_recipe
  } else if (rec == "ncorr"){
    rc <- ncorr_recipe
  } else if (rec == "pca") {
    rc <- pca_recipe
  } else if (rec == "umap") {
    rc <- umap_recipe
  } else {
    print(" ERROR! Undefined recipe!")
  }

  wf <- workflow() %>% 
    add_model(set_model(mod = mod, cores = cores)) %>% 
    add_recipe(rc)
  
  return(wf)
}

#### CREATE A WORKFLOW ####
my_wf <- set_wf(mod = opt$model, 
                rec = opt$recipe, 
                cores = opt$threads)

#### EXTRACT PARAMETERS ####
if (opt$model == "rf" | opt$model == "bt"){
  # extract settings for bayesian grid search (special case)
  param_set <- extract_parameter_set_dials(my_wf) %>%
    finalize(x = df_train %>% select(-resistance))
} else {
  # extract settings for bayesian grid search
  param_set <- extract_parameter_set_dials(my_wf)
}

#### MODEL TUNING ####
if (opt$model == "lr"){
  print("Space-filling grid search is chosen...")
  model_res <- my_wf %>%
          tune_grid(
              grid = 30,
              resamples = cv_folds,
              control = control_grid(save_pred = TRUE, save_workflow = TRUE),
              metrics = cls_metrics)
} else {
  print("Bayesian grid search is chosen...")
  model_res <- my_wf %>% 
  tune_bayes(
    resamples = cv_folds,
    # To use non-default parameter ranges
    param_info = param_set,
    # Generate N at semi-random to start
    initial = 8,
    iter = 50,
    # How to measure performance?
    metrics = metric_set(roc_auc),
    control = control_bayes(no_improve = 30, 
                            verbose = FALSE, 
                            save_pred = TRUE, 
                            save_workflow = TRUE)
  )
}


#### SAVE MODEL ####
saveRDS(object = model_res, file = opt$output)


# a script to run training and validatio of selected ML models
# using workflowsets package

# packages
library(tidymodels)
library(workflowsets)
library(bonsai)
library(themis)
library(bestNormalize)
library(baguette)
library(lightgbm)
library(colino)


#### MINOR FUNCTIONS ####
# to check if a provided DF contains column 'strain'
column_exists <- function(df, colname) {
  if (colname %in% names(df)){
    TRUE
  }  else {
    FALSE
  }
}

# to check if a provided DF has no rows
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
process_input <- function(df) {
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

# this function returns requested recipe
make_recipe <- function(name = "base", df) {
  # rec: recipe name, one of 'base', 'pca', base_yj', 'base_orq', 'ncorr',
  # 'ncorr_orq', 'ncorr_yj', 'boruta'
  # df: train data set
  if (name == "base") {
    rec <- recipe(resistance ~ ., data = df) %>%
      update_role(strain, new_role = "ID") %>%
      step_nzv(all_predictors()) %>%
      step_normalize(all_numeric_predictors()) %>%
      step_dummy(all_nominal_predictors()) %>%
      step_nzv(all_predictors()) %>%
      step_smote(resistance, over_ratio = 1, seed = 100)
  } else if (name == "pca") {
    rec <- recipe(resistance ~ ., data = df) %>%
      update_role(strain, new_role = "ID") %>%
      step_nzv(all_predictors()) %>%
      step_dummy(all_nominal_predictors()) %>%
      step_YeoJohnson(all_numeric_predictors()) %>%
      step_normalize(all_numeric_predictors()) %>%
      step_smote(resistance, over_ratio = 1, seed = 100) %>%
      step_pca(all_predictors(), num_comp = tune()) %>%
      step_normalize(all_numeric_predictors())
  } else if (name == "base_yj") {
    rec <- recipe(resistance ~ ., data = df) %>%
      update_role(strain, new_role = "ID") %>%
      step_nzv(all_predictors()) %>%
      step_dummy(all_nominal_predictors()) %>%
      step_nzv(all_predictors()) %>%
      step_YeoJohnson(all_numeric_predictors()) %>%
      step_normalize(all_numeric_predictors()) %>%
      step_smote(resistance, over_ratio = 1, seed = 100)
  } else if (name == "base_orq") {
    rec <- recipe(resistance ~ ., data = df) %>%
      update_role(strain, new_role = "ID") %>%
      step_nzv(all_predictors()) %>%
      step_dummy(all_nominal_predictors()) %>%
      step_orderNorm(all_numeric_predictors()) %>%
      step_normalize(all_numeric_predictors()) %>%
      step_smote(resistance, over_ratio = 1, seed = 100)
  } else if (name == "ncorr") {
    rec <- recipe(resistance ~ ., data = df) %>%
      update_role(strain, new_role = "ID") %>%
      step_nzv(all_predictors()) %>%
      step_normalize(all_numeric_predictors()) %>%
      step_dummy(all_nominal_predictors()) %>%
      step_nzv(all_predictors()) %>%
      step_smote(resistance, over_ratio = 1, seed = 100) %>% 
      step_corr(all_predictors(), threshold = tune("corr_tune"))
  } else if (name == "ncorr_yj") {
    rec <- recipe(resistance ~ ., data = df) %>%
      update_role(strain, new_role = "ID") %>%
      step_nzv(all_predictors()) %>%
      step_normalize(all_numeric_predictors()) %>%
      step_dummy(all_nominal_predictors()) %>%
      step_nzv(all_predictors()) %>%
      step_YeoJohnson(all_numeric_predictors()) %>%
      step_normalize(all_numeric_predictors()) %>%
      step_smote(resistance, over_ratio = 1, seed = 100) %>% 
      step_corr(all_predictors(), threshold = tune("corr_tune"))
  } else if (name == "ncorr_orq") {
    rec <- recipe(resistance ~ ., data = df) %>%
      update_role(strain, new_role = "ID") %>%
      step_nzv(all_predictors()) %>%
      step_normalize(all_numeric_predictors()) %>%
      step_dummy(all_nominal_predictors()) %>%
      step_nzv(all_predictors()) %>%
      step_orderNorm(all_numeric_predictors()) %>%
      step_normalize(all_numeric_predictors()) %>%
      step_smote(resistance, over_ratio = 1, seed = 100) %>% 
      step_corr(all_predictors(), threshold = tune("corr_tune"))
  } else if (name == "bor") {
    rec <- recipe(resistance ~ ., data = df) %>%
      update_role(strain, new_role = "ID") %>%
      step_nzv(all_predictors()) %>%
      step_normalize(all_numeric_predictors()) %>%
      step_dummy(all_nominal_predictors()) %>%
      step_nzv(all_predictors()) %>%
      step_smote(resistance, over_ratio = 1, seed = 100) %>% 
      step_select_boruta(all_predictors(), outcome = "resistance")
  } else {
    message("Provided recipe does not exist")
    stop()
  }
  return(rec)
}

# this function returns requested model specification
make_spec <- function(mod = "lr", threads = 1){
  # param mod: a two to four letter code for a model
  # param threads: number of threads
  # return: model specification
  # lr = logistic regression
  # lsvm = linear SVM
  # psvm = polynomial SVM
  # rsvm = radial basis function SVM
  # rf = random forest
  # bt = gradient boosted trees
  # mlp = multilayer perceptron
  # mlpb = multilayer perceptron with bagging
  if (mod == "lr") {
    my_mod <- logistic_reg(
      penalty = tune(),
      mixture = 1, 
      mode = "classification", 
      engine = "glmnet")
  } else if(mod == "lsvm") {
    my_mod <- svm_linear(
      cost = tune(), 
      mode = "classification", 
      engine = "kernlab",
      num.threads = threads)
  } else if(mod == "rf") {
    my_mod <- rand_forest(
      mtry = tune(),
      min_n = tune(),
      trees = 1000, 
      mode = "classification", 
      engine = "ranger",
      num.threads = threads)
  } else if (mod == "bt") {
    my_mod <- boost_tree(
      trees = tune(),
      min_n = tune(),
      tree_depth = tune(),
      learn_rate = tune(),
      loss_reduction = tune(),
      sample_size = tune(),
      stop_iter = tune(),
      engine = "lightgbm",
      mode = "classification",
      num.threads = threads) %>% 
      translate()
  } else if (mod == "psvm") {
    my_mod <- svm_poly(
      cost = tune(),
      degree = tune(),
      scale_factor = tune(),
      margin = NULL,
      engine="kernlab",
      mode = "classification",
      num.threads = threads)
  } else if (mod == "rbfsvm") {
    my_mod <- svm_rbf(
      cost = tune(),
      rbf_sigma = tune(),
      margin = tune(),
      engine = "kernlab",
      mode = "classification",
      num.threads = threads) 
  } else if (mod == "mlp") {
    my_mod <-
      mlp(hidden_units = tune(),
          penalty = tune(),
          epochs = tune(),
          engine = "nnet",
          mode = "classification")
  } else if (mod == "mlpb") {
    my_mod <- bag_mlp(
      hidden_units = tune(),
      penalty = tune(),
      epochs = tune(),
      engone = "nnet",
      mode = "classification",
      num.threads = threads) %>%
      translate()
  }
  return(my_mod)
}

# this function uses workflowsets to run training and validation
run_workflow_sets <- function(model_specs_list, recipes_list, folds, metrics_set_obj, grid_size = 30, seed = 124) {
    # param file_path: a path to file to save results to
    # param model_specs_list: list of model specs
    # param recipes_list: list of recipes
    # param folds: CV folds object
    # param metrics_set_obj: metric_set() object
    # param grid_size: size of the grid fro grid search
    # param seed: seed number
    # return: workflowset object
   models_set <-
        workflow_set(preproc = recipes_list,
                     models = model_specs_list,
                     cross = TRUE)
      set.seed(seed)
      models_set <-
        models_set %>%
        workflow_map(
          "tune_grid",
          resamples = folds,
          grid = grid_size,
          metrics = metrics_set_obj,
          verbose = TRUE,
          control = control_grid(save_pred = TRUE,
                                 save_workflow = TRUE)
        )
    return(models_set)
}

main <- function(data, model_names, recipes_names, file_path, folds, metric){


  model_specs <- lapply(model_names, \(x) make_spec(mod = x, threads = 8))
  recipes <- lapply(recipes_names, \(x) make_recipe(name = x, df = data))

  if (!file.exists(file_path)) {
     print("File does not exist. Proceeding with training and validation.")
     models_set <- run_workflow_sets(model_specs_list = model_specs,
                                     folds = cv_folds,
                                     recipes_list = recipes,
                                     metrics_set_obj = metric)
     saveRDS(object = models_set, file =  file_path)
  } else {
     print("File exists. Skipping training and validation")
  }

}

# If one wants Bayesian optimization of BT results
main_bayes <- function(data, file_path, grid_search_res_file, folds) {
  
  if (!file.exists(file_path)) {
    print("File does not exist. Proceeding with optimizing the models.")
    
    base_bt_wf <- workflow() %>%
      add_model(make_spec("bt")) %>%
      add_recipe(make_recipe(name = "base", df = data))
    
    param_set_base_bt <-
      extract_parameter_set_dials(base_bt_wf) %>%
      finalize(x = data %>% select(-resistance))
    
    models_bt <- readRDS(grid_search_res_file)
    
    base_bt_bres <-
      tune_bayes(
        base_bt_wf,
        resamples = folds,
        initial = models_bt$result[[1]],
        iter = 40,
        metrics = metric_set(roc_auc),
        param_info = param_set_base_bt,
        control = control_bayes(
          no_improve = 20,
          save_pred = TRUE,
          verbose = FALSE,
          save_workflow = TRUE
        )
      )
    saveRDS(object = base_bt_bres, file = file_path)
  } else {
    base_bt_bres <- readRDS(file_path)
  }
  return(base_bt_bres)
}

#### RUN TRAINING AND VALIDAITON ####

# prepare data sets and folds
# same seed number as in modelling.Rmd
set.seed(124)

# splitting proportion should be the same
data_strain <- read_input()
data_strain <- process_input(data_strain)
data_split <- initial_split(data_strain,
                            prop = 0.8,
                            strata = resistance)

df_train <- training(data_split)

cv_folds <- vfold_cv(df_train,
                     strata = "resistance",
                     v = 10,
                     repeats =10)

# metrics for imbalanced classes
imbalanced_metrics <- metric_set(roc_auc, j_index, mcc, pr_auc)

# LR
main(data = df_train, 
     model_names = "lr",
     recipes_names = c("base", "base_yj", "base_orq", "pca"), 
     file = "/home/andrei/Data/HeteroR/results/models/lr_resamples.rds", 
     folds = cv_folds, 
     metric = imbalanced_metrics)

# SVM
main(data = df_train, 
     model_names = c("lsvm", "psvm", "rbfsvm"), 
     recipes_names = c("ncorr", "ncorr_yj", "ncorr_orq", "pca"), 
     file = "/home/andrei/Data/HeteroR/results/models/svm_resamples.rds", 
     folds = cv_folds, 
     metric = imbalanced_metrics)

# MLP
main(data = df_train,
     model_names = c("mlp", "mlpb"), 
    recipes_names = c("ncorr", "ncorr_yj", "ncorr_orq", "pca"), 
    file = "/home/andrei/Data/HeteroR/results/models/mlp_resamples.rds", 
    folds = cv_folds, 
    metric = imbalanced_metrics)

# RF
main(data = df_train,model_names = "rf", 
    recipes_names = c("base", "bor"), 
    file = "/home/andrei/Data/HeteroR/results/models/rf_resamples.rds", 
    folds = cv_folds, 
    metric = imbalanced_metrics)

# BT
main(data = df_train,
     model_names = "bt", 
     recipes_names = c("base", "bor"), 
     file = "/home/andrei/Data/HeteroR/results/models/bt_resamples.rds", 
     folds = cv_folds, 
     metric = imbalanced_metrics)

# BT + Bayes
main_bayes(data = df_train, 
           file_path = "/home/andrei/Data/HeteroR/results/models/bt_bayes_resamples.rds", 
           grid_search_res_file = "/home/andrei/Data/HeteroR/results/models/bt_resamples.rds",
           folds = cv_folds)

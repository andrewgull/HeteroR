# this script makes a table for comparing ROC, J-index and their SE for different ML algorithms

library(readr)
library(tidymodels)
library(themis)
library(probably)
library(vip)
library(skimr)
library(stacks)
library(bestNormalize) # for ord QQ norm
library(tidyposterior) # for Bayesian ANOVA
library(finetune) # for win-loss tuning
source("functions.R")
library(colino)

# path for models
models_path <-
  "~/HeteroR/results/models/scheme12/"

rfe_model <- rand_forest(mode = "classification") %>%
  set_engine("ranger", num.threads = 8, importance = "impurity")

# function to get j index
get_j_index <-
  function(res) {
    res %>% show_best("j_index", n = 40) %>% filter(.config == (
      res %>% show_best("roc_auc", n = 1) %>% select(.config) %>%  as.character()
    )) %>%  select(mean) %>% as.numeric()
  }

get_j_index_SE <-
  function(res) {
    res %>% show_best("j_index", n = 40) %>% filter(.config == (
      res %>% show_best("roc_auc", n = 1) %>% select(.config) %>%  as.character()
    )) %>%  select(std_err) %>% as.numeric()
  }


# Read and process data

# The data sets are the same as in EDA

data_strain <-
  read_csv("data/features_strain.csv", na = c("NA", "-Inf")) %>%
  filter(strain != "DA63310")

# HR testing lables
hr_testing12 <-
  read_csv("data/heteroresistance_testing.csv",
           col_select = c(strain, resistance)) %>%
  filter(!is.na(strain))

data_strain <- data_strain %>%
  left_join(hr_testing12, by = "strain")

data_strain %>%
  group_by(resistance) %>%
  count()

#Add n.beta.lac \>4

strains <- data_strain$strain

data_strain <- data_strain %>%
  #mutate(n.beta.lac.3 = factor(ifelse(n.beta.lac > 3, "yes", "no"))) %>%
  mutate(n.beta.lac.4 = factor(ifelse(n.beta.lac > 4, "yes", "no"))) %>%
  relocate(n.beta.lac.4, .before = "n.plasmids") %>%
  filter(resistance != "R") %>%
  mutate(
    resistance = factor(resistance, levels = c("HR", "nonHR")),
    chrom.status = factor(chrom.status)
  )

data_strain$N50 <- NULL
data_strain[is.na(data_strain)] <- 0

# Data split

set.seed(124)

data_split <-
  initial_split(data_strain, prop = 0.8, strata = resistance)

df_train <- training(data_split)
df_test <- testing(data_split)

cv_folds <-
  vfold_cv(df_train,
           strata = "resistance",
           v = 10,
           repeats = 10)

# RECIPIES

base_recipe <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100)

pca_recipe <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_YeoJohnson(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>%
  step_pca(all_predictors(), num_comp = tune()) %>%
  step_normalize(all_numeric_predictors())

base_yj_recipe <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_YeoJohnson(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100)

base_orq_recipe <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_orderNorm(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100)

ncorr_recipe <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>%
  step_corr(all_predictors(), threshold = tune("corr_tune"))


ncorr_yj_recipe <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_YeoJohnson(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>%
  step_corr(all_predictors(), threshold = tune("corr_tune"))


ncorr_orq_recipe <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_orderNorm(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>%
  step_corr(all_predictors(), threshold = tune("corr_tune"))


base_boruta_recipe <- base_recipe %>%
  step_select_boruta(all_predictors(), outcome = "resistance")

# MODELS
lr_en <- logistic_reg(penalty = tune(),
                      mixture = tune()) %>%
  set_engine("glmnet")

#LR
models_llr <-
  readRDS(paste0(models_path,"models_llr.rds"))

predictors_to_remove <- c("(Intercept)", "resistance", "strain")

final_tibble <- tibble(
  Model = character(0),
  Preprocessing = character(0),
  mean_AUC = numeric(0),
  AUC_SE = numeric(0),
  mean_j_index = numeric(0),
  j_index_SE = numeric(0),
  Correlation = numeric(0),
  n_predictors = integer(0)
)


new_model_row <-
  function(model_results,
           model_name,
           preprocessing_name) {
    best_perf_hyp <- model_results %>% show_best(n = 20)
    
    mean_auc <- best_perf_hyp$mean[1]
    
    AUC_SE <- best_perf_hyp$std_err[1]
    
    # bayesian optimized BT does not have a j-index
    if (model_name == "BT_bres") {
      mean_j_index <- NA
      
      j_index_SE <- NA
    } else {
      mean_j_index <- get_j_index(model_results)
      
      j_index_SE <- get_j_index_SE(model_results)
    }
    
    corr_tune <- best_perf_hyp$corr_tune[1]
    
    
    
    # find number of predictors according to the model
    # for LR we have to fit the model since Lasso removes features.
    # for other models only the preprocessing is enough to count the number of predictors.
    
    models_except_LLR = c(
      "lSVM",
      "pSVM",
      "rbfSVM",
      "MLP_nnet",
      "BAG_MLP",
      "MARS",
      "BAG_MARS",
      "KNN",
      "RF",
      "BT",
      "BT_bres"
    )
    
    if (model_name == "LR") {
      #select preprocessing
      if (preprocessing_name == "BASE") {
        my_rec = base_recipe
      } else if (preprocessing_name == "BASE_yj") {
        my_rec = base_yj_recipe
      } else if (preprocessing_name == "BASE_orq") {
        my_rec = base_orq_recipe
      } else if (preprocessing_name == "PCA") {
        my_rec = pca_recipe
      }
      my_mod <-  logistic_reg(penalty =  best_perf_hyp$penalty[1],
                              mixture = 1) %>%
        set_engine("glmnet")
      
      my_wf <- workflow() %>%
        add_model(my_mod) %>%
        add_recipe(my_rec)
      
      best_params <- model_results %>%
        select_best("roc_auc")
      
      final_model <- finalize_workflow(my_wf, best_params) %>%
        last_fit(data_split)
      
      trained_mod <- final_model %>%
        extract_workflow() %>%
        pull_workflow_fit()
      
      coefs <- coef(trained_mod$fit, s = best_perf_hyp$penalty[1])
      
      coefficients <- coefs %>% row.names()
      n_predictors <-
        coefficients[coefs@i] %>% discard( ~ .x %in% predictors_to_remove) %>% length
    } else if (model_name %in% models_except_LLR) {
      if (preprocessing_name == "NCORR") {
        n_predictors <-  recipe(resistance ~ ., data = df_train) %>%
          update_role(strain, new_role = "ID") %>%
          step_nzv(all_predictors()) %>%
          step_normalize(all_numeric_predictors()) %>%
          step_dummy(all_nominal_predictors()) %>%
          step_nzv(all_predictors()) %>%
          step_smote(resistance, over_ratio = 1, seed = 100) %>%
          step_corr(all_predictors(), threshold = corr_tune) %>%
          prep() %>%
          juice() %>%
          colnames() %>%
          discard(~ .x %in% predictors_to_remove) %>%
          length()
      } else if (preprocessing_name == "NCORR_YJ") {
        n_predictors <-  recipe(resistance ~ ., data = df_train) %>%
          update_role(strain, new_role = "ID") %>%
          step_nzv(all_predictors()) %>%
          step_normalize(all_numeric_predictors()) %>%
          step_dummy(all_nominal_predictors()) %>%
          step_nzv(all_predictors()) %>%
          step_YeoJohnson(all_numeric_predictors()) %>%
          step_normalize(all_numeric_predictors()) %>%
          step_smote(resistance, over_ratio = 1, seed = 100) %>%
          step_corr(all_predictors(), threshold = corr_tune) %>%
          prep() %>% juice() %>% colnames()  %>% discard( ~ .x %in% predictors_to_remove) %>% length()
      } else if (preprocessing_name == "NCORRORQ") {
        n_predictors <-   recipe(resistance ~ ., data = df_train) %>%
          update_role(strain, new_role = "ID") %>%
          step_nzv(all_predictors()) %>%
          step_normalize(all_numeric_predictors()) %>%
          step_dummy(all_nominal_predictors()) %>%
          step_nzv(all_predictors()) %>%
          step_orderNorm(all_numeric_predictors()) %>%
          step_normalize(all_numeric_predictors()) %>%
          step_smote(resistance, over_ratio = 1, seed = 100) %>%
          step_corr(all_predictors(), threshold = corr_tune) %>%
          prep() %>% juice() %>% colnames()  %>% discard( ~ .x %in% predictors_to_remove) %>% length()
      } else if (preprocessing_name == "PCA") {
        n_predictors <- recipe(resistance ~ ., data = df_train) %>%
          update_role(strain, new_role = "ID") %>%
          step_nzv(all_predictors()) %>%
          step_dummy(all_nominal_predictors()) %>%
          step_YeoJohnson(all_numeric_predictors()) %>%
          step_normalize(all_numeric_predictors()) %>%
          step_smote(resistance, over_ratio = 1, seed = 100) %>%
          step_pca(all_predictors(), num_comp = best_perf_hyp$num_comp[1]) %>%
          step_normalize(all_numeric_predictors()) %>%
          prep() %>% juice() %>% colnames()  %>% discard( ~ .x %in% predictors_to_remove) %>% length()
      } else if (preprocessing_name == "BASE") {
        n_predictors <-   recipe(resistance ~ ., data = df_train) %>%
          update_role(strain, new_role = "ID") %>%
          step_nzv(all_predictors()) %>%
          step_normalize(all_numeric_predictors()) %>%
          step_dummy(all_nominal_predictors()) %>%
          step_nzv(all_predictors()) %>%
          step_smote(resistance, over_ratio = 1, seed = 100) %>% prep() %>% juice() %>% colnames()  %>% discard( ~ .x %in% predictors_to_remove) %>% length()
      }   else if (preprocessing_name == "BASE_BORUTA") {
        n_predictors <-    recipe(resistance ~ ., data = df_train) %>%
          update_role(strain, new_role = "ID") %>%
          step_nzv(all_predictors()) %>%
          step_normalize(all_numeric_predictors()) %>%
          step_dummy(all_nominal_predictors()) %>%
          step_nzv(all_predictors()) %>%
          step_smote(resistance, over_ratio = 1, seed = 100) %>%
          step_select_boruta(all_predictors(), outcome = "resistance") %>% prep() %>% juice() %>% colnames()  %>% discard( ~ .x %in% predictors_to_remove) %>% length()
        
      }
    }
    
    return(
      tibble(
        Model = model_name,
        Preprocessing = preprocessing_name,
        mean_AUC = mean_auc,
        AUC_SE = AUC_SE,
        mean_j_index = mean_j_index,
        j_index_SE = j_index_SE,
        Correlation = corr_tune,
        n_predictors =  n_predictors
      )
    )
  }


## BASE

final_tibble <-
  bind_rows(final_tibble, new_model_row(models_llr$result[[1]], "LR", "BASE"))

## baseyj_LLR

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_llr$result[[2]], "LR", "BASE_yj"))

## baseorq_LLR

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_llr$result[[3]], "LR", "BASE_orq"))

## PCA LLR

final_tibble <-
  bind_rows(final_tibble, new_model_row(models_llr$result[[4]], "LR", "PCA"))


# lSVM
models_svm <-
  readRDS(paste0(models_path,"models_svm.rds"))

## ncorr_lsvm

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_svm$result[[1]], "lSVM", "NCORR"))


## ncorryj_lsvm

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_svm$result[[4]], "lSVM", "NCORR_YJ"))

## ncorrorq_lsvm

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_svm$result[[7]], "lSVM", "NCORRORQ"))


## pca_lsvm

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_svm$result[[10]], "lSVM", "PCA"))

# pSVM

## ncorr_psvm

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_svm$result[[2]], "pSVM", "NCORR"))

## ncorryj_psvm

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_svm$result[[5]], "pSVM", "NCORR_YJ"))

## ncorrorq_psvm

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_svm$result[[8]], "pSVM", "NCORRORQ"))

## pca_psvm

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_svm$result[[11]], "pSVM", "PCA"))


# rbfSVM

## ncorr_rbfsvm

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_svm$result[[3]], "rbfSVM", "NCORR"))

## ncorryj_rbfsvm

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_svm$result[[6]], "rbfSVM", "NCORR_YJ"))

## ncorrorq_rbfsvm

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_svm$result[[9]], "rbfSVM", "NCORRORQ"))

## pca_rbfsvm

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_svm$result[[12]], "rbfSVM", "PCA"))


# MLP_nnet

models_mlp <-
  readRDS(paste0(models_path,"models_mlp.rds"))


## ncorr_mlp_nnet

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_mlp$result[[1]], "MLP_nnet", "NCORR"))


## ncorryj_mlp_nnet

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_mlp$result[[3]], "MLP_nnet", "NCORR_YJ"))

## ncorrorq_mlp_nnet

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_mlp$result[[5]], "MLP_nnet", "NCORRORQ"))

## pca_mlp_nnet

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_mlp$result[[7]], "MLP_nnet", "PCA"))

# BAG_MLP

## ncorr_bag_mlp

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_mlp$result[[2]], "BAG_MLP", "NCORR"))

## ncorryj_bag_mlp

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_mlp$result[[4]], "BAG_MLP", "NCORR_YJ"))

## ncorrorq_bag_mlp

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_mlp$result[[6]], "BAG_MLP", "NCORRORQ"))

## pca_bag_mlp

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_mlp$result[[8]], "BAG_MLP", "PCA"))



# MARS
models_mars <-
  readRDS(paste0(models_path,"models_mars.rds"))

## ncorr_mars

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_mars$result[[1]], "MARS", "NCORR"))

## ncorryj_mars

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_mars$result[[3]], "MARS", "NCORR_YJ"))

## ncorrorq_mars

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_mars$result[[5]], "MARS", "NCORRORQ"))

## pca_mars

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_mars$result[[7]], "MARS", "PCA"))

# BAG MARS

## ncorr_bag_mars

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_mars$result[[2]], "BAG_MARS", "NCORR"))

## ncorryj_bag_mars

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_mars$result[[4]], "BAG_MARS", "NCORR_YJ"))

## ncorrorq_bag_mars

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_mars$result[[6]], "BAG_MARS", "NCORRORQ"))

## pca_bag_mars

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_mars$result[[8]], "BAG_MARS", "PCA"))

# KNN

models_knn <-
  readRDS(paste0(models_path,"models_knn.rds"))

## ncorr_knn

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_knn$result[[1]], "KNN", "NCORR"))

## ncorryj_knn

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_knn$result[[2]], "KNN", "NCORR_YJ"))

## ncorrorq_knn

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_knn$result[[3]], "KNN", "NCORRORQ"))

## pca_knn

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_knn$result[[4]], "KNN", "PCA"))

# RF
models_rf <-
  readRDS(paste0(models_path,"models_rf.rds"))


## base_rf

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_rf$result[[1]], "RF", "BASE"))


## base_boruta_rf

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_rf$result[[2]], "RF", "BASE_BORUTA"))

# BT

models_bt <-
  readRDS(paste0(models_path,"models_bt.rds"))

## base_bt

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_bt$result[[1]], "BT", "BASE"))

## base_bt_bres
base_bt_bres <-
  readRDS(paste0(models_path,"models_bt_bres.rds"))


final_tibble <-
  bind_rows(final_tibble, new_model_row(base_bt_bres, "BT_bres", "BASE"))


## base_boruta_bt

final_tibble <-
  bind_rows(final_tibble,
            new_model_row(models_bt$result[[2]], "BT", "BASE_BORUTA"))


## base_boruta_bt_bres
# bt_base_boruta_bres <-
#   readRDS(paste0(models_path,"base_bt_boruta_bres.rds"))
# 
# 
# final_tibble <-
#   bind_rows(final_tibble,
#             new_model_row(bt_base_boruta_bres, "BT_bres", "BASE_BORUTA"))

# save
write_csv(final_tibble, "data/model_comparison_table.csv")

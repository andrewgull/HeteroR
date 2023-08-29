
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
models_path <- "/mnt/data/andrei/Data/HeteroR/results/models/scheme12/"

rfe_model <- rand_forest(mode = "classification") %>%
  set_engine("ranger", num.threads = 8, importance = "impurity")

# function to get j index
get_j_index <- function(res){res %>% show_best("j_index",n =40) %>% filter(.config == (
    res %>% show_best("roc_auc", n = 1) %>% select(.config) %>%  as.character()
)) %>%  select(mean) %>% as.character()}

get_j_index_SE <- function(res){res %>% show_best("j_index",n =40) %>% filter(.config == (
    res %>% show_best("roc_auc", n = 1) %>% select(.config) %>%  as.character()
)) %>%  select(std_err) %>% as.character()}


# Read and process data

# The data sets are the same as in EDA

data_strain <- read_csv("data/features_strain.csv", na = c("NA", "-Inf")) %>% 
  filter(strain != "DA63310")

# HR testing lables
hr_testing12 <- read_csv("data/heteroresistance_testing.csv", col_select = c(strain, Gr12)) %>% 
  filter(!is.na(strain)) %>% 
  dplyr::rename("resistance" = Gr12 )

data_strain <- data_strain %>% 
  left_join(hr_testing12, by = "strain")

data_strain %>% 
  group_by(resistance) %>% 
  count()

Add n.beta.lac \>4

strains <- data_strain$strain

data_strain <- data_strain %>% 
  #mutate(n.beta.lac.3 = factor(ifelse(n.beta.lac > 3, "yes", "no"))) %>% 
  mutate(n.beta.lac.4 = factor(ifelse(n.beta.lac > 4, "yes", "no"))) %>% 
  relocate(n.beta.lac.4, .before = "n.plasmids") %>% 
  filter(resistance != "R") %>% 
  mutate(resistance = factor(resistance, levels = c("HR", "nonHR")),
         chrom.status = factor(chrom.status))

data_strain$N50 <- NULL
data_strain[is.na(data_strain)] <- 0

# Data split

set.seed(124)

data_split <- initial_split(data_strain, prop = 0.8, strata = resistance)

df_train <- training(data_split)
df_test <- testing(data_split)

cv_folds <- vfold_cv(df_train, strata = "resistance", v = 10, repeats = 10) 

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
lr_en <- logistic_reg(
      penalty = tune(), 
      mixture = tune()) %>%
      set_engine("glmnet")

#LR 
models_llr <- readRDS("/mnt/data/andrei/Data/HeteroR/results/models/scheme12/models_llr.rds")

## BASE 
lr_base <- models_llr$result[[1]]

best_perf_hyp <- lr_base %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(lr_base)

j_index_SE <- get_j_index_SE(lr_base)

my_mod <-  logistic_reg(
        penalty =  best_perf_hyp$penalty[1],
        mixture = 1) %>%
      set_engine("glmnet")

my_wf <- workflow() %>%
  add_model(my_mod) %>%
  add_recipe(base_recipe)
 
best_params <- lr_base %>% 
    select_best("roc_auc")
 

final_model <- finalize_workflow(my_wf, best_params) %>%
  last_fit(data_split)

trained_mod <- final_model %>% 
  extract_workflow() %>% 
  pull_workflow_fit()

coefs <- coef(trained_mod$fit, s = best_perf_hyp$penalty[1])

predictors_to_remove <- c("(Intercept)","resistance", "strain")

coefficients <- coefs%>% row.names()
n_predictors <- coefficients[coefs@i] %>% discard(~ .x %in% predictors_to_remove) %>% length 

final_tibble <- tibble(Model = "LR", Preprocessing= "BASE", mean_AUC = mean_auc, AUC_SE = AUC_SE, mean_j_index = mean_j_index, j_index_SE=j_index_SE, Correlation = NA, n_predictors =  n_predictors)

final_tibble

## baseyj_LLR

lr_baseyj_res <-models_llr$result[[2]]

best_perf_hyp <- lr_baseyj_res %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(lr_baseyj_res)

j_index_SE <- get_j_index_SE(lr_baseyj_res)

corr_tune <-best_perf_hyp$corr_tune[1]

# Rebuild the workflow to extract the number of selected features
my_mod <-  logistic_reg(
        penalty =  best_perf_hyp$penalty[1],
        mixture = 1) %>%
      set_engine("glmnet")

 my_wf <- workflow() %>%
  add_model(my_mod) %>%
  add_recipe(base_yj_recipe)
 
best_params <- lr_baseyj_res %>% 
    select_best("roc_auc")
 

final_model <- finalize_workflow(my_wf, best_params) %>%
  last_fit(data_split)

trained_mod <- final_model %>% 
  extract_workflow() %>% 
  pull_workflow_fit()

coefs <- coef(trained_mod$fit, s = best_perf_hyp$penalty[1])

coefficients <- coefs %>% row.names()
n_predictors <- coefficients[coefs@i] %>% discard(~ .x %in% predictors_to_remove) %>% length 

tmp_tibble <- tibble(Model = "LR", Preprocessing= "BASE_yj", mean_AUC = mean_auc, AUC_SE = AUC_SE, mean_j_index = mean_j_index, j_index_SE= j_index_SE, Correlation = corr_tune, n_predictors =  n_predictors)

final_tibble <- bind_rows(final_tibble,tmp_tibble)

final_tibble

## baseorq_LLR
baseorq_LLR <- models_llr$result[[3]]

best_perf_hyp <- baseorq_LLR %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(baseorq_LLR)

j_index_SE <- get_j_index_SE(baseorq_LLR)

j_index_SE <- get_j_index_SE(lr_base)

corr_tune <-best_perf_hyp$corr_tune[1]

# Rebuild the workflow to extract the number of selected features 

my_mod <-  logistic_reg(
        penalty =  best_perf_hyp$penalty[1],
        mixture = 1) %>%
      set_engine("glmnet")

 my_wf <- workflow() %>%
  add_model(my_mod) %>%
  add_recipe(base_orq_recipe)
 
best_params <- baseorq_LLR %>% 
    select_best("roc_auc")
 

final_model <- finalize_workflow(my_wf, best_params) %>%
  last_fit(data_split)

trained_mod <- final_model %>% 
  extract_workflow() %>% 
  pull_workflow_fit()

coefs <- coef(trained_mod$fit, s = best_perf_hyp$penalty[1])

coefficients <- coefs%>% row.names()
n_predictors <- coefficients[coefs@i] %>% discard(~ .x %in% predictors_to_remove) %>% length 

tmp_tibble <- tibble(Model = "LR", Preprocessing= "BASE_orq", mean_AUC = mean_auc, AUC_SE = AUC_SE, mean_j_index = mean_j_index,j_index_SE=j_index_SE, Correlation = corr_tune, n_predictors =  n_predictors)

final_tibble <- bind_rows(final_tibble,tmp_tibble)

final_tibble

## PCA LLR

pca_LLR <- models_llr$result[[4]]

best_perf_hyp <- pca_LLR %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <- get_j_index(pca_LLR)

j_index_SE <- get_j_index_SE(pca_LLR)

corr_tune <- NA

# Rebuild the workflow to extract the number of selected features 

my_mod <-  logistic_reg(
        penalty =  best_perf_hyp$penalty[1],
        mixture = 1) %>%
      set_engine("glmnet")

pca_recipe <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_YeoJohnson(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>%
  step_pca(all_predictors(), num_comp = best_perf_hyp$num_comp[1]) %>%
  step_normalize(all_numeric_predictors())
 
 my_wf <- workflow() %>%
  add_model(my_mod) %>%
  add_recipe(pca_recipe)
 
best_params <- pca_LLR %>% 
    select_best("roc_auc")
 
final_model <- finalize_workflow(my_wf, best_params) %>%
  last_fit(data_split)

trained_mod <- final_model %>% 
  extract_workflow() %>% 
  pull_workflow_fit()

coefs <- coef(trained_mod$fit, s = best_perf_hyp$penalty[1])

coefficients <- coefs%>% row.names()
n_predictors <- coefficients[coefs@i] %>% discard(~ .x %in% predictors_to_remove) %>% length

tmp_tibble <- tibble(Model = "LR", Preprocessing= "PCA", mean_AUC = mean_auc, AUC_SE = AUC_SE, mean_j_index = mean_j_index,j_index_SE=j_index_SE, Correlation = corr_tune, n_predictors =  n_predictors)

final_tibble <- bind_rows(final_tibble,tmp_tibble)

final_tibble

# lSVM
models_svm_mlp <- readRDS("/mnt/data/andrei/Data/HeteroR/results/models/scheme12/models_svm_mlp.rds")

## ncorr_lsvm 
lsvm_NCORR <-models_svm_mlp$result[[1]]

best_perf_hyp <- lsvm_NCORR %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(lsvm_NCORR)

j_index_SE <- get_j_index_SE(lsvm_NCORR)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  discard( ~ .x %in% predictors_to_remove) %>% 
  length()

tmp_tibble <-
  tibble(
    Model = "lSVM",
    Preprocessing = "NCORR",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

final_tibble

## ncorryj_lsvm
lsvm_NCORR_YJ <-models_svm_mlp$result[[6]]

best_perf_hyp <- lsvm_NCORR_YJ %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(lsvm_NCORR_YJ)

j_index_SE <- get_j_index_SE(lsvm_NCORR_YJ)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <- tibble(Model = "lSVM", Preprocessing= "NCORR_YJ", mean_AUC = mean_auc, AUC_SE = AUC_SE, mean_j_index = mean_j_index,j_index_SE=j_index_SE, Correlation = corr_tune, n_predictors =  n_predictors)

final_tibble <- bind_rows(final_tibble,tmp_tibble)

final_tibble

## ncorrorq_lsvm
lsvm_NCORRORQ <-models_svm_mlp$result[[11]]

best_perf_hyp <- lsvm_NCORRORQ %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(lsvm_NCORRORQ)

j_index_SE <- get_j_index_SE(lsvm_NCORRORQ)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "lSVM",
    Preprocessing = "NCORRORQ",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

final_tibble

## pca_lsvm 
lsvm_pca <-models_svm_mlp$result[[16]]

best_perf_hyp <- lsvm_pca %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(lsvm_pca)

j_index_SE <- get_j_index_SE(lsvm_pca)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_YeoJohnson(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>%
  step_pca(all_predictors(), num_comp = best_perf_hyp$num_comp[1]) %>%
  step_normalize(all_numeric_predictors()) %>% 
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "lSVM",
    Preprocessing = "PCA",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

final_tibble

# pSVM 

## ncorr_psvm 

psvm_NCORR <-models_svm_mlp$result[[2]]

best_perf_hyp <- psvm_NCORR %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(psvm_NCORR)

j_index_SE <- get_j_index_SE(psvm_NCORR)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <-  recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>% 
  step_corr(all_predictors(), threshold = corr_tune) %>% prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <- tibble(Model = "pSVM", Preprocessing= "NCORR", mean_AUC = mean_auc, AUC_SE = AUC_SE, mean_j_index = mean_j_index,j_index_SE=j_index_SE, Correlation = corr_tune, n_predictors =  n_predictors)

final_tibble <- bind_rows(final_tibble,tmp_tibble)

final_tibble

## ncorryj_psvm
psvm_NCORR_YJ <-models_svm_mlp$result[[7]]

best_perf_hyp <- psvm_NCORR_YJ %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(psvm_NCORR_YJ)

j_index_SE <- get_j_index_SE(psvm_NCORR_YJ)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "pSVM",
    Preprocessing = "NCORR_YJ",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble, tmp_tibble)

## ncorrorq_psvm
psvm_NCORRORQ <-models_svm_mlp$result[[12]]

best_perf_hyp <- psvm_NCORRORQ %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(psvm_NCORRORQ)

j_index_SE <- get_j_index_SE(psvm_NCORRORQ)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <- tibble(Model = "pSVM", Preprocessing= "NCORRORQ", mean_AUC = mean_auc, AUC_SE = AUC_SE, mean_j_index = mean_j_index,j_index_SE=j_index_SE, Correlation = corr_tune, n_predictors =  n_predictors)

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## pca_psvm 
psvm_pca <-models_svm_mlp$result[[17]]

best_perf_hyp <- psvm_pca %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(psvm_pca)

j_index_SE <- get_j_index_SE(psvm_pca)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_YeoJohnson(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>%
  step_pca(all_predictors(), num_comp = best_perf_hyp$num_comp[1]) %>%
  step_normalize(all_numeric_predictors()) %>% 
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "pSVM",
    Preprocessing = "PCA",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

final_tibble

# rbfSVM 

## ncorr_rbfsvm 
rbfsvm_NCORR <-models_svm_mlp$result[[3]]

best_perf_hyp <- rbfsvm_NCORR %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(rbfsvm_NCORR)

j_index_SE <- get_j_index_SE(rbfsvm_NCORR)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <-  recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>% 
  step_corr(all_predictors(), threshold = corr_tune) %>% prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <- tibble(Model = "rbfSVM", Preprocessing= "NCORR", mean_AUC = mean_auc, AUC_SE = AUC_SE, mean_j_index = mean_j_index,j_index_SE=j_index_SE, Correlation = corr_tune, n_predictors =  n_predictors)

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## ncorryj_rbfsvm

rbfsvm_NCORR_YJ <-models_svm_mlp$result[[8]]

best_perf_hyp <- rbfsvm_NCORR_YJ %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(rbfsvm_NCORR_YJ)

j_index_SE <- get_j_index_SE(rbfsvm_NCORR_YJ)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "rbfSVM",
    Preprocessing = "NCORR_YJ",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## ncorrorq_rbfsvm
rbfsvm_NCORRORQ <-models_svm_mlp$result[[13]]

best_perf_hyp <- rbfsvm_NCORRORQ %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(rbfsvm_NCORRORQ)

j_index_SE <- get_j_index_SE(rbfsvm_NCORRORQ)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "rbfSVM",
    Preprocessing = "NCORRORQ",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## pca_rbfsvm 
rbfsvm_pca <-models_svm_mlp$result[[18]]

best_perf_hyp <- rbfsvm_pca %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(rbfsvm_pca)

j_index_SE <- get_j_index_SE(rbfsvm_pca)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_YeoJohnson(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>%
  step_pca(all_predictors(), num_comp = best_perf_hyp$num_comp[1]) %>%
  step_normalize(all_numeric_predictors()) %>% 
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "rbfSVM",
    Preprocessing = "PCA",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

# MLP_nnet

## ncorr_mlp_nnet 

mlp_nnet_NCORR <-models_svm_mlp$result[[4]]

best_perf_hyp <- mlp_nnet_NCORR %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(mlp_nnet_NCORR)

j_index_SE <- get_j_index_SE(mlp_nnet_NCORR)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <-  recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>% 
  step_corr(all_predictors(), threshold = corr_tune) %>% prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "MLP_nnet",
    Preprocessing = "NCORR",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

final_tibble

## ncorryj_mlp_nnet
mlp_nnet_NCORR_YJ <-models_svm_mlp$result[[9]]

best_perf_hyp <- mlp_nnet_NCORR_YJ %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(mlp_nnet_NCORR_YJ)

j_index_SE <- get_j_index_SE(mlp_nnet_NCORR_YJ)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "MLP_nnet",
    Preprocessing = "NCORR_YJ",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## ncorrorq_mlp_nnet
mlp_nnet_NCORRORQ <-models_svm_mlp$result[[14]]

best_perf_hyp <- mlp_nnet_NCORRORQ %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(mlp_nnet_NCORRORQ)

j_index_SE <- get_j_index_SE(mlp_nnet_NCORRORQ)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "MLP_nnet",
    Preprocessing = "NCORRORQ",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## pca_mlp_nnet 
mlp_nnet_pca <-models_svm_mlp$result[[19]]

best_perf_hyp <- mlp_nnet_pca %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(mlp_nnet_pca)

j_index_SE <- get_j_index_SE(mlp_nnet_pca)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_YeoJohnson(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>%
  step_pca(all_predictors(), num_comp = best_perf_hyp$num_comp[1]) %>%
  step_normalize(all_numeric_predictors()) %>% 
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "MLP_nnet",
    Preprocessing = "PCA",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

# BAG_MLP

## ncorr_bag_mlp 
bag_mlp_NCORR <-models_svm_mlp$result[[5]]

best_perf_hyp <- bag_mlp_NCORR %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(bag_mlp_NCORR)

j_index_SE <- get_j_index_SE(bag_mlp_NCORR)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <-  recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>% 
  step_corr(all_predictors(), threshold = corr_tune) %>% prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "BAG_MLP",
    Preprocessing = "NCORR",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## ncorryj_bag_mlp
bag_mlp_NCORR_YJ <-models_svm_mlp$result[[10]]

best_perf_hyp <- bag_mlp_NCORR_YJ %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(bag_mlp_NCORR_YJ)

j_index_SE <- get_j_index_SE(bag_mlp_NCORR_YJ)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "BAG_MLP",
    Preprocessing = "NCORR_YJ",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## ncorrorq_bag_mlp
bag_mlp_NCORRORQ <-models_svm_mlp$result[[15]]

best_perf_hyp <- bag_mlp_NCORRORQ %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(bag_mlp_NCORRORQ)

j_index_SE <- get_j_index_SE(bag_mlp_NCORRORQ)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "BAG_MLP",
    Preprocessing = "NCORRORQ",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## pca_bag_mlp 
bag_mlp_pca <-models_svm_mlp$result[[20]]

best_perf_hyp <- bag_mlp_pca %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(bag_mlp_pca)

j_index_SE <- get_j_index_SE(bag_mlp_pca)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_YeoJohnson(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>%
  step_pca(all_predictors(), num_comp = best_perf_hyp$num_comp[1]) %>%
  step_normalize(all_numeric_predictors()) %>% 
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "BAG_MLP",
    Preprocessing = "PCA",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

# MARS
models_mars_knn <- readRDS("/mnt/data/andrei/Data/HeteroR/results/models/scheme12/models_mars_knn.rds")

## ncorr_mars 
mars_NCORR <-models_mars_knn$result[[1]]

best_perf_hyp <- mars_NCORR %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(mars_NCORR)

j_index_SE <- get_j_index_SE(mars_NCORR)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <-  recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>% 
  step_corr(all_predictors(), threshold = corr_tune) %>% prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "MARS",
    Preprocessing = "NCORR",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble, tmp_tibble)

## ncorryj_mars
mars_NCORR_YJ <-models_mars_knn$result[[4]]

best_perf_hyp <- mars_NCORR_YJ %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(mars_NCORR_YJ)

j_index_SE <- get_j_index_SE(mars_NCORR_YJ)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "MARS",
    Preprocessing = "NCORR_YJ",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble, tmp_tibble)

## ncorrorq_mars
mars_NCORRORQ <-models_mars_knn$result[[7]]

best_perf_hyp <- mars_NCORRORQ %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(mars_NCORRORQ)

j_index_SE <- get_j_index_SE(mars_NCORRORQ)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "MARS",
    Preprocessing = "NCORRORQ",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## pca_mars 
mars_pca <-models_mars_knn$result[[10]]

best_perf_hyp <- mars_pca %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(mars_pca)

j_index_SE <- get_j_index_SE(mars_pca)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_YeoJohnson(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>%
  step_pca(all_predictors(), num_comp = best_perf_hyp$num_comp[1]) %>%
  step_normalize(all_numeric_predictors()) %>% 
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "MARS",
    Preprocessing = "PCA",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

# BAG MARS

## ncorr_bag_mars 
bag_mars_NCORR <-models_mars_knn$result[[3]]

best_perf_hyp <- bag_mars_NCORR %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(bag_mars_NCORR)

j_index_SE <- get_j_index_SE(bag_mars_NCORR)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <-  recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>% 
  step_corr(all_predictors(), threshold = corr_tune) %>% prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "BAG_MARS",
    Preprocessing = "NCORR",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## ncorryj_bag_mars


bag_mars_NCORR_YJ <-models_mars_knn$result[[6]]

best_perf_hyp <- bag_mars_NCORR_YJ %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(bag_mars_NCORR_YJ)

j_index_SE <- get_j_index_SE(bag_mars_NCORR_YJ)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "BAG_MARS",
    Preprocessing = "NCORR_YJ",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## ncorrorq_mars
bag_mars_NCORRORQ <-models_mars_knn$result[[9]]

best_perf_hyp <- bag_mars_NCORRORQ %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(bag_mars_NCORRORQ)

j_index_SE <- get_j_index_SE(bag_mars_NCORRORQ)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "BAG_MARS",
    Preprocessing = "NCORRORQ",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble, tmp_tibble)

## pca_mars 
bag_mars_pca <-models_mars_knn$result[[12]]

best_perf_hyp <- bag_mars_pca %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(bag_mars_pca)

j_index_SE <- get_j_index_SE(bag_mars_pca)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_YeoJohnson(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>%
  step_pca(all_predictors(), num_comp = best_perf_hyp$num_comp[1]) %>%
  step_normalize(all_numeric_predictors()) %>% 
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "BAG_MARS",
    Preprocessing = "PCA",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

# KNN 

## ncorr_knn 
knn_NCORR <-models_mars_knn$result[[2]]

best_perf_hyp <- knn_NCORR %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(knn_NCORR)

j_index_SE <- get_j_index_SE(knn_NCORR)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <-  recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>% 
  step_corr(all_predictors(), threshold = corr_tune) %>% prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "KNN",
    Preprocessing = "NCORR",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## ncorryj_knn
knn_NCORR_YJ <-models_mars_knn$result[[5]]

best_perf_hyp <- knn_NCORR_YJ %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(knn_NCORR_YJ)

j_index_SE <- get_j_index_SE(knn_NCORR_YJ)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "KNN",
    Preprocessing = "NCORR_YJ",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## ncorrorq_knn
knn_NCORRORQ <-models_mars_knn$result[[8]]

best_perf_hyp <- knn_NCORRORQ %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(knn_NCORRORQ)

j_index_SE <- get_j_index_SE(knn_NCORRORQ)

corr_tune <-best_perf_hyp$corr_tune[1]

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
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "KNN",
    Preprocessing = "NCORRORQ",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## pca_knn 
knn_pca <-models_mars_knn$result[[11]]

best_perf_hyp <- knn_pca %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(knn_pca)

j_index_SE <- get_j_index_SE(knn_pca)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_YeoJohnson(all_numeric_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>%
  step_pca(all_predictors(), num_comp = best_perf_hyp$num_comp[1]) %>%
  step_normalize(all_numeric_predictors()) %>% 
  prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "KNN",
    Preprocessing = "PCA",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

# RF
models_rf_bt <- readRDS("/mnt/data/andrei/Data/HeteroR/results/models/scheme12/models_xgb_rf_light.rds")

## base_rf 
rf_base <-models_rf_bt$result[[1]]

best_perf_hyp <- rf_base %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(rf_base)

j_index_SE <- get_j_index_SE(rf_base)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <-   recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>% prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "RF",
    Preprocessing = "BASE",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## base_boruta_rf 
rf_base_boruta <-models_rf_bt$result[[3]]

best_perf_hyp <- rf_base_boruta %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(rf_base_boruta)

j_index_SE <- get_j_index_SE(rf_base_boruta)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <-    recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>%
  step_select_boruta(all_predictors(), outcome = "resistance") %>% prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "RF",
    Preprocessing = "BASE_BORUTA",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

# BT

## base_bt
bt_base <-models_rf_bt$result[[2]]

best_perf_hyp <- bt_base %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(bt_base)

j_index_SE <- get_j_index_SE(bt_base)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <-   recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>% prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "BT",
    Preprocessing = "BASE",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## base_bt_bres
base_bt_bres <- readRDS("/mnt/data/andrei/Data/HeteroR/results/models/scheme12/base_bt_bres_light.rds")

best_perf_hyp <- base_bt_bres %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-NA

j_index_SE <- NA

get_j_index <- function(res){res %>%   show_best("j_index",n =40) %>% filter(.config == (
    res %>% show_best("roc_auc", n = 1) %>% select(.config) %>%  as.character()
)) %>%  select(mean) %>% as.character()}

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <-   recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>% prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "BT_bres",
    Preprocessing = "BASE",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

## base_boruta_bt 
bt_base_boruta <-models_rf_bt$result[[4]]

best_perf_hyp <- bt_base_boruta %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <-get_j_index(bt_base_boruta)

j_index_SE <- get_j_index_SE(bt_base_boruta)

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <-    recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>%
  step_select_boruta(all_predictors(), outcome = "resistance") %>% prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()


tmp_tibble <-
  tibble(
    Model = "BT",
    Preprocessing = "BASE_BORUTA",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)


## base_boruta_bt_bres
bt_base_boruta_bres <- readRDS("/mnt/data/andrei/Data/HeteroR/results/models/scheme12/base_bt_boruta_bres_light.rds")

best_perf_hyp <- bt_base_boruta_bres %>% show_best(n=20 )

mean_auc <- best_perf_hyp$mean[1]

AUC_SE <- best_perf_hyp$std_err[1]

mean_j_index <- NA

j_index_SE <- NA

corr_tune <-best_perf_hyp$corr_tune[1]

n_predictors <-    recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>%
  step_select_boruta(all_predictors(), outcome = "resistance") %>% prep() %>% juice() %>% colnames()  %>% discard(~ .x %in% predictors_to_remove) %>% length()

tmp_tibble <-
  tibble(
    Model = "BT_bres",
    Preprocessing = "BASE_BORUTA",
    mean_AUC = mean_auc,
    AUC_SE = AUC_SE,
    mean_j_index = mean_j_index,
    j_index_SE = j_index_SE,
    Correlation = corr_tune,
    n_predictors =  n_predictors
  )

final_tibble <- bind_rows(final_tibble,tmp_tibble)

final_tibble

# run models and save them
library(optparse)

#### CLI parsing ####
option_list <- list(
  make_option(c("-m", "--model"),
              type = "character",
              default = NULL,
              help = "A model type (should be one of the following: 'lr', 'mars', 'bag_mars', 'lsvm', 'psvm', 'rf', 'knn', 'bt', 'nb', 'bag_mlp', 'mlp_keras', 'mlp_nnet')",
              metavar = "character"),
  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "output file (.rds)",
              metavar = "character"),
  make_option(c("-t", "--threads"),
              type = "integer",
              default = 2,
              help = "number of threads for resampling stage",
              metavar = "integer"),
  make_option(c("-r", "--recipe"),
              type = "character",
              default = NULL,
              help = "recipe to use (should be one of the following: 'main', 'ncorr', 'pca', 'umap', 'ncorr-orq')",
              metavar = "character"),
  make_option(c("-s", "--search"),
              type = "character",
              default = "space",
              help = "space-filling or Bayesian grid search (should be one of the following: 'bayes' or 'space')",
              metavar = "character"),
  make_option(c("-p", "--points"),
              type = "integer",
              default = 8,
              help = "number of intial random points in Bayesian tuning, a space-filling design will be used to populate a preliminary set of results. 
              For good results, the number of initial values should be more than the number of parameters being optimized",
              metavar = "integer"),
  make_option(c("-i", "--iterations"),
              type = "integer",
              help = "number of iterations in Bayesian tuning",
              default = 50,
              metavar = "integer"),
  make_option(c("-n", "--no_improve"),
              type = "integer",
              help = "number of iterations to stop after, if improvement",
              default = 30,
              metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# check validity of options
# if (opt$search != "bayes" | opt$search != "space"){
#   print("ERROR! Unknown type of grid search!")
#   quit(status = 1)
# } 


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
  update_role(resistance, new_role = "outcome") %>% 
  step_nzv(all_predictors()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_orderNorm(all_predictors()) %>% 
  step_normalize(all_predictors()) %>%
  step_pca(all_predictors(), num_comp = tune()) %>% 
  step_normalize(all_predictors()) %>% 
  step_smote(resistance, over_ratio = 1, seed = 100)

umap_recipe <- recipe(resistance ~., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_orderNorm(all_numeric_predictors()) %>% 
  step_normalize(all_predictors()) %>% 
  step_umap(all_numeric_predictors(), 
            num_comp = 24, 
            min_dist = 0.1, 
            neighbors = 15) %>%
  step_smote(resistance, over_ratio = 1, seed = 100)

ncorq_recipe <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>% 
  step_nzv(all_predictors()) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_orderNorm(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_corr(threshold = 0.75) %>%
  step_smote(resistance, over_ratio = 1, seed = 100)

#### FOLDS & METRICS ####
cv_folds <- vfold_cv(df_train, 
                     strata = "resistance", 
                     v = 10, 
                     repeats = 10) 

cls_metrics <- metric_set(roc_auc, j_index) # metrics for imbalanced classes

#### FUNCTIONS ####
set_model <- function(mod, cores) {
  # create model specification
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
  } else if(mod == "lsvm") {
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
  } else if (mod == "psvm") {
    my_mod <- svm_poly(
      cost = tune(),
      degree = tune(),
      scale_factor = tune(),
      margin = NULL ) %>% # regression only 
      set_mode("classification") %>%
      set_engine("kernlab", num.threads = cores)
  } else if (mod == "mlp_keras"){
    my_mod <-
      mlp(hidden_units = tune(), 
          penalty = NULL, 
          epochs = tune(),
          dropout = tune(),
          activation = tune(),
          learn_rate = NULL) %>%
      set_mode("classification") %>%
      set_engine("keras", num.threads = cores)
  } else if (mod == "mlp_nnet") {
    my_mod <-
      mlp(hidden_units = tune(), 
          penalty = tune(), 
          epochs = tune()) %>%
      set_mode("classification") %>%
      set_engine("nnet", num.threads = cores)
  } else if (mod == "nb") {
    library(discrim) # for NB with engine 'klaR'
    my_mod <- naive_Bayes(
      mode = "classification",
      smoothness = tune(),
      Laplace = tune(),
      engine = "klaR"
    )
  } else if (mod == "bag_mars") {
    library(baguette) # for bag_mars
    my_mod <- bag_mars(
      mode = "classification",
      num_terms = tune(),
      prod_degree = tune(),
      prune_method = tune(),
      engine = "earth"
    )
  } else if (mod == "bag_mlp") {
    my_mod <- bag_mlp(
      mode = "classification",
      hidden_units = tune(),
      penalty = tune(),
      epochs = tune(),
      engine = "nnet"
    )
  }
  return(my_mod)
}


set_wf <- function(mod, rec, cores){
  # create a workflow
  # mod: model type, one of: lr, knn, mars, svm, rf, bt, nnet
  # rec: recipe object (one of: main, ncorr, pca, umap)
  # rec must be in GlobalEnv
  
  if (rec == "main"){
    rc <- main_recipe
  } else if (rec == "ncorr") {
    rc <- ncorr_recipe
  } else if (rec == "pca") {
    rc <- pca_recipe
  } else if (rec == "umap") {
    rc <- umap_recipe
  } else if (rec == "ncorr-orq") {
    rc <- ncorq_recipe
  } else {
    print("ERROR! Undefined recipe!")
    quit(status = 1)
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
if (opt$search == "space"){
  model_res <- my_wf %>%
          tune_grid(
              grid = opt$points,
              resamples = cv_folds,
              control = control_grid(save_pred = TRUE, 
                                     save_workflow = TRUE),
              metrics = cls_metrics)
} else if (opt$search == "bayes"){
  model_res <- my_wf %>% 
  tune_bayes(
    resamples = cv_folds,
    # To use non-default parameter ranges
    param_info = param_set,
    # Generate N at semi-random to start
    initial = opt$points,
    iter = opt$iterations,
    # How to measure performance?
    metrics = metric_set(roc_auc),
    control = control_bayes(no_improve = opt$no_improve, 
                            verbose = FALSE, 
                            save_pred = TRUE, 
                            save_workflow = TRUE))
}

# show performance
perf <- model_res %>% show_best("roc_auc", n = 5)
print(perf)

#### SAVE MODEL ####
saveRDS(object = model_res, file = opt$output)


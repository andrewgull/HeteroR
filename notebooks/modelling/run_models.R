# run models and save them
library(optparse)

#### CLI parsing ####
option_list <- list(
  make_option(
    c("-m", "--model"),
    type = "character",
    default = NULL,
    help = "A model type (should be one of the following: 'lr', 'mars', 'bag_mars', 'lsvm', 'psvm', 'rbfsvm', 'rf', 'knn', 'bt', 'nb', 'bag_mlp', 'mlp_keras', 'mlp_nnet')",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = NULL,
    help = "output file (.rds)",
    metavar = "character"
  ),
  make_option(
    c("-t", "--threads"),
    type = "integer",
    default = 2,
    help = "number of threads",
    metavar = "2"
  ),
  make_option(
    c("-r", "--recipe"),
    type = "character",
    default = NULL,
    help = "recipe to use (should be one of the following: 'main', 'ncorr', 'pca', 'ncorq', 'pcayj', 'ncoryj')",
    metavar = "character"
  ),
  make_option(
    c("-g", "--grid_search"),
    type = "character",
    default = "space",
    help = "space-filling, Bayesian, racing or racing + Bayesian grid search (should be one of the following: 'bayes', 'space', 'race', 'race+bayes')",
    metavar = "space"
  ),
  make_option(
    c("-p", "--points"),
    type = "integer",
    default = 8,
    help = "number of intial random points in Bayesian grid search, the space-filling design will be used to populate a preliminary set of results.
              for good results, the number of initial values should be more than the number of parameters being optimized",
    metavar = "8"
  ),
  make_option(
    c("-i", "--iterations"),
    type = "integer",
    help = "number of iterations in Bayesian grid search",
    default = 50,
    metavar = "50"
  ),
  make_option(
    c("-n", "--no_improve"),
    type = "integer",
    help = "number of iterations to stop after, if no improvement",
    default = 30,
    metavar = "30"
  ),
  make_option(
    c("-f", "--folds"),
    type = "integer",
    help = "number of folds in cross-validation",
    default = 10,
    metavar = "10"
  ),
  make_option(
    c("-a", "--proportion"),
    type = "double",
    help = "stratified split proportion",
    default = 0.8,
    metavar = "0.8"
  ),
  make_option(
    c("-s", "--resamples"),
    type = "integer",
    help = "number of resamples",
    default = 10,
    metavar = "10"
  ),
  make_option(
    c("-c", "--classification"),
    type = "character",
    help = "which classification scheme to choose? (one of: 12, 13, 123)",
    default = "12",
    metavar = "12"
  ),
  make_option(
    c("-z", "--corr_threshold"),
    type = "double",
    help = "correlation threshold for NOCORR recipe (currently is icluded in tuning)",
    default = 0.75,
    metavar = "0.75"
  ),
  make_option(
    c("-b", "--burnin"),
    type = "integer",
    help = "Number of initial resamples to use in racing grid search",
    default = 10,
    metavar = "10"
  ),
  make_option(
    c("-d", "--rds"),
    type = "character",
    help = "RDS file with resamples object",
    default = NULL,
    metavar = "file.rds"
  )
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
library(finetune) # for racing grid search
library(themis) # for SMOTE
library(bestNormalize) # for ORQ-norm
library(embed) # for UMAP

if (opt$model == "nb"){
  library(discrim) # for NB with engine 'klaR'
} else if (opt$model == "bag_mars" | opt$model == "bag_mlp"){
  library(baguette) # for bagged models
}


#### DATA ####
path_data <- "/home/andrei/GitProjects/HeteroR/notebooks/modelling/data/features_strain.csv"

if (opt$classification == "123") {
  path_labels <- "/home/andrei/GitProjects/HeteroR/notebooks/modelling/data/heteroresistance_testing_gr123.csv"
} else if (opt$classification == "12") {
  path_labels <- "/home/andrei/GitProjects/HeteroR/notebooks/modelling/data/heteroresistance_testing_gr12.csv"
} else if (opt$classification == "13") {
  path_labels <- "/home/andrei/GitProjects/HeteroR/notebooks/modelling/data/heteroresistance_testing_gr13.csv"
}

data_strain <- readr::read_csv(path_data,
                               na = c("NA", "-Inf"),
                               show_col_types = FALSE)

hr_testing <- readr::read_csv(path_labels,
                              show_col_types = FALSE)

data_strain <- data_strain %>%
  left_join(hr_testing, by = "strain")

data_strain <- data_strain %>%
  mutate(n.beta.lac.4 = factor(ifelse(n.beta.lac > 4, "yes", "no"))) %>%
  relocate(n.beta.lac.4, .before = "n.plasmids") %>%
  filter(resistance != "R", strain != "DA63310") %>%
  mutate(
    resistance = factor(resistance, levels = c("HR", "nonHR")),
    chrom.status = factor(chrom.status)
  )

data_strain$N50 <- NULL
data_strain[is.na(data_strain)] <- 0

#### SPLIT ####
# same seed number as in modelling.Rmd
set.seed(124)

# splitting proportion should be the same
data_split <- initial_split(data_strain, 
                            prop = opt$proportion, 
                            strata = resistance)

df_train <- training(data_split)
df_test <- testing(data_split)

#### PREPROCESSING RECIPES ####
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
  #step_corr(threshold = opt$corr_threshold) %>%
  step_corr(threshold = tune("corr_tune")) %>% 
  step_smote(resistance, over_ratio = 1, seed = 100)

pca_recipe <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  update_role(resistance, new_role = "outcome") %>% 
  step_nzv(all_predictors()) %>% 
  #step_normalize(all_numeric_predictors()) %>%
  step_orderNorm(all_numeric_predictors()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_smote(resistance, over_ratio = 1, seed = 100) %>% 
  step_pca(all_predictors(), threshold = .9)

pcayj_recipe <- 
  recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_YeoJohnson(all_numeric_predictors()) %>% 
  step_normalize(all_numeric_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>% 
  step_pca(all_predictors(), num_comp = tune()) %>% 
  step_normalize(all_numeric_predictors())

yj_recipe <- 
  recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>%
  step_nzv(all_predictors()) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_YeoJohnson(all_numeric_predictors()) %>% 
  step_normalize(all_numeric_predictors()) %>%
  step_smote(resistance, over_ratio = 1, seed = 100) %>% 
  step_corr(threshold = tune("corr_tune"))

ncorq_recipe <- recipe(resistance ~ ., data = df_train) %>%
  update_role(strain, new_role = "ID") %>% 
  step_nzv(all_predictors()) %>%
  step_orderNorm(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_corr(threshold = tune("corr_tune")) %>%
  step_smote(resistance, over_ratio = 1, seed = 100)

#### FOLDS & METRICS ####
cv_folds <- vfold_cv(df_train, 
                     strata = "resistance", 
                     v = opt$folds, 
                     repeats = opt$resamples) 
# metrics for imbalanced classes
imbalanced_metrics <- metric_set(roc_auc, j_index) 

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
  } else if (mod == "rbfsvm"){
    my_mod <- svm_rbf(
      cost = tune(),
      rbf_sigma = tune(),
      margin = tune(),
      mode = "classification"
    ) %>% 
      set_engine("kernlab", num.threads = cores)
  } else if (mod == "mlp_keras"){
    my_mod <-
      mlp(hidden_units = tune(), 
          penalty = NULL, 
          epochs = tune(),
          dropout = tune(),
          activation = "relu",
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
    my_mod <- naive_Bayes(
      mode = "classification",
      smoothness = tune(),
      Laplace = tune(),
      engine = "klaR"
    )
  } else if (mod == "bag_mars") {
    my_mod <- bag_mars(
      mode = "classification",
      num_terms = tune(),
      prod_degree = tune(),
      prune_method = tune(),
      engine = "earth"
    )
  } else if (mod == "bag_mlp") {
    my_mod <- bag_mlp(
      hidden_units = tune(),
      penalty = tune(),
      epochs = tune()) %>% 
      set_engine("nnet", num.threads = cores) %>% 
      set_mode("classification") %>% 
      translate()
  }
  return(my_mod)
}


set_rec <- function(rec, cores){
  
  # choose a recipe
  if (rec == "main"){
    rc <- main_recipe
  } else if (rec == "ncorr") {
    rc <- ncorr_recipe
  } else if (rec == "pca") {
    rc <- pca_recipe
  } else if (rec == "ncorq") {
    rc <- ncorq_recipe
  } else if (rec == "pcayj") {
    rc <- pcayj_recipe
  } else if (rec == "ncoryj") {
    rc <- yj_recipe
  } else {
    print("ERROR! Undefined recipe!")
    quit(status = 1)
  }
  return(rc)
}

#### CREATE A WORKFLOW ####
# using chosen model specification 
# and chosen recipe
my_wf <- workflow() %>% 
  add_model(set_model(mod = opt$model, cores = opt$threads)) %>% 
  add_recipe(set_rec(opt$recipe))

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
if (opt$grid_search == "space"){
  model_res <- my_wf %>%
    tune_grid(
      param_info = param_set,
      grid = opt$points,
      resamples = cv_folds,
      control = control_grid(save_pred = TRUE,
                             save_workflow = TRUE),
      metrics = imbalanced_metrics
    )
} else if (opt$grid_search == "bayes"){
  model_res <- my_wf %>%
    tune_bayes(
      resamples = cv_folds,
      # To use non-default parameter ranges
      param_info = param_set,
      initial = opt$points,
      iter = opt$iterations,
      metrics = metric_set(roc_auc),
      control = control_bayes(
        no_improve = opt$no_improve,
        verbose = FALSE,
        save_pred = TRUE,
        save_workflow = TRUE
      )
    )
} else if (opt$grid_search == "race") {
  model_res <- my_wf %>%
    tune_race_win_loss(
      param_info = param_set,
      resamples = cv_folds,
      grid = opt$points,
      metrics = imbalanced_metrics,
      control = control_race(
        verbose_elim = TRUE,
        save_pred = TRUE,
        save_workflow = TRUE,
        burn_in = 10
      )
    )
} else if (opt$grid_search == "race+bayes") {
  # explore many points and optimize the winning ones
  resamp_obj <- readRDS(opt$rds) # resamples to optimize
  
  model_res <- my_wf %>%
    tune_bayes(
      resamples = cv_folds,
      # To use non-default parameter ranges
      param_info = param_set,
      # Generate N at semi-random to start
      initial = resamp_obj,
      iter = opt$iterations,
      # How to measure performance?
      metrics = metric_set(roc_auc),
      control = control_bayes(
        no_improve = opt$no_improve,
        verbose = FALSE,
        save_pred = TRUE,
        save_workflow = TRUE
      )
    )
}

# show performance
perf <- model_res %>% show_best("roc_auc", n = 5)
print(perf)

#### SAVE MODEL RESAMPLES ####
saveRDS(object = model_res, file = opt$output)


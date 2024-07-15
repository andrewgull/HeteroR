## MAKE ROC OBJECT ##
# for autoplot()
make_roc <- function(mod_res, title, filtering_metric="roc_auc"){
  
  mod_best <- mod_res %>%
    select_best(metric = filtering_metric)
  
  mod_roc <- mod_res %>% 
    collect_predictions(parameters = mod_best) %>% 
    roc_curve(resistance, .pred_HR) %>% 
    mutate(model = title)
  
  return(mod_roc)
}

make_pr <- function(mod_res, title, filtering_metric="roc_auc"){
  
  mod_best <- mod_res %>%
    select_best(metric = filtering_metric)
  
  mod_pr <- mod_res %>% 
    collect_predictions(parameters = mod_best) %>% 
    pr_curve(resistance, .pred_HR) %>% 
    mutate(model = title)
  
  return(mod_pr)
}

# Also, there is a script 'run_models.R' to run resampling from terminal
# Bayesian ANOVA
# function to get a name of the best preprocessor from your resample
best_preprocessor <- function(res_obj) {
  res_obj %>% select_best("roc_auc") %>% pull(.config)
}

# function to make a df with AUCs from the best preproc
best_aucs <- function(res_obj, model_name){
  res_obj %>%
    collect_metrics(summarize = FALSE) %>%
    filter(.metric == "roc_auc", .config == best_preprocessor(res_obj)) %>%
    select(id, id2, {{model_name}} := .estimate)
}


# returns df of best re-sample from each model submitted
# use it for perf_mod()
res_comp_table <- function(res_list, mod_names){
  # res_list: list of resamples objects
  # mod_names: (model) names that become column names
  comp_df <-
    map2(res_list, mod_names,
         ~ best_aucs(res_obj = .x, model_name = .y)) %>% 
    reduce(inner_join, by = c("id", "id2")) %>% 
    set_names("id", "id2", mod_names) %>% 
    unite(id, id, id2)
  return(comp_df)
}

comp_plots <- function(col_name, orig_data, prep_data){
  p_orig <- 
    orig_data %>% 
    ggplot(aes_string(x = col_name)) + 
    geom_histogram(bins = 30, color = "white", fill = "blue", alpha = 1/3) + 
    ggtitle(paste0("Original data (", col_name, ")"))
  
  p_norm <- 
    prep_data %>% 
    ggplot(aes_string(x = col_name)) + 
    geom_histogram(bins = 30, color = "white", fill = "red", alpha = 1/3) + 
    ggtitle(paste0("Processed data (", col_name, ")"))
  
  return(list(p_orig, p_norm))
}

plot_validation_results <- function(dat, components, colors_vec) {
  # dat - data table
  # components - vector of components names to plot
  dat %>%
    select(-strain) %>%
    select(all_of(c(components, "resistance"))) %>% 
    # Create the scatterplot matrix
    ggplot(aes(x = .panel_x, y = .panel_y, color = resistance, fill = resistance)) +
    geom_point(alpha = 0.4, size = 1) +
    ggforce::geom_autodensity(alpha = .3) +
    ggforce::facet_matrix(vars(-resistance), layer.diag = 2) + 
    scale_color_manual(values = colors_vec) +
    scale_fill_manual(values = colors_vec) +
    theme_minimal() +
    theme(legend.title = element_blank())
}

# to read workflowsets file and return the required workflow by ID
read_wfset <- function(wf_file, wf_id, models_path="/mnt/data/andrei/Data/HeteroR/results/models/scheme12/"){
  readRDS(paste0(models_path, wf_file)) %>% 
    filter(wflow_id == wf_id) %>% 
    pull(result) %>% 
    pluck(1)
} 

maxj <- function(threshold_df) {
  # find max j-index
  max_j_index_threshold <- threshold_df %>%
    filter(.metric == "j_index") %>%
    filter(.estimate == max(.estimate)) %>%
    pull(.threshold)
  
  # max_j_index_threshold may be a vector, use its last element
  if (length(max_j_index_threshold) > 1) {
    max_j_index_threshold <-
      max_j_index_threshold[length(max_j_index_threshold)]
  }
  
  return(max_j_index_threshold)
}

get_threshold <- function(fit_obj){
  # collect sens, spec, j-index at various cut-offs
  threshold_data <-
    fit_obj %>%
    collect_predictions() %>%
    threshold_perf(resistance, .pred_HR, thresholds = seq(0.0, 1, by = 0.01)) %>%
    filter(.metric != "distance") %>%
    mutate(group = case_when(.metric == "sens" |
                               .metric == "spec" ~ "1",
                             TRUE ~ "2"))
  return(threshold_data)
}

senspec_plot <- function(fit_obj, title="") {
  # collect sens, spec, j-index at various cut-offs
  threshold_data <- get_threshold(fit_obj)
  
  # find max j-index
  max_j_index_threshold <- maxj(threshold_data)
  
  # plot metrics v cut-offs
  sens_spec_j_plot <-
    ggplot(threshold_data,
           aes(
             x = .threshold,
             y = .estimate,
             color = .metric,
             alpha = group
           )) +
    geom_line(size = 1) +
    #theme_minimal() +
    #scale_color_viridis_d(end = 0.9) +
    scale_color_brewer(palette = "Set2", guide = guide_legend(title = NULL) ) +
    scale_alpha_manual(values = c(.9, 1), guide = "none") +
    geom_vline(
      xintercept = max_j_index_threshold,
      alpha = .8,
      color = "grey30",
      linetype = "longdash"
    ) +
    labs(x = "Probability",
         y = "Metric Estimate",
         title = title) +
    theme_minimal()
  return(sens_spec_j_plot)
}

# function to pick the best config name accordong to ROC AUC
# so that you can use it to get its J-index
best_config_name <- function(resamples) {
  resamples %>% 
    show_best("roc_auc") %>% 
    pull(.config) %>% 
    pluck(1)
}

# use the function above to extract the best j-index
# for the same config that gives you best ROC
show_best_j <- function(resamples, n) {
  resamples %>% 
  show_best(metric = "j_index", n = n) %>% 
    filter(.config == best_config_name(resamples))
}

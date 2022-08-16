library(tidyverse)
library(caret)
library(optparse)
library(parallel)
library(doParallel)


### CLI parsing ###
library(optparse)

# CLI parsing
option_list <- list(
  make_option(c("-i", "--input"),
              type = "character",
              default = NULL,
              help = "A path to amp features ile",
              metavar = "character"),
  make_option(c("-o", "--output"),
              type = "character",
              default = NULL,
              help = "A path to an output RData file",
              metavar = "character"),
  make_option(c("-t", "--threads"),
              type = "integer",
              default = NULL,
              help = "number of threads",
              metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("An input must be provided", call. = FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("An output must be provided", call. = FALSE)
}
if (is.null(opt$threads)){
  print_help(opt_parser)
  stop("A number of threads must be provided", call. = FALSE)
}


### get paths ###
input_file <- opt$input
output_file <- opt$output
thread <- opt$threads

#################################
### READ AND PROCESS THE DATA ###
#################################

feat_amp <- read_csv(input_file, na = "NA")
feat_amp <- rename(feat_amp, n.plasmids="n_plasmids") %>% relocate(n.plasmids, .before="n.rep.total")


feat_amp$resistance <- as.factor(feat_amp$resistance)

feat_amp$n.beta.lac.3 <- ifelse(feat_amp$n.beta.lac > 3, 1, 0)
feat_amp <- relocate(feat_amp, n.beta.lac.3, .before="n.plasmids")

feat_amp <- feat_amp %>% select(-c(strain, AB))

##################
### IMPUTATION ###
##################

feat_amp$n.rep.tot.cen <- ifelse(is.na(feat_amp$n.rep.tot.cen), 0, feat_amp$n.rep.tot.cen)

feat_amp$n.rep.tot.non.cen <- ifelse(is.na(feat_amp$n.rep.tot.non.cen), 0, feat_amp$n.rep.tot.non.cen)

feat_amp$med.tot.rep.len <- ifelse(is.na(feat_amp$med.tot.rep.len), 0, feat_amp$med.tot.rep.len)

feat_amp$med.AR.len.cen <- ifelse(is.na(feat_amp$med.AR.len.cen), 0, feat_amp$med.AR.len.cen)

feat_amp$ampC.med.tot.rep.len <- ifelse(is.na(feat_amp$ampC.med.tot.rep.len), 0, feat_amp$ampC.med.tot.rep.len)

#############
### SCALE ###
#############

feat_amp_scaled <- as_tibble(scale(feat_amp[,c(2,4:39)]))

feat_amp_scaled$resistance <- feat_amp$resistance
feat_amp_scaled$n.beta.lac.3 <- feat_amp$n.beta.lac.3

feat_amp_scaled <- relocate(feat_amp_scaled, n.beta.lac.3, .before="n.beta.lac")
feat_amp_scaled <- relocate(feat_amp_scaled, resistance, .before="n.beta.lac.3")

######################
### SPLIT THE DATA ###
######################

set.seed(100)

in_train_all <- createDataPartition(y = feat_amp_scaled$resistance, p=0.75, list=FALSE)


training_all <- feat_amp_scaled[in_train_all,]
testing_all <- feat_amp_scaled[-in_train_all,]

x_all <- as.matrix(training_all[,-1]) 
y_all <- training_all$resistance

#######################
### CONFIG PARALLEL ###
#######################

cluster <- makeCluster(detectCores() - 8) # leave two cores
registerDoParallel(cluster)

###################
### TUNE PARAMS ###
###################
fit_ctrl_xgb_roc <- trainControl(method = "repeatedcv", 
                                 number = 5, repeats = 5, 
                                 allowParallel = T, 
                                 classProbs = T, 
                                 summaryFunction = twoClassSummary)
#############
### TRAIN ###
#############

fit_XGB <- train(x = x_all, y = y_all, 
                 method = "xgbDART", 
                 trControl = fit_ctrl_xgb_roc,
                 verbose = FALSE,
                 metric="ROC",
                 verbosity=0)
############
### TEST ###
############

predClasses <- predict(fit_XGB, newdata=testing_all)

confusionMatrix(data = predClasses, 
                reference = testing_all$resistance,
                mode="everything",
                positive="R")
############
### SAVE ###
############
save.image(f=output_file)





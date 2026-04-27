#### Random forest feature extraction
library(ranger)
library(tidyverse)
library(vip)
library(caret)
library(abc)
rm(list = ls())
source("/home/jorge/Rscripts/AuxAbc.R")

chunk <- "CORE"

names_statistics <- c("D_CHB_DEN","fD_CHB_DEN","D_EUR_DEN","fD_EUR_DEN",
                      "D_CHB_NEA","fD_CHB_NEA","D_EUR_NEA","fD_EUR_NEA",
                      "R_d_CHB_D","R_d_EUR_D","R_d_CHB_N","R_d_EUR_N",
                      "U_CHB_DEN","U_EUR_DEN","U_CHB_NEA","U_EUR_NEA",
                      "S_CHB","S_EUR","S_EAF","S_DEN","S_NEA",
                      "tw_CHB","tw_EUR","tw_EAF",
                      "tp_CHB","tp_EUR","tp_EAF",
                      "H1_CHB","H2_CHB","H12_CHB","H_CHB",
                      "H1_EUR","H2_EUR","H12_EUR","H_EUR",
                      "H1_EAF","H2_EAF","H12_EAF","H_EAF",
                      "fr_allel_CHB","fr_allel_EUR","fr_allel_EAF",
                      "fr_deni_old","fr_charg","fr_altai","fr_deni","fr_vin",
                      "fr_al_OOA","fr_al_nea",
                      "ZnS_CHB","ZnS_EUR","ZnS_EAF","ZnS_DEN","ZnS_NEA", 
                      "replicate")

dir <- "/home/jorge/PROJECTS.JORGE/SIMULATIONS"
out_dir <- paste0(dir, "/RESULTS.NEG/RF_output/")
train_path <- "/home/jorge/PROJECTS.JORGE/SLENDR/WORK_STATS/"

statistics_Deni_f <- preprocess(paste0(train_path,"DEN.NEG.TRAIN.stats.tsv"), "DENI")
statistics_Nea_f <- preprocess(paste0(train_path,"NEA.NEG.TRAIN.stats.tsv"), "NEA")
statistics_Anc_f <- preprocess(paste0(train_path,"ANC.NEG.TRAIN.stats.tsv"), "ANC")
statistics_D2n_f <- preprocess(paste0(train_path,"D2N.NEG.TRAIN.stats.tsv"), "D2N")

#### SPLIT DATA
den_indx <- sample(1:dim(statistics_Deni_f)[1],dim(statistics_Deni_f)[1]/3)
rep_den <- statistics_Deni_f[den_indx,]
lgl_den <- 1:dim(statistics_Deni_f)[1] %in% den_indx
tra_den <- statistics_Deni_f[!lgl_den,]

nea_indx <- sample(1:dim(statistics_Nea_f)[1],round(dim(statistics_Nea_f)[1]/3))
rep_nea <- statistics_Nea_f[nea_indx,]
lgl_nea <- 1:dim(statistics_Nea_f)[1] %in% nea_indx
tra_nea <- statistics_Nea_f[!lgl_nea,]

anc_indx <- sample(1:dim(statistics_Anc_f)[1],dim(statistics_Anc_f)[1]/3)
rep_anc <- statistics_Anc_f[anc_indx,]
lgl_anc <- 1:dim(statistics_Anc_f)[1] %in% anc_indx
tra_anc <- statistics_Anc_f[!lgl_anc,]

d2n_indx <- sample(1:dim(statistics_D2n_f)[1],dim(statistics_D2n_f)[1]/3)
rep_d2n <- statistics_D2n_f[d2n_indx,]
lgl_D2n <- 1:dim(statistics_D2n_f)[1] %in% d2n_indx
tra_D2n <- statistics_D2n_f[!lgl_D2n,]


training_table <- rbind(tra_den,tra_nea, tra_anc,tra_D2n) 
replicate_table <- rbind(rep_den, rep_nea, rep_anc, rep_d2n)
vector_selection <- c("S_DEN","fr_deni_old","ZnS_DEN","replicate","ZnS_NEA", "fr_charg","fr_altai","fr_vin")

# Model selection with training data
model_selection <- training_table %>% dplyr::select(-any_of(vector_selection))
model_selection$model <- as.factor(model_selection$model)

n_features <- 47

# Perm
Tree_Results_perm <- ranger(model ~., data = model_selection, write.forest = T, 
                           num.threads = 16, oob.error = T, keep.inbag = T, importance = "permutation")
default_rmse_perm <- sqrt(Tree_Results_perm$prediction.error)


hyper_grid_perm <- expand.grid(
  mtry = floor(n_features * c(.05, .25, .4)),
  min.node.size = c(1, 3, 5, 10), 
  replace = c(TRUE, FALSE),                               
  sample.fraction = c(.5, .63, .8),
  ntrees = c(500,1500,3000),
  rmse = NA                                               
)

for(i in seq_len(nrow(hyper_grid_perm))) {
  # fit model for ith hyperparameter combination
  fit <- ranger(
    formula         = model ~ ., 
    data            = model_selection, 
    num.trees       = hyper_grid_perm$ntrees[i],
    mtry            = hyper_grid_perm$mtry[i],
    min.node.size   = hyper_grid_perm$min.node.size[i],
    replace         = hyper_grid_perm$replace[i],
    sample.fraction = hyper_grid_perm$sample.fraction[i],
    verbose         = FALSE,
    num.threads = 16,
    respect.unordered.factors = 'order',
    importance = "permutation"
  )
  # export OOB error 
  hyper_grid_perm$rmse[i] <- sqrt(fit$prediction.error)
}

capture.output(hyper_grid_perm %>%
                 arrange(rmse) %>%
                 mutate(perc_gain = (default_rmse_perm - rmse) / default_rmse_perm * 100) %>%
                 head(10), file = paste0(out_dir, "Best_Params_perm.txt"))

#### BUILD RF with best grid parameters from last step and apply them to predict
rf_permutation <- ranger(model ~., data = model_selection, write.forest = T, num.trees = 3000, 
                      replace = F, mtry = 18, min.node.size = 5, sample.fraction = 0.8,
                      num.threads = 16, oob.error = T, keep.inbag = T, importance = "permutation", 
                      respect.unordered.factors = 'order')

rf_permutation$prediction.error

RF_training <- vip::vip(rf_permutation, num_features = 20, bar = FALSE) + ggtitle("Permutation RF")
ggsave(filename = paste0("RF_training.",chunk,".png"), path = out_dir, plot = RF_training, device = "png", units = "px", height = 3000, width = 3200)

cfmt_per <- confusionMatrix(data =rf_permutation$confusion.matrix ) #factor(cv$estim$tol0.1), factor(cv$true))
capture.output(cfmt_per, file = paste0(out_dir, "CFM_accuracy_permutation.txt"))

#### Replicate model
test_table <- rbind(rep_den,rep_nea, rep_d2n, rep_anc) %>% dplyr::select(-any_of(vector_selection))
rep_table_labels <- test_table$model
test_table$model <- as.factor(test_table$model)

### Same but with CV
most_important_permutation <- names(sort(rf_permutation$variable.importance, decreasing = T))

make_comparison_number_ft <- function(number) {
  cv_rf_per <- cv4postpr(rep_table_labels, test_table %>% dplyr::select(all_of(most_important_permutation[1:number])), 
                         nval=15, tol=.1, method = "mnlogistic", corr = T)
  s_rf_per <- summary(cv_rf_per)
  cfm_rf_per <- confusionMatrix(factor(cv_rf_per$estim$tol0.1), factor(cv_rf_per$true))
  capture.output(cfm_rf_per, file = paste0(out_dir, "ABC.CV.cfm.per.", number,".txt"))
  
}

make_comparison_number_ft(3)
make_comparison_number_ft(5)
make_comparison_number_ft(7) # This was the one with highest accuracy
make_comparison_number_ft(9)
make_comparison_number_ft(11)

### Comparison with NeuralNet method
cv_nnet <- cv4postpr(rep_table_labels, test_table %>% dplyr::select(-any_of(vector_selection),-model), 
                     nval=15, tol=.1, method="neuralnet")
s <- summary(cv_nnet, digit = 2)

cfm.10 <- confusionMatrix(factor(cv_nnet$estim$tol0.1), factor(cv_nnet$true))

capture.output(cfm.10, file = paste0(out_dir, "ABC.CV.cfm.nnet.0.10.txt"))

save.image("ABC_RF.4models.NEG.Rdata")


##### Regression
param_names <- c("sel_den","sel_nea","sel_chb","sel_ooa","sel_eur","onset_time","replicate","mutation","time_simulation")
parameters <- read_tsv("~/PROJECTS.JORGE/TREE_SEQUENCES/TRAINING_SIMULATIONS/D2N.MORE.NEG/twodnSimulation.Q5.parameters.tsv",
                       col_names = param_names, col_types = c("nnnnnnnn"))

make_regression <- function(stats_df, parameters, param_name, method_name = "permutation") {
  reg_den <- merge(stats_df, parameters, by = "replicate", sort = F) %>% 
  dplyr::select(-any_of(vector_selection), -any_of(param_names[!param_names == param_name ]),-model)
  
  train_reg <- reg_den[!lgl_D2n,]
  test_reg <- reg_den[lgl_D2n,]
  
  # Very usefull to use formulas inside functions
  formule <- as.formula(paste(param_name, "~ ."))
  
  # Perform the RF
  Tree_Results_reg <- ranger(formule , data = train_reg, write.forest = F, 
                             num.threads = 16, oob.error = T, keep.inbag = T, importance = method_name)
  default_rmse <- sqrt(Tree_Results_reg$prediction.error)
  
  hyper_grid <- expand.grid(
    mtry = floor(n_features * c(.05, .25, .4)),
    min.node.size = c(1, 5, 10), 
    ntrees = c(500, 1500, 3000),
    replace = c(TRUE, FALSE),                               
    sample.fraction = c(.5, .63, .8),                       
    rmse = NA                                               
  )
  
  for(i in seq_len(nrow(hyper_grid))) {
    # fit model for ith hyperparameter combination
    fit <- ranger(
      formula         = formule, 
      data            = train_reg, 
      num.trees       = hyper_grid$ntrees[i],
      mtry            = hyper_grid$mtry[i],
      min.node.size   = hyper_grid$min.node.size[i],
      replace         = hyper_grid$replace[i],
      sample.fraction = hyper_grid$sample.fraction[i],
      verbose         = FALSE,
      num.threads = 16,
      respect.unordered.factors = 'order',
      importance = method_name
    )
    # export OOB error 
    hyper_grid$rmse[i] <- sqrt(fit$prediction.error)
  }
  
  Best_hyperparams <- hyper_grid %>%
                   arrange(rmse) %>%
                   mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
                   head(10)

  capture.output(Best_hyperparams, file = paste0(out_dir, "Best_Params_", param_name, "_", method_name, ".txt"))
  
  # Use the best hyperparameters to perform accuracy tests
  Tree_Results_reg <- ranger(formule , data = train_reg, write.forest = T, 
                             mtry = Best_hyperparams$mtry[1], min.node.size = Best_hyperparams$min.node.size[1], 
                             sample.fraction = Best_hyperparams$sample.fraction[1],
                             num.trees = Best_hyperparams$ntrees[1],
                             replace = Best_hyperparams$replace[1], respect.unordered.factors = "order",
                             num.threads = 16, oob.error = T, keep.inbag = T, importance = method_name)
  
  best_variables <- names(sort(Tree_Results_reg$variable.importance, decreasing = T))
  capture.output(best_variables, file = paste0(out_dir, "Best_Variables_", param_name, "_", method_name, ".txt"))
  importance_plot <- vip::vip(Tree_Results_reg, num_features = 20, bar = FALSE) + ggtitle(param_name)
  ggsave(filename = paste0("RF_training.", param_name, ".",method_name,".png"), path = out_dir, plot = importance_plot, device = "png", units = "px", height = 3000, width = 3200)
  
  results_abc_reg <- data.frame("Mean_posterior" = rep(0, dim(test_reg)[1]), 
                                "True_param" = rep(0,dim(test_reg)[1]))
  
  # Perform the abc for each parameter
  for (i in 1:dim(test_reg)[1]) {
    abc_reg <- abc(target = test_reg[i,best_variables[1:5]], param = train_reg[,48], sumstat = train_reg[,best_variables[1:5]], method = "loclinear",tol = 0.1)
    results_abc_reg[i,1] <- mean(abc_reg$adj.values)
    results_abc_reg[i,2] <- test_reg[i,48]
  }
  
  return(results_abc_reg)
}

# Calculate factor 2 and other measures of accuracy

rs_seur <- make_regression(statistics_D2n_f,parameters, "sel_eur")
sum(rs_seur$True_param >= (rs_seur$Mean_posterior * 0.5 ) & rs_seur$True_param <= (rs_seur$Mean_posterior * 2 )) / length(rs_seur$True_param) 
sqrt(mean(rs_seur$Mean_posterior - rs_seur$True_param )^2)
cor.test(rs_seur$Mean_posterior,rs_seur$True_param, method = "spearman")

rs_nea <- make_regression(statistics_D2n_f,parameters, "sel_nea")
sum(rs_nea$True_param >= (rs_nea$Mean_posterior * 0.5 ) & rs_nea$True_param <= (rs_nea$Mean_posterior * 2 )) / length(rs_nea$True_param)
sqrt(mean(rs_nea$Mean_posterior - rs_nea$True_param )^2)
cor.test(rs_nea$Mean_posterior,rs_nea$True_param, method = "spearman")

rs_den <- make_regression(statistics_D2n_f,parameters, "sel_den")
sum(rs_den$True_param >= (rs_den$Mean_posterior * 0.5 ) & rs_den$True_param <= (rs_den$Mean_posterior * 2 )) / length(rs_den$True_param)
sqrt(mean(rs_den$Mean_posterior - rs_den$True_param )^2)
cor.test(rs_den$Mean_posterior,rs_den$True_param, method = "spearman")

rs_chb <- make_regression(statistics_D2n_f,parameters, "sel_chb")
sum(rs_chb$True_param >= (rs_chb$Mean_posterior * 0.5 ) & rs_chb$True_param <= (rs_chb$Mean_posterior * 2 )) / length(rs_chb$True_param)
sqrt(mean(rs_chb$Mean_posterior - rs_chb$True_param )^2)
cor.test(rs_chb$Mean_posterior,rs_chb$True_param, method = "spearman")

rs_time <- make_regression(statistics_D2n_f,parameters, "onset_time")
sum(rs_time$True_param >= (rs_time$Mean_posterior * 0.5 ) & rs_time$True_param <= (rs_time$Mean_posterior * 2 )) / length(rs_time$True_param)
sqrt(mean(rs_time$Mean_posterior - rs_time$True_param )^2)
cor.test(rs_time$Mean_posterior,rs_time$True_param, method = "spearman")

rs_ooa <- make_regression(statistics_D2n_f,parameters, "sel_ooa")
sum(rs_ooa$True_param >= (rs_ooa$Mean_posterior * 0.5 ) & rs_ooa$True_param <= (rs_ooa$Mean_posterior * 2 )) / length(rs_ooa$True_param)
sqrt(mean(rs_ooa$Mean_posterior - rs_ooa$True_param )^2)
cor.test(rs_ooa$Mean_posterior,rs_ooa$True_param, method = "spearman")

plot(rs_ooa$Mean_posterior, rs_ooa$True_param)
plot(rs_time$Mean_posterior, rs_time$True_param)
plot(rs_chb$Mean_posterior, rs_chb$True_param)
plot(rs_ooa$Mean_posterior, rs_ooa$True_param)


###### Neural net parameter estimation
reg_den <- merge(statistics_D2n_f, parameters, by = "replicate", sort = F) %>% 
  dplyr::select(-any_of(vector_selection),-model)

abc_nnet_tra <- cv4abc(param = reg_den[,48:53], sumstat = reg_den[,1:47], nval = 100, tols = 0.1, method = "neuralnet")

process_abc <- function(cvabc_out, parameter_name, indx) {
  trve <- abc_nnet_tra$true[,parameter_name] 
  post <- abc_nnet_tra$estim[[indx]][,parameter_name]
  factor_2 <- sum(trve >= (post * 0.5) & trve <= (post * 2)  ) / length(trve)
  rmse <- sqrt(mean(post - trve) ^2)
  Spearm <- cor.test(post,trve, method = "spearman")
  return(list(factor_2,rmse,Spearm))
}

seur_stats_01 <- process_abc(abc_nnet_tra, "sel_eur", 1)
schb_stats_01 <- process_abc(abc_nnet_tra, "sel_chb", 1)
sden_stats_01 <- process_abc(abc_nnet_tra, "sel_den", 1)
snea_stats_01 <- process_abc(abc_nnet_tra, "sel_nea", 1)
time_stats_01 <- process_abc(abc_nnet_tra, "onset_time", 1)
ooa_stats_01 <- process_abc(abc_nnet_tra, "sel_ooa", 1) 

# END
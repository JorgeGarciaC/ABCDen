### Auxiliary functions for the ABC
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
                      "fr_deni_old","fr_charg","fr_altai","fr_deni","fr_vin"
                      ,"ZnS_CHB","ZnS_EUR","ZnS_EAF","ZnS_DEN","ZnS_NEA", 
                      "replicate")

param_names <- c("sel_den","sel_nea","sel_chb","sel_eur","onset_time","replicate","mutation","time_simulation")

#Preprocess Reference tables
preprocess <- function(tsv_path, POP) {
  statistics <- read_tsv(tsv_path, col_names = names_statistics )
  statistics_f <- statistics[statistics$fr_allel_EUR > 0 & statistics$fr_allel_CHB > 0,]
  statistics_f$model <- POP
  return(statistics_f)
}

#Posterior estimation using RF
posterior_parameter_RF <- function(param,training,obs, vector_stats_remove=vector_selection) {
  training <- training %>% dplyr::select(c(names_statistics[-53], param, -any_of(vector_stats_remove)))
  model_regRF <- regAbcrf(formula = as.formula(paste(param,".", sep = " ~ ")),
                          data = training,
                          ntree = 500, ncores = 10, paral = T)
  posterior <- predict(object = model_regRF,
                       obs = obs,
                       training = training,
                       paral = T,
                       rf.weights = T)
  results_ABCRF <-list(model_regRF,posterior)
  return(results_ABCRF)
}


# making PCA
to_select <- c("fD_CHB_DEN","R_d_CHB_D","D_CHB_DEN","U_CHB_DEN","fD_EUR_NEA","fr_allel_EAF")
make_pca <- function (train,obs, subsetting = F, vector_to_select = to_select) {
  aux <- obs
  aux$model <- "obs"
  aux <- rbind(train,aux)
  labels_pca <- as.character(aux$model)

  if (subsetting == T) {
    aux <- aux %>% dplyr::select(all_of(vector_to_select))
  }
  
  pca_obj <- prcomp(data.frame(aux[,-length(aux)]), center = T, scale. = T)
  pca_df <- data.frame(pca_obj$x[,1:2])
  pca_df$model <- labels_pca
  
  p_pca <- ggplot(pca_df, aes(PC1,PC2, col=model)) + geom_point(shape=4) + 
    geom_point( aes(pca_df[length(labels_pca),1],pca_df[length(labels_pca),2]), size = 7,shape=19) + 
    scale_color_manual(values=c("#B1ABC3", "#867BA0","#3C1871","black")) +
    theme_bw() + guides(color=guide_legend(override.aes=list(shape=19), title = "Model"))
  return(p_pca)
}

# Get goodness of fit for different statistics for different populations gfit
get_adjust <- function(stats_observed,model_selection, POP, select_stats, rpl=100){
  den_gfit <- gfit(stats_observed[select_stats], 
                   model_selection[model_selection$model == POP,select_stats],
                   statistic=mean,nb.replicate=rpl)
  name <- colnames(model_selection)[i]
  p <- plot(den_gfit, main = paste0("Histogramme of the null distribution for statistic ", name, " in ",POP))
  #return(summary(den_gfit))
  return(den_gfit)
}



# Compare if distributions of parameters are updated after the ABC
make_histogram_parameter <- function(reference_and_param, results_posterior,parameter) {
  df_prior <- data.frame(param = reference_and_param[[parameter]], type = "Prior")
  df_posterior <- data.frame(param = reference_and_param[[parameter]], weight = results_posterior[[2]]$weights, type = "Posterior")
  df_combined <- rbind(df_prior, df_posterior[,-2])
  p <- ggplot(df_combined, aes(x = param, fill = type)) +
    geom_histogram(data = df_prior, aes(y = after_stat(density)), bins = 30, alpha = 0.5, color = "black") + 
    geom_histogram(data = df_posterior, aes(y = after_stat(density), weight = weight), bins = 30, alpha = 0.5, color = "black") + 
    scale_fill_manual(values = c("Prior" = "#B1ABC3", "Posterior" = "#1E1247")) + 
    #labs(x = "s in eur", y = "Probability Density", fill = "Distribution") +
    theme_minimal()
  return(p)
  
}


# The same but for normal abc output
make_histogram_abcclassic <- function(reference_parameters,results_posterior, paramer) {
  df_prior <- data.frame(param = reference_parameters[[paramer]], type = "Prior")
  df_posterior <- data.frame(param = results_posterior, type = "Posterior")
  p <- ggplot(df_prior, aes(x = param, fill = type)) +
    geom_histogram(data = df_prior, aes(y = after_stat(density)), bins = 30, alpha = 0.5, color = "black") + 
    geom_histogram(data = df_posterior, aes(y = after_stat(density)), bins = 30, alpha = 0.5, color = "black") + 
    scale_fill_manual(values = c("Prior" = "#B1ABC3", "Posterior" = "#1E1247")) + 
    #labs(x = "s in eur", y = "Probability Density", fill = "Distribution") +
    theme_classic()
  
  return(p)
}

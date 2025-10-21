##### Final ABC DEN, NEA, ANC

library(tidyverse)
library(abc)
library(stats)
library(gridExtra)
library(caret)
rm(list=ls())

set.seed(123)

source("/home/jorge/Rscripts/AuxAbc.R")

chunk <- "CORE"
vector_selection <- c("S_DEN","fr_deni_old","ZnS_DEN","replicate","ZnS_NEA", 
                        "fr_charg","fr_altai","fr_vin")

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
out_dir <- paste0(dir, "/RESULTS.NEG/")
train_path <- "/home/jorge/PROJECTS.JORGE/SLENDR/WORK_STATS/"

statistics_Deni_f <- preprocess(paste0(train_path,"DEN.NEG.TEST.stats.tsv"), "DENI")
statistics_Nea_f <- preprocess(paste0(train_path,"NEA.NEG.TEST.stats.tsv"), "NEA")
statistics_Anc_f <- preprocess(paste0(train_path,"ANC.NEG.TEST.stats.tsv"), "ANC")
statistics_D2n_f <- preprocess(paste0(train_path,"D2N.NEG.TEST.stats.tsv"), "D2N")

reference_table <- rbind(statistics_Deni_f,statistics_Nea_f, statistics_Anc_f,statistics_D2n_f) 

param_names <- c("sel_den","sel_nea","sel_chb","sel_ooa","sel_eur","onset_time","replicate","mutation","time_simulation")
parameters <- read_tsv(paste0(train_path, "parameters.D2N.NEG.TEST.stats.tsv"),
                         col_names = param_names, skip = 1)
# With so low samples sizes and small chunks, S may be 0 for NEA and then ZnS is NA
stats_observed <- read_tsv("/home/jorge/PROJECTS.JORGE/SIMULATIONS/WITH.ALLEL.FREQUENCIES/obs.CORE.stats.separated.70kb.tsv", 
                           col_names  = names_statistics) %>% dplyr::select(-any_of(vector_selection))

#most_important_impurity <- c("fr_al_OOA","fr_deni", "fD_EUR_DEN", "fD_CHB_DEN" ,"R_d_EUR_D" )#,"H1_CHB", "H12_CHB")
most_important_permutation <- c("fr_al_OOA", "fr_deni","fD_EUR_DEN", "fD_CHB_DEN","R_d_CHB_D", "R_d_EUR_D", "fr_al_nea", "D_EUR_DEN", "D_CHB_DEN")#, "S_NEA","D_CHB_NEA")

# Model discrimination
model_selection <- reference_table[,-53] %>% dplyr::select(-any_of(vector_selection))
model_selection$model <- as.factor(model_selection$model)

p_nosubset <- make_pca(model_selection,stats_observed, subsetting = F) + scale_color_manual(values=c("darkgrey","#3C1871","purple","black","lightgreen"))
#p_subset <- make_pca(model_selection,stats_observed, subsetting = T, vector_to_select = most_important_impurity) + scale_color_manual(values=c("darkgrey","#3C1871","purple","black","lightgreen"))
p_RF_perm <- make_pca(model_selection,stats_observed, subsetting = T, vector_to_select = most_important_permutation) + scale_color_manual(values=c("darkgrey","#3C1871","purple","black","lightgreen"))

lay <- rbind(c(1),
             c(2))

abc_rf <- grid.arrange(p_nosubset,p_RF_perm, layout_matrix= lay) # Bien pero mejorar 
ggsave(filename = paste0("PCA.results.",chunk,".png"), path = out_dir, plot = abc_rf, device = "png", units = "px", height = 3000, width = 3200)

# Posterior probabilities
set.seed(123)
post_pr_rfl <- postpr(target = stats_observed  %>% dplyr::select(all_of(most_important_permutation)), 
                      index = reference_table$model, 
                      sumstat = model_selection  %>% dplyr::select(all_of(most_important_permutation)),
                      tol = 0.1, method = "mnlogistic", corr = T)
summ_rfl <- summary(post_pr_rfl, digits = 2)

### Because -Inf so it does not allow to make regression so we discard it
post_pr_rfl <- postpr(target = stats_observed  %>% dplyr::select(all_of(most_important_permutation)), 
                      index = reference_table$model[!(reference_table$model == "DENI" | reference_table$model == "ANC")], 
                      sumstat = model_selection[!(reference_table$model == "DENI" | reference_table$model == "ANC"),]  %>% dplyr::select(all_of(most_important_permutation)),
                      tol = 0.1, method = "mnlogistic", corr = T)
summ_rfl <- summary(post_pr_rfl,digits = 2)

# Adjusted values of stats
pdf(file=paste0(out_dir,"/adjust.values.", chunk, ".pdf"))
for (i in 1:(dim(model_selection)[2]-1)) {
  par(mfrow = c(4,1))
  try(den_fdchb <- get_adjust(stats_observed, model_selection, "DENI",c(i)))
  try(nea_fdchb <- get_adjust(stats_observed, model_selection, "NEA",c(i)))
  try(anc_fdchb <- get_adjust(stats_observed, model_selection, "ANC",c(i)))
  try(anc_fdchb <- get_adjust(stats_observed, model_selection, "D2N",c(i)))
  
}
dev.off()

pdf(file=paste0(out_dir,"/Boxplots.simulations.",chunk, ".pdf"))
for (i in 1:(dim(model_selection)[2]-1)) {
  statistic_name <- colnames(model_selection)[i]
  print(statistic_name)
  md_longer <- pivot_longer(model_selection[,c(i,length(model_selection))], cols = any_of(statistic_name))
  boxplots <- ggplot(md_longer, aes(x=model, y=value,fill=model)) + 
    geom_boxplot(position = position_dodge(1)) + 
    geom_hline(yintercept = as.numeric(stats_observed[,i]),linetype = "dashed") +
    theme_classic() + labs(x = "Importance", y = statistic_name)
  print(boxplots)
}
dev.off()

###### Posterior parameter 
ref_infer <- merge(statistics_D2n_f,parameters, by="replicate") #%>% select(-S_DEN,-fr_deni_old,-ZnS_DEN)

make_plot_with_KS <- function(prior_ref_table,results_from_abc, parameter_name) {
  KS_r <- ks.test(prior_ref_table[,parameter_name],results_from_abc$adj.values[,parameter_name])
  y <- make_histogram_abcclassic(prior_ref_table,results_from_abc$adj.values[,parameter_name],parameter_name)
  y <- y + labs(x = parameter_name, y = "Probability Density", fill = "Distribution", title = paste0("KS_D =", round(KS_r$statistic,2),"; p-value = ", round(KS_r$p.value,2)))
  return(y)
}

best_seur <- c("fr_allel_EUR","S_EUR","R_d_EUR_N","R_d_EUR_N","D_EUR_DEN")
abc_RF_seur <- abc(stats_observed %>% dplyr::select(any_of(best_seur)), 
                   ref_infer[,"sel_eur"], 
                   ref_infer %>% dplyr::select(any_of(best_seur)),
                   method = "loclinear", tol = 0.1, corr = F)
y <- make_histogram_abcclassic(ref_infer,abc_RF_seur$adj.values,"sel_eur")  
KS_seur <- ks.test(ref_infer[,"sel_eur"],abc_RF_seur$adj.values)
y <- y + labs(x = "s in eur", y = "Probability Density", fill = "Distribution",title = paste0("KS_D = ", round(KS_seur$statistic,2),"; p-value = ", round(KS_seur$p.value,2)))

best_schb <- c("fr_allel_CHB", "tw_EUR", "S_CHB", "tw_CHB", "R_d_CHB_D")
abc_RF_schb <- abc(stats_observed %>% dplyr::select(any_of(best_schb)), 
                   ref_infer[,"sel_chb"], 
                   ref_infer %>% dplyr::select(any_of(best_schb)),
                   method = "loclinear", tol = 0.1, corr = F)
x <- make_histogram_abcclassic(ref_infer,abc_RF_schb$adj.values,"sel_chb")
KS_schb <- ks.test(ref_infer[,"sel_chb"],abc_RF_schb$adj.values)
x <- x + labs(x = "s in chb", y = "Probability Density", fill = "Distribution", title = paste0("KS_D = ", round(KS_schb$statistic,2),"; p-value = ", round(KS_schb$p.value,2)))

best_sden <- c("S_CHB","tw_CHB", "fr_allel_CHB", "R_d_CHB_N", "tw_EUR")
abc_RF_sden <- abc(stats_observed %>% dplyr::select(any_of(best_sden)), 
                   ref_infer[,"sel_den"], 
                   ref_infer %>% dplyr::select(any_of(best_sden)),
                   method = "loclinear", tol = 0.1, corr = F)
f <- make_histogram_abcclassic(ref_infer,abc_RF_sden$adj.values,"sel_den")
KS_sden <- ks.test(ref_infer[,"sel_den"],abc_RF_sden$adj.values)
f <- f + labs(x = "s in den", y = "Probability Density", fill = "Distribution", title = paste0("KS_D = ", round(KS_sden$statistic,2),"; p-value = ", round(KS_sden$p.value,2)))

best_snea <- c("fr_al_nea", "R_d_CHB_N","fD_CHB_DEN","fD_EUR_NEA","D_CHB_NEA")
abc_RF_snea <- abc(stats_observed %>% dplyr::select(any_of(best_snea)), 
                   ref_infer[,"sel_nea"], 
                   ref_infer %>% dplyr::select(any_of(best_snea)),
                   method = "loclinear", tol = 0.1, corr = F)
j <- make_histogram_abcclassic(ref_infer,abc_RF_snea$adj.values,"sel_nea")
KS_snea <- ks.test(ref_infer[,"sel_nea"],abc_RF_snea$adj.values)
j <- j + labs(x = "s in nea", y = "Probability Density", fill = "Distribution", title = paste0("KS_D = ", round(KS_snea$statistic,2),"; p-value = ", round(KS_snea$p.value,2)))

best_time <- c("fD_EUR_DEN", "fD_CHB_DEN", "R_d_CHB_D", "D_CHB_DEN","R_d_EUR_D")
abc_RF_time <- abc(stats_observed %>% dplyr::select(any_of(best_time)), 
                   ref_infer[,"onset_time"], 
                   ref_infer %>% dplyr::select(any_of(best_time)),
                   method = "loclinear", tol = 0.1, corr = F)
e <- make_histogram_abcclassic(ref_infer,abc_RF_time$adj.values,"onset_time")
KS_time <- ks.test(ref_infer[,"onset_time"],abc_RF_time$adj.values)
e <- e + labs(x = "onset time", y = "Probability Density", fill = "Distribution", title = paste0("KS_D = ", round(KS_time$statistic,2),"; p-value = ", round(KS_time$p.value,2)))

best_ooa <- c("fr_allel_OOA","H1_EUR", "tw_CHB","S_CHB","R_d_CHB_D")
abc_RF_ooa <- abc(stats_observed %>% dplyr::select(any_of(best_ooa)), 
                  ref_infer[,"sel_ooa"], 
                  ref_infer %>% dplyr::select(any_of(best_ooa)),
                   method = "loclinear", tol = 0.1, corr = F)
g <- make_histogram_abcclassic(ref_infer,abc_RF_ooa$adj.values,"sel_ooa")
KS_ooa <- ks.test(ref_infer[,"sel_ooa"],abc_RF_ooa$adj.values)
g <- g + labs(x = "s in OOA", y = "Probability Density", fill = "Distribution", title = paste0("KS_D = ", round(KS_time$statistic,2),"; p-value = ", round(KS_time$p.value,2)))


Parameter_estimation_RF <- grid.arrange(y,x,j,f,e,g, layout_matrix = rbind(c(1,2),c(3,4),c(5,6) ))
ggsave(filename = paste0("ABC_RF.parameter_estimation.results.KS",chunk,".png"), path = out_dir, plot = Parameter_estimation_RF, device = "png", units = "px", height = 3000, width = 3200)

dispersion_statistics <- function(distribution){
  return(c(mean(distribution),median(distribution),sd(distribution),
           quantile(distribution, probs = c(0.025,0.975))))
           
}

st_ooa <- dispersion_statistics(abc_RF_ooa$adj.values)
st_eur <- dispersion_statistics(abc_RF_seur$adj.values)
st_chb <- dispersion_statistics(abc_RF_schb$adj.values)
st_den <- dispersion_statistics(abc_RF_sden$adj.values)
st_nea <- dispersion_statistics(abc_RF_snea$adj.values)
st_time <- dispersion_statistics(abc_RF_time$adj.values)


###### Reinference from posterior parameters

post_stats_2dn <- preprocess(paste0(train_path,"post_d2n.stats.tsv"), "D2N")
post_params_2dn <- read_tsv(paste0(train_path, "post_d2n.params.tsv"),
                       col_names = param_names)

ref_infer_post <- merge(post_stats_2dn,post_params_2dn, by="replicate") #%>% select(-S_DEN,-fr_deni_old,-ZnS_DEN)

abc_RF_seur_post <- abc(stats_observed %>% dplyr::select(any_of(best_seur)), 
                   ref_infer_post[,"sel_eur"], 
                   ref_infer_post %>% dplyr::select(any_of(best_seur)),
                   method = "loclinear", tol = 0.1, corr = F)
y_post <- make_histogram_abcclassic(post_params_2dn,abc_RF_seur_post$adj.values,"sel_eur")  
KS_seur_post <- ks.test(post_params_2dn[,"sel_eur"],abc_RF_seur_post$adj.values)
y_post <- y_post + labs(x = "s in eur", y = "Probability Density", fill = "Distribution",title = paste0("KS_D = ", round(KS_seur_post$statistic,2),"; p-value = ", round(KS_seur_post$p.value,2)))

abc_RF_schb_post <- abc(stats_observed %>% dplyr::select(any_of(best_schb)), 
                   ref_infer_post[,"sel_chb"], 
                   ref_infer_post %>% dplyr::select(any_of(best_schb)),
                   method = "loclinear", tol = 0.1, corr = F)
x_post <- make_histogram_abcclassic(post_params_2dn,abc_RF_schb_post$adj.values,"sel_chb")
KS_schb_post <- ks.test(post_params_2dn[,"sel_chb"],abc_RF_schb_post$adj.values)
x_post <- x_post + labs(x = "s in chb", y = "Probability Density", fill = "Distribution", title = paste0("KS_D = ", round(KS_schb_post$statistic,2),"; p-value = ", round(KS_schb_post$p.value,2)))

abc_RF_sden_post <- abc(stats_observed %>% dplyr::select(any_of(best_sden)), 
                   ref_infer_post[,"sel_den"], 
                   ref_infer_post %>% dplyr::select(any_of(best_sden)),
                   method = "loclinear", tol = 0.1, corr = F)
f_post <- make_histogram_abcclassic(post_params_2dn,abc_RF_sden_post$adj.values,"sel_den")
KS_sden_post <- ks.test(post_params_2dn[,"sel_den"],abc_RF_sden_post$adj.values)
f_post <- f_post + labs(x = "s in den", y = "Probability Density", fill = "Distribution", title = paste0("KS_D = ", round(KS_sden_post$statistic,2),"; p-value = ", round(KS_sden_post$p.value,2)))


abc_RF_snea_post <- abc(stats_observed %>% dplyr::select(any_of(best_snea)), 
                   ref_infer_post[,"sel_nea"], 
                   ref_infer_post %>% dplyr::select(any_of(best_snea)),
                   method = "loclinear", tol = 0.1, corr = F)
j_post <- make_histogram_abcclassic(post_params_2dn,abc_RF_snea_post$adj.values,"sel_nea")
KS_snea_post <- ks.test(post_params_2dn[,"sel_nea"],abc_RF_snea_post$adj.values)
j_post <- j_post + labs(x = "s in nea", y = "Probability Density", fill = "Distribution", title = paste0("KS_D = ", round(KS_snea_post$statistic,2),"; p-value = ", round(KS_snea_post$p.value,2)))

abc_RF_time_post <- abc(stats_observed %>% dplyr::select(any_of(best_time)), 
                   ref_infer_post[,"onset_time"], 
                   ref_infer_post %>% dplyr::select(any_of(best_time)),
                   method = "loclinear", tol = 0.1, corr = F)
e_post <- make_histogram_abcclassic(post_params_2dn,abc_RF_time_post$adj.values,"onset_time")
KS_time_post <- ks.test(post_params_2dn[,"onset_time"],abc_RF_time_post$adj.values)
e_post <- e_post + labs(x = "onset time", y = "Probability Density", fill = "Distribution", title = paste0("KS_D = ", round(KS_time_post$statistic,2),"; p-value = ", round(KS_time_post$p.value,2)))

abc_RF_ooa_post <- abc(stats_observed %>% dplyr::select(any_of(best_ooa)), 
                  ref_infer_post[,"sel_ooa"], 
                  ref_infer_post %>% dplyr::select(any_of(best_ooa)),
                  method = "loclinear", tol = 0.1, corr = F)
g_post <- make_histogram_abcclassic(post_params_2dn,abc_RF_ooa_post$adj.values,"sel_ooa")
KS_ooa_post <- ks.test(post_params_2dn[,"sel_ooa"],abc_RF_ooa_post$adj.values)
g_post <- g_post + labs(x = "s in OOA", y = "Probability Density", fill = "Distribution", title = paste0("KS_D = ", round(KS_ooa_post$statistic,2),"; p-value = ", round(KS_ooa_post$p.value,2)))

Parameter_estimation_RF_post <- grid.arrange(y_post,x_post,j_post,f_post,e_post,g_post, layout_matrix = rbind(c(1,2),c(3,4),c(5,6) ))
ggsave(filename = paste0("ABC_RF.parameter_estimation.results.KS_post_2",chunk,".png"), path = out_dir, plot = Parameter_estimation_RF_post, device = "png", units = "px", height = 3000, width = 3200)

st_ooa_post <- dispersion_statistics(abc_RF_ooa_post$adj.values)
st_eur_post <- dispersion_statistics(abc_RF_seur_post$adj.values)
st_chb_post <- dispersion_statistics(abc_RF_schb_post$adj.values)
st_den_post <- dispersion_statistics(abc_RF_sden_post$adj.values)
st_nea_post <- dispersion_statistics(abc_RF_snea_post$adj.values)
st_time_post <- dispersion_statistics(abc_RF_time_post$adj.values)

## Poster figures

p_nosubset <- make_pca(model_selection,stats_observed, subsetting = F) + scale_color_manual(values=c("#B1ABC3", "#1e1247", "#867BA0","#3C1871", "#00C853")) 
p_subset <- make_pca(model_selection,stats_observed, subsetting = T, vector_to_select = most_important_permutation) + scale_color_manual(values=c("#B1ABC3", "#1e1247", "#867BA0","#3C1871", "#00C853"))
lgd <- cowplot::get_legend(p_nosubset)

plot_poster <- cowplot::plot_grid(p_nosubset + theme(legend.position = "none") , p_subset + theme(legend.position = "none"), lgd, ncol = 3, rel_widths = c(3,3,1))
ggsave(filename = paste0("PCA.poster",chunk,".png"), path = out_dir, plot = plot_poster, device = "png", units = "px", height = 1200, width = 2200)

## paper figure
parameters_legend <- cowplot::get_legend(x_post)
cowplot::plot_grid(y_post + labs(title ="") + theme(legend.position = "none"),
                   x_post + labs(title ="") + theme(legend.position = "none"),
                   g_post + labs(title ="") + theme(legend.position = "none"),
                   parameters_legend,
                   f_post + labs(title ="") + theme(legend.position = "none"),
                   j_post + labs(title ="") + theme(legend.position = "none"),
                   e_post + labs(title ="") + theme(legend.position = "none"), 
                    nrow = 2, rel_widths = c(3,3,3,1,3,3,3), labels = c("A","B","C","","D","E","F"))
ggsave(filename = paste0("ABC_RF.parameter_estimation.labels_post_2",chunk,".png"), path = out_dir, plot = last_plot(), device = "png", units = "px", height = 2400, width = 3200)

# Posteriors pal poster
cowplot::plot_grid(y + labs(title ="") + theme(legend.position = "none"),
                   x + labs(title ="", x = "s in eas") + theme(legend.position = "none"),
                   parameters_legend, nrow = 1, rel_widths = c(3,3,1))
ggsave(filename = paste0("ABC_RF.poster",chunk,".png"), path = out_dir, plot = last_plot(), device = "png", units = "px", height = 1200, width = 3200)

save.image(paste0(out_dir,"ABC_FINAL.4models.R"))

load(paste0(out_dir,"ABC_FINAL.4models.R"))

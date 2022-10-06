library(dplyr)
library(car)
library(eulerr)

df_all_data = read.csv(file = "path/to/this/file/data_for_reproduction.csv")
outcome_vars <- c('FEV1_utah', 
                  "SGRQ_scoreTotal", 
                  "distwalked", 
                  "FEV1_FVC_utah", 
                  "FRC_Thirona",
                  "FRC_TLC_ratio_Thirona")

###  Code to reproduce regression modeling results in Tables 3 and 4, and Figures 5 and 6
fits_all_quant <-  bind_rows(lapply(1:length(outcome_vars), function(idx_var){
  voi_sub <- outcome_vars[idx_var]
  dat_fit <- df_all_data
  dat_fit$voi = df_all_data[, voi_sub]
  
  dat_fit$lemph = logit(dat_fit$ovr_pLAA)
  dat_fit$lemph = (dat_fit$lemph - mean(dat_fit$lemph, na.rm = T)) / sd(dat_fit$lemph, na.rm = T)
  dat_fit$d_mle = (dat_fit$d_mle - mean(dat_fit$d_mle, na.rm = T)) / sd(dat_fit$d_mle, na.rm = T)
  dat_fit$mean_avg_cluster_size = (dat_fit$mean_avg_cluster_size - mean(dat_fit$mean_avg_cluster_size, na.rm = T)) / sd(dat_fit$mean_avg_cluster_size, na.rm = T)
  fit_sub <- lm(voi ~ as.factor(gender) + age_visit + BMI + Height_CM + SmokCigNow + lemph + mean_avg_cluster_size + d_mle , data = dat_fit)
  fit_red <- lm(voi ~ as.factor(gender) + age_visit + BMI + Height_CM + SmokCigNow + lemph +  d_mle , data = dat_fit)
  fit_red2 <- lm(voi ~ as.factor(gender) + age_visit + BMI + Height_CM + SmokCigNow + mean_avg_cluster_size  , data = dat_fit)
  ssub <- summary(fit_sub)
  ssub_red = summary(fit_red)
  ssub_red2 = summary(fit_red2)
  asub = anova(fit_red, fit_sub, test = "LRT")
  asub2 = anova(fit_red2, fit_sub, test = "LRT")
  ret <- data.frame(voi = voi_sub, 
                    estimate_pLAA = as.numeric(ssub$coefficients[7, 1]), 
                    se_pLAA = as.numeric(ssub$coefficients[7, 2]), 
                    p_val_pLAA = as.numeric(ssub$coefficients[7, 4]),
                    estimate_d_mle = as.numeric(ssub$coefficients[9, 1]), 
                    se_d_mle = as.numeric(ssub$coefficients[9, 2]), 
                    p_val_d_mle = as.numeric(ssub$coefficients[9, 4]), 
                    estimate_MA = as.numeric(ssub$coefficients[8, 1]), 
                    se_MA = as.numeric(ssub$coefficients[8, 2]), 
                    p_val_MA = as.numeric(ssub$coefficients[8, 4]), 
                    r2_all = as.numeric(ssub$adj.r.squared),
                    r2_red = as.numeric(ssub_red$adj.r.squared),
                    r2_red2 = as.numeric(ssub_red2$adj.r.squared),
                    p_aov = as.numeric(asub$`Pr(>Chi)`[2]),
                    p_aov2 = as.numeric(asub2$`Pr(>Chi)`[2]))
  
  return(ret)
}))

print(fits_all_quant)

###   Correlations presented in Table 2
cor(df_all_data[, c(26, 25, 23, 20, 22)], method = "pearson")

###   Code for reproducing p-values presented in Figure 6 panel associated with visual assessment
fits_all_vis <-  bind_rows(lapply(1:length(outcome_vars), function(idx_var){
  voi_sub <- outcome_vars[idx_var]
  dat_fit <- df_all_data
  dat_fit$voi = df_all_data[, voi_sub]
  
  dat_fit$mean_avg_cluster_size = (dat_fit$mean_avg_cluster_size - mean(dat_fit$mean_avg_cluster_size, na.rm = T)) / sd(dat_fit$mean_avg_cluster_size, na.rm = T)
  
  fit_sub <- lm(voi ~ as.factor(gender) + age_visit + BMI + Height_CM + SmokCigNow + as.factor(CT_Visual_Emph_Severity) + as.factor(CT_Visual_Emph_Paraseptal) + mean_avg_cluster_size, data = dat_fit)
  fit_red <- lm(voi ~ as.factor(gender) + age_visit + BMI + Height_CM + SmokCigNow + as.factor(CT_Visual_Emph_Severity) + as.factor(CT_Visual_Emph_Paraseptal), data = dat_fit)
  fit_red2 <- lm(voi ~ as.factor(gender) + age_visit + BMI + Height_CM + SmokCigNow + mean_avg_cluster_size  , data = dat_fit)
  ssub <- summary(fit_sub)
  ssub_red = summary(fit_red)
  ssub_red2 = summary(fit_red2)
  asub = anova(fit_red, fit_sub, test = "LRT")
  asub2 = anova(fit_red2, fit_sub, test = "LRT")
  ret <- data.frame(voi = voi_sub,
                    r2_all = as.numeric(ssub$adj.r.squared),
                    r2_red = as.numeric(ssub_red$adj.r.squared),
                    r2_red2 = as.numeric(ssub_red2$adj.r.squared),
                    p_aov = as.numeric(asub$`Pr(>Chi)`[2]),
                    p_aov2 = as.numeric(asub2$`Pr(>Chi)`[2]))
  
  return(ret)
}))

print(fits_all_vis)

RBM_data_both_ways_SJ <- read.delim("path/to/this/file/RBM_data_both_ways_SJ.txt", 
                                    stringsAsFactors=FALSE)

###   Code for reproducing regression modeling results related to plasma protein array
res_lrts <- bind_rows(lapply(2:110, function(idx_col){
  bm_sub <- names(RBM_data_both_ways_SJ)[idx_col]
  dat_fit <- RBM_data_both_ways_SJ[, c(1, idx_col)] %>% 
    right_join(df_all_data) 
  names(dat_fit)[2] <- 'expr'
  
  if(idx_col < 97){
    fit_full <- lm(expr ~ logit(ovr_pLAA) + d_mle + mean_avg_cluster_size + 
                     mean_num_clust +
                     BMI + SmokCigNow + age_visit + as.factor(gender), data = dat_fit)
    fit_old = lm(expr ~ logit(ovr_pLAA) + d_mle +  BMI + SmokCigNow + age_visit + as.factor(gender), data = dat_fit)
    fit_new = lm(expr ~ mean_avg_cluster_size + 
                   mean_num_clust +
                   BMI + SmokCigNow + age_visit + as.factor(gender), data = dat_fit)
    fit_red <- lm(expr ~ BMI + SmokCigNow + age_visit + as.factor(gender), data = dat_fit)
    asub = anova(fit_red, fit_full, test = "LRT")
    asub2 = anova(fit_red, fit_old, test = "LRT")
    asub3 = anova(fit_red, fit_new, test = "LRT")
    ret = data.frame(biomarker = bm_sub,
                     p_aov_full = asub$`Pr(>Chi)`[2],
                     p_aov_old = asub2$`Pr(>Chi)`[2],
                     p_aov_new = asub3$`Pr(>Chi)`[2])
    return(ret)
  }
  else{
    dat_fit$expr <- 1 * (dat_fit$expr == 2)
    fit_full <- glm(expr ~ logit(ovr_pLAA) + d_mle + mean_avg_cluster_size + 
                      mean_num_clust +
                      BMI + SmokCigNow + age_visit + as.factor(gender), data = dat_fit, family = 'binomial')
    fit_old <- glm(expr ~ logit(ovr_pLAA) + d_mle + BMI + SmokCigNow + age_visit + as.factor(gender), data = dat_fit, family = 'binomial')
    fit_new <- glm(expr ~ mean_avg_cluster_size + 
                     mean_num_clust +
                     BMI + SmokCigNow + age_visit + as.factor(gender), data = dat_fit, family = 'binomial')
    fit_red <- glm(expr ~ BMI + SmokCigNow + age_visit + as.factor(gender), data = dat_fit, family = 'binomial')
    asub = anova(fit_red, fit_full, test = "LRT")
    asub2 = anova(fit_red, fit_old, test = "LRT")
    asub3 = anova(fit_red, fit_new, test = "LRT")
    ret = data.frame(biomarker = bm_sub,
                     p_aov_full = asub$`Pr(>Chi)`[2],
                     p_aov_old = asub2$`Pr(>Chi)`[2],
                     p_aov_new = asub3$`Pr(>Chi)`[2])
    return(ret)
  }
}))

head(res_lrts)

res_lrts_adj = res_lrts %>% 
  mutate(adj_p_full = p.adjust(p_aov_full, method = "BH"),
         adj_p_old = p.adjust(p_aov_old, method = "BH"),
         adj_p_new = p.adjust(p_aov_new, method = "BH")) %>% 
  mutate(new_better = adj_p_full < adj_p_old)

colSums(res_lrts_adj[, -1] < .10)

res_sig = res_lrts_adj %>% filter(adj_p_full < .10) 
table(res_sig$new_better)
res_sig %>% arrange(adj_p_full) %>% filter(adj_p_old > .10)

df_vd = cbind(pLAA_D = res_lrts_adj$adj_p_old < .10, SNCP = res_lrts_adj$adj_p_new < .10, all = res_lrts_adj$adj_p_full < .10)
colSums(df_vd)

res_sig_total = res_lrts_adj %>% filter(adj_p_full < .10 | adj_p_old < .10) 
res_sig_out = res_sig_total[, c(1, 5,6)] %>% arrange(adj_p_full)

names(res_sig_out) = c("Protein", "FDR_new", "FDR_old")
head(res_sig_out)

### Genes/proteins that go into Table 5
res_sig_out_sub = res_sig_out %>% filter(FDR_new < 0.10 & FDR_old > 0.10)
res_sig_out_sub

res_sig_out %>% filter(FDR_new < 0.10 & FDR_old > 0.10)
res_sig_out %>% filter(FDR_new > 0.10 & FDR_old < 0.10)

pval_frmt <- function(pvals, digits = 3){
  ret2 <- do.call(c, lapply(1:length(pvals), function(i){
    #i = 1
    psub <- pvals[i]
    if(is.na(psub)){
      return(" - ")
    }
    if(psub > 10^(-digits)){
      ret <- formatC(psub, format = 'f', digits = digits)
    }
    else{
      ret <- formatC(psub, format = 'e', digits = 1)
      ret <- strsplit(ret, 'e')
      ret <- paste0(paste0(ret[[1]], collapse = " x $10^{"), "}$")
    }
    return(ret)
  }))
  return(ret2)
}
res_sig_total_out = res_sig_out
res_sig_total_out[, 2:3] = apply(res_sig_total_out[, c(2, 3)], 
                                 MARGIN = 2, 
                                 FUN =  function(x){
                                   return(pval_frmt(x, digits = 3))
                                 })
###   Genes/proteins that go into Supplementary Table 1
res_sig_total_out



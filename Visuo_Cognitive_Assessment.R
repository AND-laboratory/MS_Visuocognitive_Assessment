# _______________________ LOAD PACKAGES __________________________________________

library(dplyr)
library(tidyr)
library(readxl)
library(writexl)
library(lme4)
library(lmerTest)
library(emmeans)
library(parameters)
library(performance)
library(sjPlot)
library(r2mlm)
library(lattice)
library(logistf)
library(ordinal)
library(brms)
library(posterior)
library(pROC)
library(rstatix)
library(car)
library(corrr)
library(psych)
library(effsize)
library(dunn.test)
library(bestNormalize)
library(ggplot2)
library(ggpubr)
library(visreg)
library(naniar)
library(finalfit)

#________________________ ANALYSIS 1 (NP Task Trajectories) ______________________

A1 <-read_excel("Longitudina_Data_Manuscript3.xlsx")

#extract baseline/V1 predictors
variables <- c("RNFL", "VG_Latency", "AS_Latency", "AS_Error", "AS_FEP", "VDI", 
               "NART","Age", "Sex", "Grp_Coded", 
               "SDMT_Imp", "PASAT_Imp", "CVLT_Imp")

visit1_data <- A1 %>%
  filter(Visit == 0) %>%
  dplyr::select(ID, all_of(variables))

visit1_data <- visit1_data %>%
  rename_with(~ paste0(.x, "_V1"), -ID)

A1 <- A1 %>%
  left_join(visit1_data, by = "ID")

### S D M T   L M M ###

#visit only model
SDMT_visit <- lmer(SDMT_z ~ Visit + (Visit | ID), data = A1)
summary(SDMT_visit)

#final model
SDMT_final <- lmer(SDMT_z ~ Visit * scale(VG_Latency_V1) + scale(AS_Error_V1) + scale(NART_V1) +  (Visit | ID), data = A1)
summary(SDMT_final)
check_collinearity(SDMT_final)
confint(SDMT_final, method = "Wald")
icc(SDMT_final)
model_performance(SDMT_final)
standardize_parameters(SDMT_final)

#final model diagnostics
#OVERALL
plot_model(SDMT_final, type = "diag")

#NORMALITY 
qqmath(SDMT_final) 

#HOMOSCEDASTICITY 
plot(fitted(SDMT_final), residuals(SDMT_final), xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs Fitted Values")
abline(h = 0, col = "red")

### P A S A T   L M M ###

#transformation to improve homoscedasticity
yeo_johnson_result <- yeojohnson(A1$PASAT_zc)
A1$yeo_PASAT <- yeo_johnson_result$x.t
hist(A1$yeo_PASAT)

#visit only model
PASAT_visit <- lmer(yeo_PASAT ~ Visit + (Visit | ID), data = A1)
summary(PASAT_visit)

#final model
PASAT_final <- lmer(yeo_PASAT ~ Visit + scale(AS_Error_V1) + PASAT_Imp_V1 + scale(NART_V1) + (Visit | ID), data = A1)
summary(PASAT_final)
confint(PASAT_final, method = "Wald")
check_collinearity(PASAT_final)
icc(PASAT_final)
model_performance(PASAT_final)
standardize_parameters(PASAT_final)

##final model diagnostics
#OVERALL
plot_model(PASAT_final , type = "diag")

#NORMALITY 
qqmath(PASAT_final ) 

#HOMOSCEDASTICITY 
plot(fitted(PASAT_final ), residuals(PASAT_final), xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs Fitted Values")
abline(h = 0, col = "red")

### C V L T    L M M ###

##random intercept only models (variance of random slope near 0) 

#visit only model
CVLT_visit <- lmer(CVLT_z ~ Visit + (1 | ID), data = A1)
summary(CVLT_visit)

#final model
CVLT_final <- lmer(CVLT_z ~ Visit + scale(VDI_V1) + CVLT_Imp_V1 + scale(NART_V1) + (1 | ID), data = A1)
summary(CVLT_final)
check_collinearity(CVLT_final)
confint(CVLT_final, method = "Wald")
icc(CVLT_final)
model_performance(CVLT_final)
standardize_parameters(CVLT_final)

##final model diagnostics
#OVERALL
plot_model(CVLT_final , type = "diag")

#NORMALITY 
qqmath(CVLT_final ) 

#HOMOSCEDASTICITY 
plot(fitted(CVLT_final ), residuals(CVLT_final), xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs Fitted Values")
abline(h = 0, col = "red")

# ________________________________________ ANALYSIS 2 (Cognitive Impairment) ____________________________________

A2 <-read_excel("Longitudina_Data_Manuscript3_wide.xlsx")

### L O G I S T I C   R E G R E S S I O N ###

table(A2$Imp_Binary)

#covariate model (no sig covariates)
firth_COV <-logistf(Imp_Binary ~ NART, data = A2, family = binomial)  
summary(firth_COV)

#visuo-cognitive model (no sig VP predictor)
firth_VP <- logistf(Imp_Binary ~ scale(VDI), data = A2, family = binomial)  
summary(firth_VP)

### O R D I N A L   R E G R E S S I O N ###

A2$Imp_Category <- factor(A2$Imp_Category, ordered = TRUE)
table(A2$Imp_Category)

#covariate model (no sig covariates)
ordinal_COV <- clm(Imp_Category ~ NART, data = A2)  
summary(ordinal_COV) 

#visuo-cognitive model
ordinal_VP <- clm(Imp_Category ~ scale(VDI), data = A2)
summary(ordinal_VP)
confint(ordinal_VP)

#OR 
coef_ordinal <- coef(ordinal_VP) 
conf_ordinal <- confint(ordinal_VP) 
odds_ratios <- exp(coef_ordinal) 
ci_lower <- exp(conf_ordinal[, 1])  
ci_upper <- exp(conf_ordinal[, 2])  
results <- data.frame(
  Predictor = names(coef_ordinal),
  Coefficient = coef_ordinal,
  Odds_Ratio = odds_ratios,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper)
print(results)

##final model diagnostics 
#proportional odds
prop_odds_test <- nominal_test(ordinal_VP)
print(prop_odds_test)

#ROC curves
#categories 0 vs 1
A2$binary_01 <- ifelse(A2$Imp_Category == 0, 0, ifelse(A2$Imp_Category == 1, 1, NA))
A2b_clean <- A2[!is.na(A2$binary_01) & !is.na(A2$VDI), ]
model_01 <- glm(binary_01 ~ VDI, data = A2b_clean, family = binomial)
prob_01 <- predict(model_01, type = "response")
roc_01 <- roc(A2b_clean$binary_01, prob_01)

print(auc(roc_01))
plot(roc_01, col = "blue", lwd = 2)
auc_value <- auc(roc_01)
abline(a = 0, b = 1, col = "red", lty = 2)
legend("bottomright", legend = paste("AUC =", round(auc_value, 3)), 
       col = "blue", lwd = 2)

#categories 1 vs 2
A2$binary_12 <- ifelse(A2$Imp_Category == 1, 0, ifelse(A2$Imp_Category == 2, 1, NA))
A2b_clean_12 <- A2[!is.na(A2$binary_12) & !is.na(A2$VDI), ]
model_12 <- glm(binary_12 ~ VDI, data = A2b_clean_12, family = binomial)
prob_12 <- predict(model_12, type = "response")
roc_12 <- roc(A2b_clean_12$binary_12, prob_12)

print(auc(roc_12))
plot(roc_12, col = "blue", lwd = 2)
auc_value <- auc(roc_12)
abline(a = 0, b = 1, col = "red", lty = 2)
legend("bottomright", legend = paste("AUC =", round(auc_value, 3)), 
       col = "blue", lwd = 2)

#categories 0 vs 2
A2$binary_02 <- ifelse(A2$Imp_Category == 0, 0, ifelse(A2$Imp_Category == 2, 1, NA))
A2b_clean_02 <- A2[!is.na(A2$binary_02) & !is.na(A2$VDI), ]
model_02 <- glm(binary_02 ~ VDI, data = A2b_clean_02, family = binomial)
prob_02 <- predict(model_02, type = "response")
roc_02 <- roc(A2b_clean_02$binary_02, prob_02)
print(auc(roc_02))

plot(roc_02, col = "blue", lwd = 2)
auc_value <- auc(roc_02)
abline(a = 0, b = 1, col = "red", lty = 2)
legend("bottomright", legend = paste("AUC =", round(auc_value, 3)), 
       col = "blue", lwd = 2)

#______________________________ ANALYSIS 3a (MRI) _________________________________

A3a <- read_excel("Manuscript3_MRI.xlsx")

### L E S I O N   F R A C T I O N  ### 

#transformation to improve homoscedasticity
yeo_johnson_result <- yeojohnson(A3a$lesion_fraction)
A3a$yeo_LF <- yeo_johnson_result$x.t
hist(A3a$yeo_LF)

#visit and covariate models
lf_visit <- lmer(yeo_LF ~ Visit + (Visit | ID), data = A3a) 
summary(lf_visit)

lf_cov <- lmer(yeo_LF ~ Visit + Grp_Coded_B + (Visit | ID), data = A3a) 
summary(lf_cov)

#neuropsych LMM (no sig NP predictor)
lf_NP <- lmer(yeo_LF  ~ Visit +  Grp_Coded_B + SDMT_z_B +  (Visit | ID), data = A3a) 
summary(lf_NP)
check_collinearity(lf_NP)

#visuo-cognitive LMM (Grp_Coded_B removed as no longer sig in final model)
lf_final <- lmer(yeo_LF  ~ Visit + scale(VG_Latency_B) +  (Visit | ID), data = A3a) 
summary(lf_final)
icc(lf_final)
model_performance(lf_final)
confint(lf_final , method = "Wald")
standardize_parameters(lf_final)

##final model diagnostics
#OVERALL
plot_model(lf_final, type = "diag")

#NORMALITY 
qqmath(lf_final ) 

#HOMOSCEDASTICITY 
plot(fitted(lf_final), residuals(lf_final), xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs Fitted Values")
abline(h = 0, col = "red")

### B P F ###

#visit and covariate models
bpf_visit <- lmer(BPF ~ Visit + (Visit | ID), data = A3a) 
summary(bpf_visit)

bpf_cov <- lmer(BPF ~ Visit + scale(Age_B) + (Visit | ID), data = A3a) 
summary(bpf_cov)

#neuropsych LMM (no sig NP predictor)
bpf_NP <- lmer(yeo_LF  ~ Visit + scale(Age_B) + SDMT_z_B +  (Visit | ID), data = A3a) 
summary(bpf_NP)
check_collinearity(bpf_NP)

#visuo-cognitive LMM
bpf_final <- lmer(BPF ~ Visit * scale(VG_Latency_B) + scale(Age_B) + (Visit | ID), data = A3a) 
summary(bpf_final)
check_collinearity(bpf_final)
icc(bpf_final)
model_performance(bpf_final)
confint(bpf_final , method = "Wald")
standardize_parameters(bpf_final)

#final model diagnostics
#OVERALL
plot_model(bpf_final , type = "diag")

#NORMALITY 
qqmath(bpf_final) 

#HOMOSCEDASTICITY 
plot(fitted(bpf_final ), residuals(bpf_final ), xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs Fitted Values")
abline(h = 0, col = "red")

#______________________________ ANALYSIS 3b (EDSS) ________________________________

A3b <- read_excel("Manuscript3_EDSS.xlsx")

#extract baseline/V1 predictors
variables <- c("RNFL","VG_Latency", "AS_Latency", "AS_Error", "AS_FEP", "VDI", 
              "SDMT_z",  "PASAT_zc", "CVLT_z",  
               "NART","Age", "Sex", "Grp_Coded")

visit1_edss <- A3b  %>%
  filter(Visit == 0) %>%
  dplyr::select(ID, all_of(variables))

visit1_edss <- visit1_edss %>%
  rename_with(~ paste0(.x, "_V1"), -ID)

A3b  <- A3b  %>%
  left_join(visit1_edss, by = "ID")

#convert EDSS to 3-level ordinal variable
A3b$EDSS_grp_3 <- cut(
  A3b$EDSS,
  breaks = c(-Inf, 0, 1.5 , Inf),  
  labels = c(0, 1, 2),           
  right = TRUE)

A3b$EDSS_grp_3 <- as.ordered(A3b$EDSS_grp_3)

#covariate ordinal MM 
A3b$Age_V1_scaled <- scale(A3b$Age_V1)

EDSS_ordinal_COV <- brm(
  formula = EDSS_grp_3 ~ Visit + Age_V1_scaled + (1 | ID), 
  data = A3b, 
  family = cumulative, seed = 123)
summary(EDSS_ordinal_COV)

#neuropsych ordinal MM (no sig NP predictor) 
EDSS_ordinal_NP <- brm(
  formula = EDSS_grp_3 ~ Visit + SDMT_z_V1 + Age_V1_scaled + (1 | ID), 
  data = A3b, 
  family = cumulative, seed = 123)
summary(EDSS_ordinal_NP)

#visuo-cognitive ordinal MM
A3b$RNFL_V1_scaled <- scale(A3b$RNFL_V1)

EDSS_ordinal <- brm(
  formula = EDSS_grp_3 ~ Visit + RNFL_V1_scaled + Age_V1_scaled  + (1 | ID),
  data = A3b,
  family = cumulative(),
  save_pars = save_pars(all = TRUE), seed = 123)
print(summary(EDSS_ordinal), digits = 6)

#OR
posterior_summary <- posterior_summary(EDSS_ordinal)
print(posterior_summary, digits = 6)

#DIAGNOSTICS
draws <- as_draws_array(EDSS_ordinal)  
bulk_ess <- ess_bulk(draws)  
print(bulk_ess)
tail_ess <- ess_tail(draws) 
print(tail_ess)

posterior_draws <- as_draws_df(EDSS_ordinal)
pp_check_density <- pp_check(EDSS_ordinal, type = "dens_overlay")
print(pp_check_density) 

mcmc_plot(EDSS_ordinal, type = "rhat")

#RNFL 
histogram_plot <- mcmc_plot(EDSS_ordinal, type = "hist", variable = "b_RNFL_V1_scaled")
histogram_plot <- histogram_plot +
  labs(
    x = "Posterior Estimates for Baseline RNFL Thickness",
    y = "Frequency")
print(histogram_plot)

trace_plot <- mcmc_plot(EDSS_ordinal, type = "trace", variable = "b_RNFL_V1_scaled")
trace_plot <- trace_plot +
  labs(
    x = "Posterior Draws",
    y = "Posterior Estimates for Baseline RNFL Thickness")
print(trace_plot)

#AGE  
histogram_plot <- mcmc_plot(EDSS_ordinal, type = "hist", variable = "b_Age_V1_scaled")
histogram_plot <- histogram_plot +
  labs(
    x = "Posterior Estimates for Age",
    y = "Frequency")
print(histogram_plot)

trace_plot <- mcmc_plot(EDSS_ordinal, type = "trace", variable = "b_Age_V1_scaled")
trace_plot <- trace_plot +
  labs(
    x = "Posterior Draws",
    y = "Posterior Estimates for Age")
print(trace_plot)

#_____________________________ MISSINGNESS ANALYSIS _______________________________

miss <- read_excel("Longitudina_Data_Manuscript3_out.xlsx")

#extract baseline/V1 predictors
variables <- c( "RNFL", "VG_Latency", "AS_Latency", "AS_Error", "AS_FEP", "VDI",
                "Grp_Coded", "Age", "Sex", "NART", 
               "SDMT_z", "PASAT_zc", "CVLT_z") 

visit1_data <- miss %>%
  filter(Visit == 0) %>%
  dplyr::select(ID, all_of(variables))

visit1_data <- visit1_data %>%
  rename_with(~ paste0(.x, "_V1"), -ID)

miss  <- miss  %>%
  left_join(visit1_data, by = "ID")

#subset longitudinal outcomes 
selected_vars_long <- c("CVLT", "SDMT", "PASAT", "EDSS", "BPF", "LF" )
subset_data <- miss[selected_vars_long]

#visualise missing data
vis_miss(subset_data) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  
    axis.text.x = element_text(angle = 60, hjust = 0, size = 12),  
    axis.text.y = element_text(size = 12),  
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 12, margin = margin(r = 20)), 
    legend.title = element_text(size = 12),  
    legend.text = element_text(size = 10)    
  ) +
  scale_y_continuous(
    limits = c(1, 351),  
    breaks = c(59, 176, 294),  
    labels = c("Baseline", "Follow-Up 1", "Follow-Up 2"), 
    name = "Individual Observations"  
  ) +
  scale_fill_manual(
    values = c("white", "grey"),  
    labels = c("", "Missing Data"),  
    name = ""  
  ) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = 115, alpha = 0.1, fill = "green") +  
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 116, ymax = 235, alpha = 0.1, fill = "blue") + 
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 236, ymax = 351, alpha = 0.1, fill = "red")  

## MISSING STUDY VISITS

#baseline predictors
wilcox.test(RNFL_V1 ~ Missing_Visit, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(VG_Latency_V1)) %>%  # Exclude missing Age_V1 values
  dplyr::group_by(Missing_Visit) %>%
  dplyr::summarize(
    Median = median(VG_Latency_V1, na.rm = TRUE),
    IQR = IQR(VG_Latency_V1, na.rm = TRUE))
print(summary_stats)

wilcox.test(AS_Latency_V1 ~ Missing_Visit, data = miss, na.action = na.omit)

wilcox.test(AS_Error_V1 ~ Missing_Visit, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(AS_Error_V1)) %>%  # Exclude missing Age_V1 values
  dplyr::group_by(Missing_Visit) %>%
  dplyr::summarize(
    Median = median(AS_Error_V1, na.rm = TRUE),
    IQR = IQR(AS_Error_V1, na.rm = TRUE))
print(summary_stats)

wilcox.test(AS_FEP_V1 ~ Missing_Visit, data = miss, na.action = na.omit)

wilcox.test(VDI_V1 ~ Missing_Visit, data = miss, na.action = na.omit)

wilcox.test(Age_V1 ~ Missing_Visit, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(Age_V1)) %>%  # Exclude missing Age_V1 values
  dplyr::group_by(Missing_Visit) %>%
  dplyr::summarize(
    Median = median(Age_V1, na.rm = TRUE),
    IQR = IQR(Age_V1, na.rm = TRUE))
print(summary_stats)

wilcox.test(NART_V1 ~ Missing_Visit, data = miss, na.action = na.omit)

table(miss$Grp_Coded_V1,miss$Missing_Visit )
contingency_table_grp <- table(miss$Missing_Visit, miss$Grp_Coded_V1)
chi_squared_test_grp <- chisq.test(contingency_table_grp)
print(chi_squared_test_grp)

table(miss$Sex_V1,miss$Missing_Visit )
contingency_table_sex <- table(miss$Missing_Visit, miss$Sex_V1)
chi_squared_test_sex <- chisq.test(contingency_table_sex)
print(chi_squared_test_sex)

#longitudinal outcomes 
wilcox.test(SDMT_z ~ Missing_Visit, data = miss, na.action = na.omit)

wilcox.test(CVLT_z ~ Missing_Visit, data = miss, na.action = na.omit)

wilcox.test(PASAT_zc ~ Missing_Visit, data = miss, na.action = na.omit)

wilcox.test(BPF ~ Missing_Visit, data = miss, na.action = na.omit)

wilcox.test(LF ~ Missing_Visit, data = miss, na.action = na.omit)

wilcox.test(EDSS~ Missing_Visit, data = miss, na.action = na.omit)

## MISSING NP DATA (ANALYSIS 2)

#baseline predictors 
wilcox.test(RNFL_V1 ~ Missing_NP, data = miss, na.action = na.omit)

wilcox.test(VG_Latency_V1 ~ Missing_NP, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(VG_Latency_V1)) %>%  
  dplyr::group_by(Missing_NP) %>%
  dplyr::summarize(
    Median = median(VG_Latency_V1, na.rm = TRUE),
    IQR = IQR(VG_Latency_V1, na.rm = TRUE))
print(summary_stats)

wilcox.test(AS_Latency_V1 ~ Missing_NP, data = miss, na.action = na.omit)

wilcox.test(AS_Error_V1 ~ Missing_NP, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(AS_Error_V1)) %>%  
  dplyr::group_by(Missing_NP) %>%
  dplyr::summarize(
    Median = median(AS_Error_V1, na.rm = TRUE),
    IQR = IQR(AS_Error_V1, na.rm = TRUE))
print(summary_stats)

wilcox.test(AS_FEP_V1 ~ Missing_NP, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(AS_FEP_V1)) %>%  
  dplyr::group_by(Missing_NP) %>%
  dplyr::summarize(
    Median = median(AS_FEP_V1, na.rm = TRUE),
    IQR = IQR(AS_FEP_V1, na.rm = TRUE))
print(summary_stats)

wilcox.test(VDI_V1 ~ Missing_NP, data = miss, na.action = na.omit)

wilcox.test(Age_V1 ~ Missing_NP, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(Age_V1)) %>% 
  dplyr::group_by(Missing_NP) %>%
  dplyr::summarize(
    Median = median(Age_V1, na.rm = TRUE),
    IQR = IQR(Age_V1, na.rm = TRUE))
print(summary_stats)

wilcox.test(NART_V1 ~ Missing_NP, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(NART_V1)) %>%  
  dplyr::group_by(Missing_NP) %>%
  dplyr::summarize(
    Median = median(NART_V1, na.rm = TRUE),
    IQR = IQR(NART_V1, na.rm = TRUE))
print(summary_stats)

table(miss$Grp_Coded_V1,miss$Missing_NP )
contingency_table_grp <- table(miss$Missing_NP, miss$Grp_Coded_V1)
chi_squared_test_grp <- chisq.test(contingency_table_grp)
print(chi_squared_test_grp)

table(miss$Sex_V1,miss$Missing_NP )
contingency_table_sex <- table(miss$Missing_NP, miss$Sex_V1)
chi_squared_test_sex <- chisq.test(contingency_table_sex)
print(chi_squared_test_sex)

#longitudinal outcome 
wilcox.test(SDMT_z ~ Missing_NP, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(SDMT_z)) %>%  
  dplyr::group_by(Missing_NP) %>%
  dplyr::summarize(
    Median = median(SDMT_z, na.rm = TRUE),
    IQR = IQR(SDMT_z, na.rm = TRUE))
print(summary_stats)

wilcox.test(CVLT_z ~ Missing_NP, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(CVLT_z)) %>%  
  dplyr::group_by(Missing_NP) %>%
  dplyr::summarize(
    Median = median(CVLT_z, na.rm = TRUE),
    IQR = IQR(CVLT_z, na.rm = TRUE))
print(summary_stats)

wilcox.test(PASAT_zc ~ Missing_NP, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(PASAT_zc)) %>%  
  dplyr::group_by(Missing_NP) %>%
  dplyr::summarize(
    Median = median(PASAT_zc, na.rm = TRUE),
    IQR = IQR(PASAT_zc, na.rm = TRUE))
print(summary_stats)

wilcox.test(BPF ~ Missing_NP, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(BPF)) %>%  
  dplyr::group_by(Missing_NP) %>%
  dplyr::summarize(
    Median = median(BPF, na.rm = TRUE),
    IQR = IQR(BPF, na.rm = TRUE))
print(summary_stats)

wilcox.test(LF ~ Missing_NP, data = miss, na.action = na.omit)

wilcox.test(EDSS~ Missing_NP, data = miss, na.action = na.omit)

## MISSING MRI DATA (ANALYSIS 3a)

#baseline predictors
wilcox.test(RNFL_V1 ~ Missing_MRI, data = miss, na.action = na.omit)

wilcox.test(VG_Latency_V1 ~ Missing_MRI, data = miss, na.action = na.omit)

wilcox.test(AS_Latency_V1 ~ Missing_MRI, data = miss, na.action = na.omit)

wilcox.test(AS_Error_V1 ~ Missing_MRI, data = miss, na.action = na.omit)

wilcox.test(AS_FEP_V1 ~ Missing_MRI, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(AS_FEP_V1)) %>%  
  dplyr::group_by(Missing_MRI) %>%
  dplyr::summarize(
    Median = median(AS_FEP_V1, na.rm = TRUE),
    IQR = IQR(AS_FEP_V1, na.rm = TRUE))
print(summary_stats)

wilcox.test(VDI_V1 ~ Missing_MRI, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(VDI_V1)) %>%  
  dplyr::group_by(Missing_MRI) %>%
  dplyr::summarize(
    Median = median(VDI_V1, na.rm = TRUE),
    IQR = IQR(VDI_V1, na.rm = TRUE))
print(summary_stats)

wilcox.test(Age_V1 ~ Missing_MRI, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(Age_V1)) %>%  
  dplyr::group_by(Missing_MRI) %>%
  dplyr::summarize(
    Median = median(Age_V1, na.rm = TRUE),
    IQR = IQR(Age_V1, na.rm = TRUE))
print(summary_stats)

wilcox.test(NART_V1 ~ Missing_MRI, data = miss, na.action = na.omit)

table(miss$Grp_Coded_V1,miss$Missing_MRI)
contingency_table_grp <- table(miss$Missing_MRI, miss$Grp_Coded_V1)
chi_squared_test_grp <- chisq.test(contingency_table_grp)
print(chi_squared_test_grp)

table(miss$Sex_V1,miss$Missing_MRI )
contingency_table_sex <- table(miss$Missing_MRI, miss$Sex_V1)
chi_squared_test_sex <- chisq.test(contingency_table_sex)
print(chi_squared_test_sex)

#longitudinal outcomes 
wilcox.test(SDMT_z ~ Missing_MRI, data = miss, na.action = na.omit)

wilcox.test(CVLT_z ~ Missing_MRI, data = miss, na.action = na.omit)

wilcox.test(PASAT_zc ~ Missing_MRI, data = miss, na.action = na.omit)

wilcox.test(EDSS~ Missing_MRI, data = miss, na.action = na.omit)


## MISSING EDSS DATA (ANALYSIS 3b)

#baseline predictors
wilcox.test(RNFL_V1 ~ Missing_EDSS, data = miss, na.action = na.omit)

wilcox.test(VG_Latency_V1 ~ Missing_EDSS, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(VG_Latency_V1)) %>%  
  dplyr::group_by(Missing_EDSS) %>%
  dplyr::summarize(
    Median = median(VG_Latency_V1, na.rm = TRUE),
    IQR = IQR(VG_Latency_V1, na.rm = TRUE))
print(summary_stats)

wilcox.test(AS_Latency_V1 ~ Missing_EDSS, data = miss, na.action = na.omit)

wilcox.test(AS_Error_V1 ~ Missing_EDSS, data = miss, na.action = na.omit)

wilcox.test(AS_FEP_V1 ~ Missing_EDSS, data = miss, na.action = na.omit)

wilcox.test(VDI_V1 ~ Missing_EDSS, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(VDI_V1)) %>%  
  dplyr::group_by(Missing_EDSS) %>%
  dplyr::summarize(
    Median = median(VDI_V1, na.rm = TRUE),
    IQR = IQR(VDI_V1, na.rm = TRUE))
print(summary_stats)

wilcox.test(Age_V1 ~ Missing_EDSS, data = miss, na.action = na.omit)

wilcox.test(NART_V1 ~ Missing_EDSS, data = miss, na.action = na.omit)

table(miss$Grp_Coded_V1,miss$Missing_EDSS )
contingency_table_grp <- table(miss$Missing_EDSS, miss$Grp_Coded_V1)
chi_squared_test_grp <- chisq.test(contingency_table_grp)
print(chi_squared_test_grp)

table(miss$Sex_V1,miss$Missing_EDSS )
contingency_table_sex <- table(miss$Missing_EDSS, miss$Sex_V1)
chi_squared_test_sex <- chisq.test(contingency_table_sex)
print(chi_squared_test_sex)

#longitudinal outcomes 
wilcox.test(SDMT_z ~ Missing_EDSS, data = miss, na.action = na.omit)

wilcox.test(CVLT_z ~ Missing_EDSS, data = miss, na.action = na.omit)

wilcox.test(PASAT_zc ~ Missing_EDSS, data = miss, na.action = na.omit)

wilcox.test(EDSS~ Missing_EDSS, data = miss, na.action = na.omit)

wilcox.test(BPF ~ Missing_EDSS, data = miss, na.action = na.omit)
summary_stats <- miss %>%
  filter(!is.na(BPF)) %>%  
  dplyr::group_by(Missing_EDSS) %>%
  dplyr::summarize(
    Median = median(BPF, na.rm = TRUE),
    IQR = IQR(BPF, na.rm = TRUE))
print(summary_stats)

wilcox.test(LF ~ Missing_EDSS, data = miss, na.action = na.omit)

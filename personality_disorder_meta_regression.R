#####################################################################
#### Estimation of personality disorder prevalence and mortlaity ####
#####################################################################

library(data.table)
library(openxlsx)
library(msm)
library(metafor)
library(ggplot2)
source("/FILEPATH/get_ids.R")
source("/FILEPATH/get_location_metadata.R")
source("/FILEPATH/custom_functions.R")
source("/FILEPATH/get_crosswalk_version.R")
source("/FILEPATH/get_population.R")
library(reticulate)
library(msm)


reticulate::use_python("/FILEPATH/python")
mr <- import("mrtool")

# Load age-metadata -------------------------------------------------------
age_ids <- get_ids('age_group')[age_group_id %in%  c(2:3, 388:389, 238, 34,  6:20, 30:32, 235),]
suppressWarnings(age_ids[, `:=` (age_start = as.numeric(unlist(strsplit(age_group_name, " "))[1]), age_end = as.numeric(unlist(strsplit(age_group_name, " "))[3])), by = "age_group_id"])
age_ids[age_group_id %in% c(2, 3, 388, 389), `:=` (age_start = 0, age_end = 0)]
age_ids[age_group_id %in% c(238), `:=` (age_start = 1, age_end = 1)]
age_ids[age_start == 95, age_end := 99]
age_ids[, mid_age := (age_start+age_end)/2]

locations <- get_location_metadata(9, release_id = 9)
sex <- get_ids ("sex")
release <- get_ids("release")


# Load in and prep data ---------------------------------------------------

dataset <- data.table(read.xlsx("/FILEPATH/systematic_review_extractions.xlsx", sheet = "extraction"))
dataset <- dataset[mean != 0,]
dataset[, sex_id := ifelse(sex ==  "Male", 1, ifelse(sex == "Female", 2, 3))]
dataset[age_end == 99, age_end := 74]
deltamethod <- msm::deltamethod

#transform prevalence data
dataset[measure == "prevalence", `:=`(mean_logit = logit(mean), lower_logit = logit(lower), upper_logit = logit(upper))]
dataset[measure == "prevalence" & lower > 0, `:=`(se_logit = (upper_logit - lower_logit) / 3.92)]
dataset[measure == "prevalence" & is.na(se_logit), `:=` (se_logit = deltamethod(~log(x1/(1-x1)), mean, standard_error^2)), by = c("row_id")]

#transform mortality data
dataset[measure %in% c("mtstandard", "relrisk"), `:=`(mean_log = log(mean), lower_log = log(lower), upper_log = log(upper))]
dataset[measure %in% c("mtstandard", "relrisk") & lower > 0, `:=`(se_log = (upper_log - lower_log) / 3.92)]
dataset[measure %in% c("mtstandard", "relrisk") & is.na(se_log), `:=` (se_log = deltamethod(~log(x1), mean, standard_error^2)), by = c("row_id")]


#if you want to model with age as a covariate for prevalence
dataset[measure == 'prevalence', age := ifelse(
  (age_end - age_start) < 60,
  (age_start + age_end) / 2,
  ifelse(!is.na(age_mean), age_mean, mean(age_mean, na.rm = TRUE))
)]

mean_age_prev <- dataset[measure == 'prevalence', mean(age)]
dataset[measure == 'prevalence', m_age_prev := age - mean_age_prev]

#if you want to model with year as a covariate for prevalence
dataset[, year := (year_start + year_end)/2]

mean_year_prev <- dataset[measure == 'prevalence', mean(year)]
dataset[measure == 'prevalence', m_year_prev := year - 2024]

#if you want to model with age/sex/year as a covariate for mortality
dataset[measure %in% c('relrisk', 'mtstandard'), age := ifelse(
  (age_end - age_start) < 60,
  (age_start + age_end) / 2,
  ifelse(!is.na(age_mean), age_mean, mean(age_mean, na.rm = TRUE))
)]

# Compute mean age for mortality rows
mean_age_mort <- dataset[measure %in% c('relrisk', 'mtstandard'), mean(age)]

# Center age for mortality
dataset[measure %in% c('relrisk', 'mtstandard'), m_age_mort := age - mean_age_mort]

mean_year_mort <- dataset[measure %in% c('relrisk', 'mtstandard'), mean(year)]
dataset[measure %in% c('relrisk', 'mtstandard'), m_year_mort := year - 2024]

dataset[sex == "Male", percent_female := 0]
dataset[sex == "Female", percent_female := 1]
dataset[sex == "Both" & is.na(percent_female), percent_female := 0.5]
dataset[,m_percent_female:= percent_female - 0.5]

dataset <- merge(dataset, locations[,.(location_id, region_id, region_name, super_region_id, super_region_name, path_to_top_parent)], all.x=TRUE, by = 'location_id')

# Manually set values for location_id 95
dataset[location_id == 95, `:=` (region_id = 73, region_name = "Western Europe", super_region_id = 64, super_region_name = "High-income")]
dataset[location_id == 4624, `:=` (region_id = 73, region_name = "Western Europe", super_region_id = 64, super_region_name = "High-income")]
dataset[location_id == 4626, `:=` (region_id = 73, region_name = "Western Europe", super_region_id = 64, super_region_name = "High-income")]
dataset[location_id == 434, `:=` (region_id = 73, region_name = "Western Europe", super_region_id = 64, super_region_name = "High-income")]

# fill in location specific covariates
dataset[, `:=`(
  cv_high_income = ifelse(region_name == "High-income", 0, 1),
  cv_north_america = ifelse(region_name == "High-income North America", 1, 0),
  cv_low_middle_income = ifelse(super_region_name == "High-income", 0, 1),
  cv_australasia = ifelse(region_name == "Australasia", 1, 0),
  cv_south_asia = ifelse(region_name == "South Asia", 1, 0),
  cv_western_europe = ifelse(region_name == "Western Europe", 1, 0),
  cv_south_and_east_asia = ifelse(region_name %in% c("South Asia", "East Asia"), 1, 0),
  cv_high_income_other = ifelse(region_name %in% c("Australasia", "Southern Latin America"), 1, 0),
  cv_low_income_other = ifelse(region_name %in% c("Southern Sub-Saharan Africa", "Western Sub-Saharan Africa", "North Africa and Middle East", "Central Latin America"), 1, 0)
)]
dataset[, `:=`(cv_high_income_other = ifelse(super_region_name == "High-income" & cv_north_america == 0 & cv_western_europe == 0, 1, 0))]


# Develop meta-regression model for prevalence:-------------------------------------------

#trimming <- 0.95 # for sensitivity testing
trimming <- 1

# Develop meta-regression model for mortality:
#For prevalence: the following covariates were tested while controlling for ASPD and BPD via their respecitve covariates, and they were found to have a p-value >0.1 and therefore not included in the meta-regression: m_year_prev, cv_pd_lay, cv_aspd_lay, cv_wmhs, cv_2_stage, cv_western_europe
potential_cvs <- c("m_percent_female","m_age_prev", "cv_icd", "cv_pd_nonclin", "cv_aspd_nonclin", "cv_bpd_nonclin","cv_antisocial_only", "cv_borderline_only", "cv_high_income_other", "cv_low_middle_income")

mr_dataset <- mr$MRData()
mr_dataset$load_df(
  data = dataset[measure == "prevalence",], 
  col_obs = "mean_logit", col_obs_se = "se_logit",
  col_covs = as.list(c(potential_cvs, "row_id")), col_study_id = "study_id")


## Step 1
selected_covs <- potential_cvs
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas1 <- betas[,.(model = "Model 1", cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]
betas1 # worst is cv_high_income_other
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas1[, `:=` (aic = aic_current, aic_diff = NA, bic = bic_current, bic_diff = NA)]
aic_previous <- aic_current
bic_previous <- bic_current
model1 <- model


## Step 2
selected_covs <- selected_covs[selected_covs != "cv_high_income_other"]
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas2 <- betas[,.(model = 'Model 2', cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]
betas2 # worst is cv_aspd_nonclin
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas2[, `:=` (aic = aic_current, aic_diff = aic_current - aic_previous, bic = bic_current, bic_diff = bic_current - bic_previous)]
aic_previous <- aic_current
bic_previous <- bic_current
model2 <- model

## Step 3
selected_covs <- selected_covs[selected_covs != "cv_aspd_nonclin"]
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas3 <- betas[,.(model = 'Model 3', cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]
betas # all significant but worst is cv_low_middle_income
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas3[, `:=` (aic = aic_current, aic_diff = aic_current - aic_previous, bic = bic_current, bic_diff = bic_current - bic_previous)]
aic_previous <- aic_current
bic_previous <- bic_current
model3 <- model

# Step 4
selected_covs <- selected_covs[selected_covs != "cv_low_middle_income"]
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas4 <- betas[,.(model = "Model 4", cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = p)]
betas4 # worst is cv_2_stage
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas4[, `:=` (aic = aic_current, aic_diff = aic_current - aic_previous, bic = bic_current, bic_diff = bic_current - bic_previous)]
aic_previous <- aic_current
bic_previous <- bic_current
model4 <- model


betas_table <- rbind(betas1, betas2, betas3, betas4)

betas_table <- betas_table[,.(Model = model, Covariate = cov, Coefficient = beta, `Lower UI` = lower, `Upper UI` = upper, p = round(p, 4), AIC = aic, BIC = bic)]
betas_table[Covariate == 'intercept', Covariate := 'Intercept']
betas_table[Covariate == 'm_percent_female', Covariate := 'Percent female']
betas_table[Covariate == 'm_age_prev', Covariate := 'Age']
betas_table[Covariate == 'cv_icd', Covariate := 'ICD']
betas_table[Covariate == 'cv_pd_nonclin', Covariate := 'Mental health professional administered tool - PD']
betas_table[Covariate == 'cv_aspd_nonclin', Covariate := 'Mental health professional administered tool - ASPD']
betas_table[Covariate == 'cv_bpd_nonclin', Covariate := 'Mental health professional administered tool - BPD']
betas_table[Covariate == 'cv_antisocial_only', Covariate := 'Antisocial personality disorder']
betas_table[Covariate == 'cv_borderline_only', Covariate := 'Borderline personality disorder']
betas_table[Covariate == 'cv_high_income_other', Covariate := 'Australasia and Southern Latin America']
betas_table[Covariate == 'cv_low_middle_income', Covariate := 'Low and middle income countries']
betas_table

# Summary of results ------------------------------------------------------

#model 1 had the best AIC
model <- model1
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]

## Funnel plot & Egger test
funnel_plot_mrbrt_bydata(model)
dev.off() # Turn the PDF device off

funnel_plot_mrbrt_bystudy(model)
dev.off() # Turn the PDF device off

egger_mr_brt_pval(model)


## Summary prevalence ##

summary_prev<- rbind(data.table(intercept = 1, m_age_prev = 0, m_percent_female = 0, cv_high_income_other = 0, cv_low_middle_income = 0, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 0, m_percent_female = 0.5, cv_high_income_other = 0, cv_low_middle_income = 0, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 0, m_percent_female = -0.5, cv_high_income_other = 0, cv_low_middle_income = 0, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 0, m_percent_female = 0, cv_high_income_other = 0, cv_low_middle_income = 1, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 0, m_percent_female = 0.5, cv_high_income_other = 0, cv_low_middle_income = 1, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 0, m_percent_female = -0.5, cv_high_income_other = 0, cv_low_middle_income = 1, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     
                     data.table(intercept = 1, m_age_prev = 20 - mean_age_prev, m_percent_female = 0, cv_high_income_other = 0, cv_low_middle_income = 0, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 20 - mean_age_prev, m_percent_female = 0.5, cv_high_income_other = 0, cv_low_middle_income = 0, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 20 - mean_age_prev, m_percent_female = -0.5, cv_high_income_other = 0, cv_low_middle_income = 0, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 20 - mean_age_prev, m_percent_female = 0, cv_high_income_other = 0, cv_low_middle_income = 1, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 20 - mean_age_prev, m_percent_female = 0.5, cv_high_income_other = 0, cv_low_middle_income = 1, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 20 - mean_age_prev, m_percent_female = -0.5, cv_high_income_other = 0, cv_low_middle_income = 1, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     
                     data.table(intercept = 1, m_age_prev = 40 - mean_age_prev, m_percent_female = 0, cv_high_income_other = 0, cv_low_middle_income = 0, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 40 - mean_age_prev, m_percent_female = 0.5, cv_high_income_other = 0, cv_low_middle_income = 0, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 40 - mean_age_prev, m_percent_female = -0.5, cv_high_income_other = 0, cv_low_middle_income = 0, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 40 - mean_age_prev, m_percent_female = 0, cv_high_income_other = 0, cv_low_middle_income = 1, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 40 - mean_age_prev, m_percent_female = 0.5, cv_high_income_other = 0, cv_low_middle_income = 1, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 40 - mean_age_prev, m_percent_female = -0.5, cv_high_income_other = 0, cv_low_middle_income = 1, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     
                     data.table(intercept = 1, m_age_prev = 60 - mean_age_prev, m_percent_female = 0, cv_high_income_other = 0, cv_low_middle_income = 0, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 60 - mean_age_prev, m_percent_female = 0.5, cv_high_income_other = 0, cv_low_middle_income = 0, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 60 - mean_age_prev, m_percent_female = -0.5, cv_high_income_other = 0, cv_low_middle_income = 0, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 60 - mean_age_prev, m_percent_female = 0, cv_high_income_other = 0, cv_low_middle_income = 1, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 60 - mean_age_prev, m_percent_female = 0.5, cv_high_income_other = 0, cv_low_middle_income = 1, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 60 - mean_age_prev, m_percent_female = -0.5, cv_high_income_other = 0, cv_low_middle_income = 1, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     
                     data.table(intercept = 1, m_age_prev = 80 - mean_age_prev, m_percent_female = 0, cv_high_income_other = 0, cv_low_middle_income = 0, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 80 - mean_age_prev, m_percent_female = 0.5, cv_high_income_other = 0, cv_low_middle_income = 0, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 80 - mean_age_prev, m_percent_female = -0.5, cv_high_income_other = 0, cv_low_middle_income = 0, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 80 - mean_age_prev, m_percent_female = 0, cv_high_income_other = 0, cv_low_middle_income = 1, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 80 - mean_age_prev, m_percent_female = 0.5, cv_high_income_other = 0, cv_low_middle_income = 1, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0),
                     data.table(intercept = 1, m_age_prev = 80 - mean_age_prev, m_percent_female = -0.5, cv_high_income_other = 0, cv_low_middle_income = 1, cv_icd = 0, cv_pd_nonclin =0, cv_aspd_nonclin =0, cv_bpd_nonclin = 0, cv_antisocial_only = 0, cv_borderline_only = 0))

#prevalence prediction
predict_data <- mr$MRData()
predict_data$load_df(data = summary_prev, col_covs=as.list(model$cov_names))
beta_samples <- mr$core$other_sampling$sample_simple_lme_beta(500L, model)
gamma_outer_samples <- matrix(rep(model$gamma_soln, each = 500L), nrow = 500L)

draws <- model$create_draws(predict_data,
                            beta_samples = beta_samples,
                            gamma_samples = gamma_outer_samples,
                            random_study = F)

summary_prev$pred <- model$predict(data = predict_data)
summary_prev$pred_lo <- apply(draws, 1, function(x) quantile(x, 0.025))
summary_prev$pred_hi <- apply(draws, 1, function(x) quantile(x, 0.975))


summary_prev[, sex := ifelse(m_percent_female == 0, "Both", ifelse(m_percent_female < 0, "Male", "female"))]
summary_prev[, region := ifelse(cv_low_middle_income == 1, "Low- and middle-income","High-income")]
summary_prev <- summary_prev[, .(sex, age = round(m_age_prev + mean_age_prev, 1), prev = round(100*rlogit(pred), 1), lower = round(100*rlogit(pred_lo), 1), upper = round(100*rlogit(pred_hi), 1))]

summary_prev

# Develop meta-regression model for mortality:-------------------------------------------

#trimming <- 0.95 # for sensitivity testing
trimming <- 1

#for mortality: the following covariates were tested in a univariate model, and found to have a p-value >0.1 and therefore were not included in the meta-regression: m_year_mort, cv_pd_broad
potential_cvs <- c("m_percent_female", "m_age_mort","cv_inpatient",  "cv_outpatient")

mr_dataset <- mr$MRData()
mr_dataset$load_df(
  data = dataset[measure %in% c('relrisk', 'mtstandard'),], 
  col_obs = "mean_log", col_obs_se = "se_log",
  col_covs = as.list(c(potential_cvs, "row_id")), col_study_id = "study_id")

## Step 1
selected_covs <- potential_cvs
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas1 <- betas[,.(model = "Model 1", cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]
betas1 # all significant but worst is cv_outpatient
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas1[, `:=` (aic = aic_current, aic_diff = NA, bic = bic_current, bic_diff = NA)]
aic_previous <- aic_current
bic_previous <- bic_current
model1 <- model


## Step 2
selected_covs <- selected_covs[selected_covs != "cv_outpatient"]
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas2 <- betas[,.(model = 'Model 2', cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]
betas2 # all significant
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas2[, `:=` (aic = aic_current, aic_diff = aic_current - aic_previous, bic = bic_current, bic_diff = bic_current - bic_previous)]
aic_previous <- aic_current
bic_previous <- bic_current
model2 <- model

betas_table <- rbind(betas1, betas2)

betas_table <- betas_table[,.(Model = model, Covariate = cov, Coefficient = beta, `Lower UI` = lower, `Upper UI` = upper, p = round(p, 4), AIC = aic, BIC = bic)]
betas_table[Covariate == 'intercept', Covariate := 'Intercept']
betas_table[Covariate == 'm_percent_female', Covariate := 'Percent female']
betas_table[Covariate == 'm_age_mort', Covariate := 'Age']
betas_table[Covariate == 'cv_inpatient', Covariate := 'Inpatient only']
betas_table[Covariate == 'cv_outpatient', Covariate := 'Outpatient only']
betas_table


# Summary of results ------------------------------------------------------

#model 1 had the best AIC
model <- model1
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]

## Funnel plot & Egger test
funnel_plot_mrbrt_bydata(model)
dev.off() # Turn the PDF device off

funnel_plot_mrbrt_bystudy(model)
dev.off() # Turn the PDF device off

egger_mr_brt_pval(model)


## Summary mortality ##

summary_mort<- rbind(data.table(intercept = 1, m_age_mort = 0, cv_inpatient = 0, cv_outpatient = 0, m_percent_female = 0.5), 
                     data.table(intercept = 1, m_age_mort = 0, cv_inpatient = 0, cv_outpatient = 0, m_percent_female = -0.5),
                     data.table(intercept = 1, m_age_mort = 0, cv_inpatient = 1, cv_outpatient = 0, m_percent_female = 0),
                     data.table(intercept = 1, m_age_mort = 0, cv_inpatient = 1, cv_outpatient = 0, m_percent_female = 0.5), 
                     data.table(intercept = 1, m_age_mort = 0, cv_inpatient = 1, cv_outpatient = 0, m_percent_female = -0.5),
                     data.table(intercept = 1, m_age_mort = 0, cv_inpatient = 0, cv_outpatient = 1, m_percent_female = 0),
                     data.table(intercept = 1, m_age_mort = 0, cv_inpatient = 0, cv_outpatient = 1, m_percent_female = 0.5), 
                     data.table(intercept = 1, m_age_mort = 0, cv_inpatient = 0, cv_outpatient = 1, m_percent_female = -0.5),
                     
                     data.table(intercept = 1, m_age_mort = 20 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 0, m_percent_female = 0),
                     data.table(intercept = 1, m_age_mort = 20 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 0, m_percent_female = 0.5), 
                     data.table(intercept = 1, m_age_mort = 20 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 0, m_percent_female = -0.5),
                     data.table(intercept = 1, m_age_mort = 20 - mean_age_mort, cv_inpatient = 1, cv_outpatient = 0, m_percent_female = 0),
                     data.table(intercept = 1, m_age_mort = 20 - mean_age_mort, cv_inpatient = 1, cv_outpatient = 0, m_percent_female = 0.5), 
                     data.table(intercept = 1, m_age_mort = 20 - mean_age_mort, cv_inpatient = 1, cv_outpatient = 0, m_percent_female = -0.5),
                     data.table(intercept = 1, m_age_mort = 20 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 1, m_percent_female = 0),
                     data.table(intercept = 1, m_age_mort = 20 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 1, m_percent_female = 0.5), 
                     data.table(intercept = 1, m_age_mort = 20 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 1, m_percent_female = -0.5),
                     
                     data.table(intercept = 1, m_age_mort = 40 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 0, m_percent_female = 0),
                     data.table(intercept = 1, m_age_mort = 40 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 0, m_percent_female = 0.5), 
                     data.table(intercept = 1, m_age_mort = 40 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 0, m_percent_female = -0.5),
                     data.table(intercept = 1, m_age_mort = 40 - mean_age_mort, cv_inpatient = 1, cv_outpatient = 0, m_percent_female = 0),
                     data.table(intercept = 1, m_age_mort = 40 - mean_age_mort, cv_inpatient = 1, cv_outpatient = 0, m_percent_female = 0.5), 
                     data.table(intercept = 1, m_age_mort = 40 - mean_age_mort, cv_inpatient = 1, cv_outpatient = 0, m_percent_female = -0.5),
                     data.table(intercept = 1, m_age_mort = 40 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 1, m_percent_female = 0),
                     data.table(intercept = 1, m_age_mort = 40 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 1, m_percent_female = 0.5), 
                     data.table(intercept = 1, m_age_mort = 40 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 1, m_percent_female = -0.5),
                     
                     data.table(intercept = 1, m_age_mort = 60 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 0, m_percent_female = 0),
                     data.table(intercept = 1, m_age_mort = 60 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 0, m_percent_female = 0.5), 
                     data.table(intercept = 1, m_age_mort = 60 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 0, m_percent_female = -0.5),
                     data.table(intercept = 1, m_age_mort = 60 - mean_age_mort, cv_inpatient = 1, cv_outpatient = 0, m_percent_female = 0),
                     data.table(intercept = 1, m_age_mort = 60 - mean_age_mort, cv_inpatient = 1, cv_outpatient = 0, m_percent_female = 0.5), 
                     data.table(intercept = 1, m_age_mort = 60 - mean_age_mort, cv_inpatient = 1, cv_outpatient = 0, m_percent_female = -0.5),
                     data.table(intercept = 1, m_age_mort = 60 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 1, m_percent_female = 0),
                     data.table(intercept = 1, m_age_mort = 60 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 1, m_percent_female = 0.5), 
                     data.table(intercept = 1, m_age_mort = 60 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 1, m_percent_female = -0.5),
                     
                     data.table(intercept = 1, m_age_mort = 80 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 0, m_percent_female = 0),
                     data.table(intercept = 1, m_age_mort = 80 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 0, m_percent_female = 0.5), 
                     data.table(intercept = 1, m_age_mort = 80 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 0, m_percent_female = -0.5),
                     data.table(intercept = 1, m_age_mort = 80 - mean_age_mort, cv_inpatient = 1, cv_outpatient = 0, m_percent_female = 0),
                     data.table(intercept = 1, m_age_mort = 80 - mean_age_mort, cv_inpatient = 1, cv_outpatient = 0, m_percent_female = 0.5), 
                     data.table(intercept = 1, m_age_mort = 80 - mean_age_mort, cv_inpatient = 1, cv_outpatient = 0, m_percent_female = -0.5),
                     data.table(intercept = 1, m_age_mort = 80 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 1, m_percent_female = 0),
                     data.table(intercept = 1, m_age_mort = 80 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 1, m_percent_female = 0.5), 
                     data.table(intercept = 1, m_age_mort = 80 - mean_age_mort, cv_inpatient = 0, cv_outpatient = 1, m_percent_female = -0.5))


#mortality prediction
predict_data <- mr$MRData()
predict_data$load_df(data = summary_mort, col_covs=as.list(model$cov_names))
beta_samples <- mr$core$other_sampling$sample_simple_lme_beta(500L, model)
gamma_outer_samples <- matrix(rep(model$gamma_soln, each = 500L), nrow = 500L)

draws <- model$create_draws(predict_data,
                            beta_samples = beta_samples,
                            gamma_samples = gamma_outer_samples,
                            random_study = F)

summary_mort$pred <- model$predict(data = predict_data)
summary_mort$pred_lo <- apply(draws, 1, function(x) quantile(x, 0.025))
summary_mort$pred_hi <- apply(draws, 1, function(x) quantile(x, 0.975))

summary_mort[, sample := ifelse(cv_inpatient == 1, "Inpatient only", ifelse(cv_outpatient == 1, "Outpatient only", "In and out patient"))]

summary_mort[, sex := ifelse(m_percent_female == 0, "Both", ifelse(m_percent_female < 0, "Male", "female"))]

summary_mort <- summary_mort[, .(sex, sample, age = round(m_age_mort + mean_age_mort, 1), mort = round(exp(pred), 1), lower = round(exp(pred_lo), 1), upper = round(exp(pred_hi), 1))]

summary_mort


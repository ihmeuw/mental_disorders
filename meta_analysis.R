#########################################################
#### Estimation of RRs of ASD to self-harm mortality ####
#########################################################

library(data.table)
library(openxlsx)
library(msm)
library(metafor)
library(ggplot2)
source("/FILEPATH/get_draws.R")
source("/FILEPATH/interpolate.R")
source("/FILEPATH/get_population.R")
source("/FILEPATH/get_ids.R")
source("/FILEPATH/custom_functions.R")
library(reticulate)
reticulate::use_python("/FILEPATH/python")
mr <- import("mrtool")

cod_correction_version <- 393 
como_version <- 1471

# Load age-metadata -------------------------------------------------------
age_ids <- get_ids('age_group')[age_group_id %in%  c(2:3, 388:389, 238, 34,  6:20, 30:32, 235),]
suppressWarnings(age_ids[, `:=` (age_start = as.numeric(unlist(strsplit(age_group_name, " "))[1]), age_end = as.numeric(unlist(strsplit(age_group_name, " "))[3])), by = "age_group_id"])
age_ids[age_group_id %in% c(2, 3, 388, 389), `:=` (age_start = 0, age_end = 0)]
age_ids[age_group_id %in% c(238), `:=` (age_start = 1, age_end = 1)]
age_ids[age_start == 95, age_end := 99]
age_ids[, mid_age := (age_start+age_end)/2]

# Load in and prep data ---------------------------------------------------
dataset <- data.table(read.xlsx("/filepath/systematic_review_extractions.xlsx", sheet = "Extraction"))
dataset[, `:=` (Notes = NULL, row_id = seq_len(.N), cv_attempt = ifelse(outcome == "Suicide attempts", 1, 0))]
suppressWarnings(dataset[, mean := as.numeric(mean)])
dataset[, sex_id := ifelse(sex ==  "Male", 1, ifelse(sex == "Female", 2, 3))]

rr_se_escalc <- data.table(escalc(measure="RR", ai =outcome_asd_yes, bi=outcome_asd_no, ci = outcome_control_yes, di=outcome_control_no, data=dataset[is.na(lower) & (measure %in% c("RR", "HR")),]))
dataset <- merge(dataset, rr_se_escalc[,.(row_id, yi_rr = yi, vi_rr = vi)], by = 'row_id', all = T)
irr_se_escalc <- data.table(escalc(measure="IRR", x1i =outcome_asd_yes, t1i=sample_asd, x2i = outcome_control_yes, t2i=sample_control, data=dataset[is.na(lower) & measure == "IRR",]))
dataset <- merge(dataset, irr_se_escalc[,.(row_id, yi_irr = yi, vi_irr = vi)], by = 'row_id', all = T)

dataset[measure == "RR" & !is.na(vi_rr), `:=` (log_rr = yi_rr, log_rr_se = sqrt(vi_rr))]
dataset[measure == "HR" & !is.na(vi_rr), `:=` (log_rr = yi_irr, log_rr_se = sqrt(vi_rr))]
dataset[measure == "IRR" & !is.na(vi_irr), `:=` (log_rr = yi_irr, log_rr_se = sqrt(vi_irr))]
dataset[, c("yi_rr", "vi_rr", "yi_irr", "vi_irr") := NULL]

dataset[is.na(lower), `:=` (lower = exp(log_rr - 1.96*log_rr_se), upper = exp(log_rr + 1.96*log_rr_se))]

# Cybulski et al., 2022 missing data on % female, and so sourced from Table 1 of following study which uses same data source (Clinical Practice Research Datalink)
# by estimating  number of female cases between 2003 and 2018 and diving this by the total cases between 2003 and 2018
# Russell, G., Stapley, S., Newlove-Delgado, T., Salmon, A., White, R., Warren, F., Pearson, A. and Ford, T. (2022), Time trends in autism diagnosis over 20 years: a UK population-based cohort study. J Child Psychol Psychiatr, 63: 674-682. https://doi-org.ezproxy.library.uq.edu.au/10.1111/jcpp.13505
dataset[figure_citation =="Cybulski et al., 2022", `:=` (percent_female_asd = 0.194760953)]

# Convert ORs to RRs ------------------------------------------------------
dataset[measure == "OR" & is.na(outcome_control_no), outcome_ratio_control := (mean * outcome_asd_no) / outcome_asd_yes]
dataset[measure == "OR" & is.na(outcome_control_no), `:=` (outcome_control_yes = sample_control / (1 + (outcome_ratio_control)))]
dataset[measure == "OR" & is.na(outcome_control_no), `:=` (outcome_control_no = sample_control - outcome_control_yes)]
dataset[measure == "OR", `:=` (mean = or_2_rr(mean, (outcome_asd_yes + outcome_control_yes) / (sample_asd + sample_control), sample_asd / (sample_asd + sample_control)),
                               lower = or_2_rr(lower, (outcome_asd_yes + outcome_control_yes) / (sample_asd + sample_control), sample_asd / (sample_asd + sample_control)),
                               upper = or_2_rr(upper, (outcome_asd_yes + outcome_control_yes) / (sample_asd + sample_control), sample_asd / (sample_asd + sample_control)),
                               measure = "RR"), by = "row_id"]

# Conduct remaining data transformations and age-sex-ID splitting -------------------

dataset[!is.na(lower) & is.na(log_rr_se), `:=` (log_rr = log(mean), log_rr_se = (log(upper) - log(lower))/(qnorm(0.975,0,1)*2))]
dataset[is.na(lower) & is.na(log_rr_se), `:=` (log_rr = log(mean), log_rr_se= deltamethod(~log(x1), mean, se^2)), by = c("row_id")]

## Age-sex-ID split Kõlves et al. 2021 study
dataset[, `:=` (min_age_start = min(age_start), max_age_end = max(age_end)), by = "citation"]
dataset[, widest_age := ifelse(age_start == min_age_start & age_end == max_age_end, 1, 0)]

for(attempt in c(0, 1)){
  dataset_kolves_total <- dataset[figure_citation == "Kõlves et al. 2021" & sex == "Both" & widest_age == 1 & id == "Both" & cv_attempt == attempt, .(log_rr, log_rr_se)]
  dataset_kolves_male <- dataset[figure_citation == "Kõlves et al. 2021" & sex == "Male" & widest_age == 1 & id == "Both" & cv_attempt == attempt, .(log_rr, log_rr_se)]
  dataset_kolves_female <- dataset[figure_citation == "Kõlves et al. 2021" & sex == "Female" & widest_age == 1 & id == "Both" & cv_attempt == attempt, .(log_rr, log_rr_se)]
  dataset_kolves_idpres <- dataset[figure_citation == "Kõlves et al. 2021" & sex == "Both" & widest_age == 1 & percent_id_asd == 1 & cv_attempt == attempt, .(log_rr, log_rr_se)]
  dataset_kolves_idabs <- dataset[figure_citation == "Kõlves et al. 2021" & sex == "Both" & widest_age == 1 & percent_id_asd == 0 & cv_attempt == attempt, .(log_rr, log_rr_se)]
  male_xwalk <- as.numeric(dataset_kolves_total$log_rr) - as.numeric(dataset_kolves_male$log_rr)
  male_xwalk_se <- sqrt(as.numeric(dataset_kolves_total$log_rr_se)^2 + as.numeric(dataset_kolves_male$log_rr_se)^2)
  female_xwalk <- as.numeric(dataset_kolves_total$log_rr) - as.numeric(dataset_kolves_female$log_rr)
  female_xwalk_se <- sqrt(as.numeric(dataset_kolves_total$log_rr_se)^2 + as.numeric(dataset_kolves_female$log_rr_se)^2)
  idpres_xwalk <-  as.numeric(dataset_kolves_total$log_rr) - as.numeric(dataset_kolves_idpres$log_rr)
  idpres_xwalk_se <- sqrt(as.numeric(dataset_kolves_total$log_rr_se)^2 + as.numeric(dataset_kolves_idpres$log_rr_se)^2)
  idabs_xwalk <-  as.numeric(dataset_kolves_total$log_rr) - as.numeric(dataset_kolves_idabs$log_rr)
  idabs_xwalk_se <- sqrt(as.numeric(dataset_kolves_total$log_rr_se)^2 + as.numeric(dataset_kolves_idabs$log_rr_se)^2)

  data_temp_male_idpres <- dataset[figure_citation == "Kõlves et al. 2021" & widest_age == 0 & cv_attempt == attempt,]
  data_temp_male_idpres[, `:=` (sex = "Male", percent_female_asd = 0, id = "ID present", percent_id_asd = 1, log_rr = log_rr - male_xwalk - idpres_xwalk, log_rr_se = sqrt(log_rr_se^2 + male_xwalk_se^2 + idpres_xwalk_se^2))]

  data_temp_female_idpres <- dataset[figure_citation == "Kõlves et al. 2021" & widest_age == 0 & cv_attempt == attempt,]
  data_temp_female_idpres[, `:=` (sex = "Female", percent_female_asd = 1, id = "ID present", percent_id_asd = 1, log_rr = log_rr - female_xwalk - idpres_xwalk, log_rr_se = sqrt(log_rr_se^2 + female_xwalk_se^2 + idpres_xwalk_se^2))]

  data_temp_male_idabs <- dataset[figure_citation == "Kõlves et al. 2021" & widest_age == 0 & cv_attempt == attempt,]
  data_temp_male_idabs[, `:=` (sex = "Male", percent_female_asd = 0, id = "ID absent", percent_id_asd = 0, log_rr = log_rr - male_xwalk - idabs_xwalk, log_rr_se = sqrt(log_rr_se^2 + male_xwalk_se^2 + idabs_xwalk_se^2))]

  data_temp_female_idabs <- dataset[figure_citation == "Kõlves et al. 2021" & widest_age == 0 & cv_attempt == attempt,]
  data_temp_female_idabs[, `:=` (sex = "Female", percent_female_asd = 1, id = "ID absent", percent_id_asd = 0, log_rr = log_rr - female_xwalk - idabs_xwalk, log_rr_se = sqrt(log_rr_se^2 + female_xwalk_se^2 + idabs_xwalk_se^2))]

  dataset <-  rbind(dataset[!(figure_citation == "Kõlves et al. 2021" & cv_attempt == attempt),], data_temp_male_idpres, data_temp_female_idpres, data_temp_male_idabs, data_temp_female_idabs)
}

## Sex-ID split Tsai et al. 2023 study
dataset_tsai_total <- dataset[figure_citation == "Tsai et al., 2023" & sex == "Both" & widest_age == 1 & id == "Both", .(log_rr, log_rr_se)]
dataset_tsai_male <- dataset[figure_citation == "Tsai et al., 2023" & sex == "Male" & widest_age == 1 & id == "Both", .(log_rr, log_rr_se)]
dataset_tsai_female <- dataset[figure_citation == "Tsai et al., 2023" & sex == "Female" & widest_age == 1 & id == "Both", .(log_rr, log_rr_se)]
dataset_tsai_idpres <- dataset[figure_citation == "Tsai et al., 2023" & sex == "Both" & widest_age == 1 & percent_id_asd == 1, .(log_rr, log_rr_se)]
dataset_tsai_idabs <- dataset[figure_citation == "Tsai et al., 2023" & sex == "Both" & widest_age == 1 & percent_id_asd == 0, .(log_rr, log_rr_se)]
male_xwalk <- as.numeric(dataset_tsai_total$log_rr) - as.numeric(dataset_kolves_male$log_rr)
male_xwalk_se <- sqrt(as.numeric(dataset_tsai_total$log_rr_se)^2 + as.numeric(dataset_kolves_male$log_rr_se)^2)
female_xwalk <- as.numeric(dataset_tsai_total$log_rr) - as.numeric(dataset_kolves_female$log_rr)
female_xwalk_se <- sqrt(as.numeric(dataset_tsai_total$log_rr_se)^2 + as.numeric(dataset_kolves_female$log_rr_se)^2)

data_temp_male <- dataset[figure_citation == "Tsai et al., 2023" & id != "Both",]
data_temp_male[, `:=` (sex = "Male", percent_female_asd = 0, log_rr = log_rr - male_xwalk, log_rr_se = sqrt(log_rr_se^2 + male_xwalk_se^2))]

data_temp_female <- dataset[figure_citation == "Tsai et al., 2023" & id != "Both",]
data_temp_female[, `:=` (sex = "Female", percent_female_asd = 1, log_rr = log_rr - female_xwalk, log_rr_se = sqrt(log_rr_se^2 + female_xwalk_se^2))]

dataset <-  rbind(dataset[!(figure_citation == "Tsai et al., 2023"),], data_temp_male, data_temp_female)

gbd2021_prop_id <- 0.356917979

dataset[, followup_years := ifelse(!is.na(followup_years_mean), followup_years_mean, followup_years_mid)]

dataset[, age := (age_start + age_end)/2]

dataset[!is.na(enrolment_age_mean_asd) & !is.na(followup_years_mean), age := enrolment_age_mean_asd + followup_years_mean/2]
dataset[!is.na(enrolment_age_mean_asd) & is.na(followup_years_mean) & !is.na(followup_years_mid), age := enrolment_age_mean_asd + followup_years_mid/2]
dataset[is.na(enrolment_age_mean_asd) & is.na(followup_years_mean) & is.na(followup_years_mid) & !is.na(followup_age_mean), age := followup_age_mean - followup_years_mid/2]

mean_age <- dataset[, mean(age)]
mean_age_log <- dataset[, mean(log(age))]

dataset[is.na(percent_id_asd), percent_id_asd := gbd2021_prop_id] # assume two studies missing ID data that their distribution matches global distribution from GBD 2020 -- have checked for alternative papers using these data sources reporting ID among ASD on 30/08/2022 with no luck

dataset[, `:=` (m_percent_female_asd = percent_female_asd - 0.5, m_percent_id_asd = percent_id_asd - gbd2021_prop_id, m_age = age - mean_age, m_age_log = log(age) - mean_age_log)]

dataset[, `:=` (int_attempt_female = cv_attempt * m_percent_female_asd, int_attempt_id = cv_attempt * m_percent_id_asd, int_attempt_age = cv_attempt * m_age, int_attempt_age_log = cv_attempt * m_age_log)]
dataset[, `:=` (int_suicide_female = (1-cv_attempt) * m_percent_female_asd, int_suicide_id = (1-cv_attempt) * m_percent_id_asd, int_suicide_age = (1-cv_attempt) * m_age)]

dataset[, study := figure_citation]

dataset[, row_id := seq_len(.N)] # redo row_id for ease

# Preliminary descriptive statistics of studies ---------------------------

unique(dataset[,.(country, figure_citation)])
unique(dataset[,.(country)])

dataset[mean == min(mean), .(figure_citation, country, sex, id, age_start, age_end, mean, lower, upper)]
dataset[mean == max(mean), .(figure_citation, country, sex, id, age_start, age_end, mean, lower, upper)]


# Develop meta-regression model -------------------------------------------

#trimming <- 0.95 # for sensitivity testing
trimming <- 1

potential_cvs <- c("cv_attempt", "m_percent_female_asd", "m_percent_id_asd", "m_age", "int_attempt_female", "int_attempt_id", "int_attempt_age")

mr_dataset <- mr$MRData()
mr_dataset$load_df(
  data = dataset,
  col_obs = "log_rr", col_obs_se = "log_rr_se",
  col_covs = as.list(c(potential_cvs, "row_id")), col_study_id = "study" )

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
betas1 # worst is m_age
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas1[, `:=` (aic = aic_current, aic_diff = NA, bic = bic_current, bic_diff = NA)]
aic_previous <- aic_current
bic_previous <- bic_current

## Step 2
selected_covs <- selected_covs[selected_covs != "m_age"]
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas2 <- betas[,.(model = 'Model 2', cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]
betas2 # worst is int_attempt_female 
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas2[, `:=` (aic = aic_current, aic_diff = aic_current - aic_previous, bic = bic_current, bic_diff = bic_current - bic_previous)]
aic_previous <- aic_current
bic_previous <- bic_current
model2 <- model

## Step 3
selected_covs <- selected_covs[selected_covs != "int_attempt_female"]
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas3 <- betas[,.(model = 'Model 3', cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]
betas # all are significant but worst is cv_attempt p =  0.04056558
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas3[, `:=` (aic = aic_current, aic_diff = aic_current - aic_previous, bic = bic_current, bic_diff = bic_current - bic_previous)]
aic_previous <- aic_current
bic_previous <- bic_current
model3 <- model

# Step 4
selected_covs <- selected_covs[selected_covs != "cv_attempt"]
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas4 <- betas[,.(model = "Model 4", cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = p)]
betas4 # all are significant p < 0.05
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas4[, `:=` (aic = aic_current, aic_diff = aic_current - aic_previous, bic = bic_current, bic_diff = bic_current - bic_previous)]
aic_previous <- aic_current
bic_previous <- bic_current
model4 <- model

betas_table <- rbind(betas1, betas2, betas3, betas4)

betas_table <- betas_table[,.(Model = model, Covariate = cov, Coefficient = beta, `Lower UI` = lower, `Upper UI` = upper, p = round(p, 4), AIC = aic, BIC = bic)]
betas_table[Covariate == 'intercept', Covariate := 'Intercept']
betas_table[Covariate == 'cv_attempt', Covariate := 'Suicide attempt (vs suicide)']
betas_table[Covariate == 'm_percent_female_asd', Covariate := 'Percent female']
betas_table[Covariate == 'm_percent_id_asd', Covariate := 'Percent with ID']
betas_table[Covariate == 'm_age', Covariate := 'Mid age']
betas_table[Covariate == 'int_attempt_female', Covariate := 'Interaction: Suicide attempt (vs suicide) X Percent female']
betas_table[Covariate == 'int_attempt_id', Covariate := 'Interaction: Suicide attempt (vs suicide) X Percent with ID']
betas_table[Covariate == 'int_attempt_age', Covariate := 'Interaction: Suicide attempt (vs suicide) X Mid age']

write.csv(betas_table, "/FILEPATH/model_building_table.csv", row.names = F)

# Summary of results ------------------------------------------------------
model <- model3
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]

## Summary RRs
summary_rrs_deaths <- rbind(data.table(intercept = 1, cv_attempt = 0, m_percent_female_asd = 0, m_percent_id_asd = 0,  int_attempt_id = 0, int_attempt_age = 0, m_age = 0),
                        data.table(intercept = 1, cv_attempt = 0, m_percent_female_asd = -0.5, m_percent_id_asd = 0,  int_attempt_id = 0, int_attempt_age = 0, m_age = 0),
                        data.table(intercept = 1, cv_attempt = 0, m_percent_female_asd = 0.5, m_percent_id_asd = 0,  int_attempt_id = 0, int_attempt_age = 0, m_age = 0),
                        data.table(intercept = 1, cv_attempt = 0, m_percent_female_asd = 0, m_percent_id_asd = 1-gbd2020_prop_id,  int_attempt_id = 0, int_attempt_age = 0, m_age = 0),
                        data.table(intercept = 1, cv_attempt = 0, m_percent_female_asd = 0,  m_percent_id_asd = -gbd2020_prop_id,  int_attempt_id = 0, int_attempt_age = 0, m_age = 0))
summary_rrs_attempts <- rbind(data.table(intercept = 1, cv_attempt = 1, m_percent_female_asd = 0, m_percent_id_asd = 0,  int_attempt_id = 0, int_attempt_age = 0, m_age = 0),
                            data.table(intercept = 1, cv_attempt = 1, m_percent_female_asd = -0.5, m_percent_id_asd = 0,  int_attempt_id = 0, int_attempt_age = 0, m_age = 0),
                            data.table(intercept = 1, cv_attempt = 1, m_percent_female_asd = 0.5, m_percent_id_asd = 0,  int_attempt_id = 0, int_attempt_age = 0, m_age = 0),
                            data.table(intercept = 1, cv_attempt = 1, m_percent_female_asd = 0, m_percent_id_asd = 1-gbd2020_prop_id,  int_attempt_id = 1-gbd2020_prop_id, int_attempt_age = 0, m_age = 0),
                            data.table(intercept = 1, cv_attempt = 1, m_percent_female_asd = 0,  m_percent_id_asd = -gbd2020_prop_id,  int_attempt_id = -gbd2020_prop_id, int_attempt_age = 0, m_age = 0))

summary_rrs <- rbind(summary_rrs_deaths, summary_rrs_attempts)

predict_data <- mr$MRData()
predict_data$load_df(data = summary_rrs, col_covs=as.list(model$cov_names))
beta_samples <- mr$core$other_sampling$sample_simple_lme_beta(500L, model)
gamma_outer_samples <- matrix(rep(model$gamma_soln, each = 500L), nrow = 500L)

draws <- model$create_draws(predict_data,
                            beta_samples = beta_samples,
                            gamma_samples = gamma_outer_samples,
                            random_study = F)

summary_rrs$pred <- model$predict(data = predict_data)
summary_rrs$pred_lo <- apply(draws, 1, function(x) quantile(x, 0.025))
summary_rrs$pred_hi <- apply(draws, 1, function(x) quantile(x, 0.975))


summary_rrs[, sex := ifelse(m_percent_female_asd == 0, "Both", ifelse(m_percent_female_asd < 0, "Male", "female"))]
summary_rrs[, ID := ifelse(m_percent_id_asd == 0, "Both", ifelse(m_percent_id_asd < 0, "ID absent", "ID present"))]
summary_rrs[, attempt := ifelse(cv_attempt == 0, "Suicide", "Attempt")]

intercept_for_plot_mortality <- summary_rrs[sex == "Both" & ID == "Both" & attempt == "Suicide",.(pred, pred_lo, pred_hi)]
intercept_for_plot_attempts <- summary_rrs[sex == "Both" & ID == "Both" & attempt == "Attempt",.(pred, pred_lo, pred_hi)]


summary_rrs <- summary_rrs[,.(sex, ID, attempt, rr = round(exp(pred), 2), lower = round(exp(pred_lo),2), upper = round(exp(pred_hi), 2))]

summary_rrs

# Create model plots ------------------------------------------------------
dataset[, figure_citation :=gsub("et al. ", "et al., ", figure_citation)]
dataset[, figure_citation :=gsub("\\(", "", figure_citation)]
dataset[, figure_citation :=gsub(")", "", figure_citation)]

dataset[, figure_citation := trimws(figure_citation)]

setorder(dataset, cols = -cv_attempt, -figure_citation, -percent_female_asd, -percent_id_asd, -age)

dataset[, ord := seq_len(.N)]

## Forest plot
pdf(file='/FILEPATH/forestplot.pdf', width = 11, height = 11)
forest(x = dataset[, log_rr], sei = dataset[, log_rr_se], atransf = exp, slab = dataset[, figure_citation], xlab = "",
       ilab=dataset[, .(round(percent_female_asd*100), round(percent_id_asd*100), round_c(age, 1))], ilab.xpos=c(-5.5,-4.5,-3.5),
       rows = c(6:31, 38:69), ylim=c(3, 74), xlim=c(-9, 6), cex=0.8, psize=1, header = c("Study", "Relative risk [95% UI]"))
op <- par(cex=0.75, font=2)
op <- par(cex=0.8, font=2)
text(c(-5.5,-4.5,-3.5), 73, c("% female", "% ID", "Age (years)"))
par(font=4)
text(-9, c(33.5,70.5), pos=4, c("Suicide attempt", "Suicide mortality"))
addpoly(x = intercept_for_plot_mortality$pred, ci.lb = intercept_for_plot_mortality$pred_lo, ci.ub = intercept_for_plot_mortality$pred_hi, row=36, cex=1.1, atransf=exp, mlab="Pooled Relative Risk")
addpoly(x = intercept_for_plot_attempts$pred, ci.lb = intercept_for_plot_attempts$pred_lo, ci.ub = intercept_for_plot_attempts$pred_hi, row=3, cex=1.1, atransf=exp, mlab="Pooled Relative Risk")
dev.off() # Turn the PDF device off

## Funnel plot
pdf(file='/FILEPATH/funnelplot_bydata.pdf', width = 10, height = 9)
funnel_plot_mrbrt_bydata(model)
dev.off() # Turn the PDF device off

pdf(file='/FILEPATH/funnelplot_bystudy.pdf', width = 10, height = 9)
funnel_plot_mrbrt_bystudy(model)
dev.off() # Turn the PDF device off

egger_mr_brt_pval(model)

# Create prediction matrix ------------------------------------------------

predict_matrix <- rbind(data.table(intercept = 1, cv_attempt = 0, m_percent_female_asd = -0.5, m_percent_id_asd = 0,  int_attempt_id = 0, int_attempt_age = 0),
                        data.table(intercept = 1, cv_attempt = 0, m_percent_female_asd = 0.5, m_percent_id_asd = 0,  int_attempt_id = 0, int_attempt_age = 0),
                        data.table(intercept = 1, cv_attempt = 0, m_percent_female_asd = -0.5, m_percent_id_asd = 1-gbd2020_prop_id, int_attempt_id = 0, int_attempt_age = 0),
                        data.table(intercept = 1, cv_attempt = 0, m_percent_female_asd = 0.5, m_percent_id_asd = 1-gbd2020_prop_id,  int_attempt_id = 0, int_attempt_age = 0),
                        data.table(intercept = 1, cv_attempt = 0, m_percent_female_asd = -0.5, m_percent_id_asd = -gbd2020_prop_id, int_attempt_id = 0, int_attempt_age = 0),
                        data.table(intercept = 1, cv_attempt = 0, m_percent_female_asd = 0.5, m_percent_id_asd = -gbd2020_prop_id,  int_attempt_id = 0, int_attempt_age = 0))

predict_data <- mr$MRData()
predict_data$load_df(data = predict_matrix, col_covs=as.list(model$cov_names))
beta_samples <- mr$core$other_sampling$sample_simple_lme_beta(500L, model)
gamma_outer_samples <- matrix(rep(model$gamma_soln, each = 500L), nrow = 500L)
draws <- model$create_draws(predict_data,
                                beta_samples = beta_samples,
                                gamma_samples = gamma_outer_samples,
                                random_study = F)

predict_matrix$pred <- model$predict(data = predict_data)
predict_matrix$pred_lo <- apply(draws, 1, function(x) quantile(x, 0.025))
predict_matrix$pred_hi <- apply(draws, 1, function(x) quantile(x, 0.975))

# Bind and save draws -----------------------------------------------------
predict_matrix <- cbind(predict_matrix, draws)

predict_matrix <- melt.data.table(predict_matrix[m_percent_id_asd != 0,], id.vars = names(predict_matrix)[!(names(predict_matrix) %like% "V")], value.name="rr", variable.name="draw")
predict_matrix[, draw := gsub("V", "draw_", draw)]
predict_matrix[draw == "draw_500", draw := "draw_0"]
predict_matrix[, `:=` (sex_id = m_percent_female_asd + 1.5, intel_dis = ifelse(m_percent_id_asd < 0, 0, 1), rr = exp(rr))]
predict_matrix <- unique(predict_matrix[,.(sex_id, intel_dis, rr, draw)])

write.csv(predict_matrix, "/FILEPATH/rr_matrix.csv", row.names = F)

# Sensitivity analyses #1 - 5% trimming -----------------------------------
trimming <- 0.95

potential_cvs <- c("cv_attempt", "m_percent_female_asd", "m_percent_id_asd", "m_age", "int_attempt_female", "int_attempt_id", "int_attempt_age")

mr_dataset <- mr$MRData()
mr_dataset$load_df(
  data = dataset[age != 49.5],
  col_obs = "log_rr", col_obs_se = "log_rr_se",
  col_covs = as.list(c(potential_cvs, "row_id")), col_study_id = "study" )

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
betas1
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas1[, `:=` (aic = aic_current, aic_diff = NA, bic = bic_current, bic_diff = NA)]
aic_previous <- aic_current
bic_previous <- bic_current

## Step 2
selected_covs <- selected_covs[selected_covs != "int_attempt_id"]
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas2 <- betas[,.(model = 'Model 2', cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]
betas2
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas2[, `:=` (aic = aic_current, aic_diff = aic_current - aic_previous, bic = bic_current, bic_diff = bic_current - bic_previous)]
aic_previous <- aic_current
bic_previous <- bic_current

## Step 3
selected_covs <- selected_covs[selected_covs != "cv_attempt"]
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas3 <- betas[,.(model = 'Model 3', cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]
betas # all are significant but worst is cv_attempt 
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas3[, `:=` (aic = aic_current, aic_diff = aic_current - aic_previous, bic = bic_current, bic_diff = bic_current - bic_previous)]
aic_previous <- aic_current
bic_previous <- bic_current
model3 <- model

betas_table <- rbind(betas1, betas2, betas3)#, betas4)

betas_table <- betas_table[,.(Model = model, Covariate = cov, Coefficient = beta, `Lower UI` = lower, `Upper UI` = upper, p = round(p, 4), AIC = aic, BIC = bic)]
betas_table[Covariate == 'intercept', Covariate := 'Intercept']
betas_table[Covariate == 'cv_attempt', Covariate := 'Suicide attempt (vs suicide)']
betas_table[Covariate == 'm_percent_female_asd', Covariate := 'Percent female']
betas_table[Covariate == 'm_percent_id_asd', Covariate := 'Percent with ID']
betas_table[Covariate == 'm_age', Covariate := 'Mid age']
betas_table[Covariate == 'int_attempt_female', Covariate := 'Interaction: Suicide attempt (vs suicide) X Percent female']
betas_table[Covariate == 'int_attempt_id', Covariate := 'Interaction: Suicide attempt (vs suicide) X Percent with ID']
betas_table[Covariate == 'int_attempt_age', Covariate := 'Interaction: Suicide attempt (vs suicide) X Mid age']

# Sensitivity analyses #2 - 100% ID assumption -----------------------------------

dataset_100_id <- copy(dataset)
dataset_100_id[percent_id_asd == gbd2020_prop_id, percent_id_asd := 1]
dataset_100_id[, `:=` (m_percent_id_asd = percent_id_asd - gbd2020_prop_id)]
dataset_100_id[, `:=` (int_attempt_id = cv_attempt * m_percent_id_asd)]

trimming <- 1

potential_cvs <- c("cv_attempt", "m_percent_female_asd", "m_percent_id_asd", "m_age", "int_attempt_female", "int_attempt_id", "int_attempt_age")

mr_dataset <- mr$MRData()
mr_dataset$load_df(
  data = dataset_100_id,
  col_obs = "log_rr", col_obs_se = "log_rr_se",
  col_covs = as.list(c(potential_cvs, "row_id")), col_study_id = "study" )

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
betas1 # worst is  m_age 
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas1[, `:=` (aic = aic_current, aic_diff = NA, bic = bic_current, bic_diff = NA)]
aic_previous <- aic_current
bic_previous <- bic_current

## Step 2
selected_covs <- selected_covs[selected_covs != "m_age"]
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas2 <- betas[,.(model = 'Model 2', cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]
betas2 # worst is int_attempt_female
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas2[, `:=` (aic = aic_current, aic_diff = aic_current - aic_previous, bic = bic_current, bic_diff = bic_current - bic_previous)]
aic_previous <- aic_current
bic_previous <- bic_current

## Step 3
selected_covs <- selected_covs[selected_covs != "int_attempt_female"]
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas3 <- betas[,.(model = 'Model 3', cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]
betas # all are significant but worst is cv_attempt 
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas3[, `:=` (aic = aic_current, aic_diff = aic_current - aic_previous, bic = bic_current, bic_diff = bic_current - bic_previous)]
aic_previous <- aic_current
bic_previous <- bic_current
model3 <- model

# Step 4
selected_covs <- selected_covs[selected_covs != "cv_attempt"]
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas4 <- betas[,.(model = "Model 4", cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = p)]
betas4 # all are significant p < 0.001
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas4[, `:=` (aic = aic_current, aic_diff = aic_current - aic_previous, bic = bic_current, bic_diff = bic_current - bic_previous)]
aic_previous <- aic_current
bic_previous <- bic_current
model4 <- model

betas_table <- rbind(betas1, betas2, betas3, betas4)

betas_table <- betas_table[,.(Model = model, Covariate = cov, Coefficient = beta, `Lower UI` = lower, `Upper UI` = upper, p = round(p, 4), AIC = aic, BIC = bic)]
betas_table[Covariate == 'intercept', Covariate := 'Intercept']
betas_table[Covariate == 'cv_attempt', Covariate := 'Suicide attempt (vs suicide)']
betas_table[Covariate == 'm_percent_female_asd', Covariate := 'Percent female']
betas_table[Covariate == 'm_percent_id_asd', Covariate := 'Percent with ID']
betas_table[Covariate == 'm_age', Covariate := 'Mid age']
betas_table[Covariate == 'int_attempt_female', Covariate := 'Interaction: Suicide attempt (vs suicide) X Percent female']
betas_table[Covariate == 'int_attempt_id', Covariate := 'Interaction: Suicide attempt (vs suicide) X Percent with ID']
betas_table[Covariate == 'int_attempt_age', Covariate := 'Interaction: Suicide attempt (vs suicide) X Mid age']

# Sensitivity analyses #3 - 0% ID assumption -----------------------------------

dataset_0_id <- copy(dataset)
dataset_0_id[percent_id_asd == gbd2020_prop_id, percent_id_asd := 0]
dataset_0_id[percent_id_asd == gbd2020_prop_id, percent_id_asd := 1]
dataset_0_id[, `:=` (m_percent_id_asd = percent_id_asd - gbd2020_prop_id)]
dataset_0_id[, `:=` (int_attempt_id = cv_attempt * m_percent_id_asd)]

trimming <- 1

potential_cvs <- c("cv_attempt", "m_percent_female_asd", "m_percent_id_asd", "m_age", "int_attempt_female", "int_attempt_id", "int_attempt_age")

mr_dataset <- mr$MRData()
mr_dataset$load_df(
  data = dataset_0_id,
  col_obs = "log_rr", col_obs_se = "log_rr_se",
  col_covs = as.list(c(potential_cvs, "row_id")), col_study_id = "study" )

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
betas1 # worst is  m_age  
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas1[, `:=` (aic = aic_current, aic_diff = NA, bic = bic_current, bic_diff = NA)]
aic_previous <- aic_current
bic_previous <- bic_current

## Step 2
selected_covs <- selected_covs[selected_covs != "m_age"]
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas2 <- betas[,.(model = 'Model 2', cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]
betas2 # worst is int_attempt_female
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas2[, `:=` (aic = aic_current, aic_diff = aic_current - aic_previous, bic = bic_current, bic_diff = bic_current - bic_previous)]
aic_previous <- aic_current
bic_previous <- bic_current

## Step 3
selected_covs <- selected_covs[selected_covs != "int_attempt_female"]
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas3 <- betas[,.(model = 'Model 3', cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]
betas # all are significant but worst is cv_attempt 
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas3[, `:=` (aic = aic_current, aic_diff = aic_current - aic_previous, bic = bic_current, bic_diff = bic_current - bic_previous)]
aic_previous <- aic_current
bic_previous <- bic_current
model3 <- model

# Step 4
selected_covs <- selected_covs[selected_covs != "cv_attempt"]
cov_list <- list(mr$LinearCovModel('intercept', use_re = T))
for(c in selected_covs){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)
model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas[p == betas[cov != "intercept", max(p) ] & cov != "intercept", .(cov, p)]
betas4 <- betas[,.(model = "Model 4", cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = p)]
betas4 # all are significant p < 0.001
aic_current <- mrbrt_aic(model)
bic_current <- mrbrt_bic(model)
betas4[, `:=` (aic = aic_current, aic_diff = aic_current - aic_previous, bic = bic_current, bic_diff = bic_current - bic_previous)]
aic_previous <- aic_current
bic_previous <- bic_current
model4 <- model

betas_table <- rbind(betas1, betas2, betas3, betas4)

betas_table <- betas_table[,.(Model = model, Covariate = cov, Coefficient = beta, `Lower UI` = lower, `Upper UI` = upper, p = round(p, 4), AIC = aic, BIC = bic)]
betas_table[Covariate == 'intercept', Covariate := 'Intercept']
betas_table[Covariate == 'cv_attempt', Covariate := 'Suicide attempt (vs suicide)']
betas_table[Covariate == 'm_percent_female_asd', Covariate := 'Percent female']
betas_table[Covariate == 'm_percent_id_asd', Covariate := 'Percent with ID']
betas_table[Covariate == 'm_age', Covariate := 'Mid age']
betas_table[Covariate == 'int_attempt_female', Covariate := 'Interaction: Suicide attempt (vs suicide) X Percent female']
betas_table[Covariate == 'int_attempt_id', Covariate := 'Interaction: Suicide attempt (vs suicide) X Percent with ID']
betas_table[Covariate == 'int_attempt_age', Covariate := 'Interaction: Suicide attempt (vs suicide) X Mid age']

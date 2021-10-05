rm(list=ls())
library(openxlsx)
library(data.table)
library(ggplot2)
library(crosswalk, lib.loc = "/FILEPATH/")
library(mrbrt001, lib.loc = '/FILEPATH/')
source("/FILEPATH/get_population.R")
source("/FILEPATH/get_location_metadata.R")
source("/FILEPATH/get_ids.R")
rlogit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}

dataset <- data.table(read.xlsx("/FILEPATH/covid19_extraction_sheet.xlsx"))[!1]
  
# data prep ---------------------------------------------------------------
dataset[is.na(cohort), cohort := gsub(",", "", unlist(strsplit(field_citation_value, " "))[1]),  by = "field_citation_value"]

numeric_columns <- c("mean_covid", "mean_ref", "se_covid", "se_ref", "sample_size_covid", "sample_size_ref", "lower_covid", "lower_ref", "upper_covid", "upper_ref", "age_start", "age_end", "date_ref_start", "date_ref_end", "date_covid_start", "date_covid_end", "social_mobility")
dataset[, (numeric_columns):=lapply(.SD, as.numeric), .SD=numeric_columns]

dataset[, `:=` (covid_cases = mean_covid * sample_size_covid, ref_cases = mean_ref * sample_size_ref)]
dataset[!is.na(lower_covid), se_covid := (upper_covid - lower_covid)/(qnorm(0.975,0, 1)*2)]
dataset[!is.na(lower_ref), se_ref := (upper_ref - lower_ref)/(qnorm(0.975,0, 1)*2)]
dataset[is.na(se_covid), se_covid := sqrt((mean_covid*(1-mean_covid))/sample_size_covid)]
dataset[is.na(se_ref), se_ref := sqrt((mean_ref*(1-mean_ref))/sample_size_ref)]
dataset[, `:=` (cv_dep = ifelse(case_diagnostics %in% c("SMFQ", "PHQ-8", "PHQ-9", "PHQ-2", "MFQ", "CES-D", "DASS-D", "WHO-5", "HADS - Depression"), 1, 0), cv_anx = ifelse(case_diagnostics %in% c("DASS-A", "GAD-7", "HBQ", "HAD - Anxiety", "SCARED - Generalised Anxiety subscale", "STAI", "HADS - Anxiety"), 1, 0))]
table(dataset[,.(case_diagnostics, cv_dep)])
table(dataset[,.(case_diagnostics, cv_anx)])
dataset[case_definition == "Major depressive disorder", cv_dep := 1]
dataset[case_definition == "Major depressive episode", cv_dep := 1]
dataset[case_definition == "Anxiety disorder", cv_anx := 1]
dataset[case_definition == "Anxiety disorders", cv_anx := 1]
table(dataset[,.(case_diagnostics, cv_dep)])
table(dataset[,.(case_diagnostics, cv_anx)])

dep_studies <- dataset[cv_dep == 1, field_citation_value]
anx_studies <- dataset[cv_anx == 1, field_citation_value]
both_studies <- dep_studies[dep_studies %in% anx_studies]

# Apply logit transformations
dataset[, row_id := seq_len(.N)]
dataset[, `:=` (logit_covid = delta_transform(mean = mean_covid, sd = se_covid, transformation = "linear_to_logit")[1],
                logit_covid_se = delta_transform(mean = mean_covid, sd = se_covid, transformation = "linear_to_logit")[2],
                logit_ref = delta_transform(mean = mean_ref, sd = se_ref, transformation = "linear_to_logit")[1],
                logit_ref_se = delta_transform(mean = mean_ref, sd = se_ref, transformation = "linear_to_logit")[2]), by = "row_id"]

## Subset data to most-unique estimates
dataset[, unique_label := paste0(cohort, "_", case_definition, "_", location_id, "_", date_covid_start, "_", date_covid_end)]
dataset[, `:=` (widest_age = ifelse(age_start == min(age_start) & age_end == max(age_end), 1, 0)), by = "unique_label"]

## make list of studies with both-sex estimates
both_sex_studies <- unique(dataset[sex == "Both", unique_label])
sex_specific_studies <- unique(dataset[sex != "Both", unique_label])
only_both_sex_studies <- both_sex_studies[!(both_sex_studies %in% sex_specific_studies)]

all_age_studies <- unique(dataset[widest_age == 1, unique_label])
age_specific_studies <- unique(dataset[widest_age == 0, unique_label])
only_all_age_studies <- all_age_studies[!(all_age_studies %in% age_specific_studies)]

## Both sex, all age
dataset_bothsex_allage <- dataset[(unique_label %in% only_both_sex_studies) & (unique_label %in% only_all_age_studies),]

## Sex specific, all age
dataset_sexspecific_allage <- dataset[(unique_label %in% sex_specific_studies) & (unique_label %in% only_all_age_studies) & sex != "Both",]

## both sex, age specific
dataset_bothsex_agespecific <- dataset[(unique_label %in% only_both_sex_studies) & (unique_label %in% age_specific_studies) & widest_age == 0,]

## age and sex specific estimates
dataset_sexspecific_agespecific <- dataset[sex != "Both" & widest_age == 0,]

## Bind all data that doesn't need splitting
dataset_raw <- rbind(dataset_bothsex_allage, dataset_sexspecific_allage, dataset_bothsex_agespecific, dataset_sexspecific_agespecific)

# Age-sex split where applicable
dataset_agesexsplitting <- dataset[!(unique_label %in% dataset_raw$unique_label),]

sex_ratios_both <- dataset_agesexsplitting[sex == "Both" & widest_age == 1, .(unique_label, location_id, date_covid_start, date_covid_end, logit_covid_both = logit_covid, logit_covid_both_se = logit_covid_se, logit_ref_both = logit_ref, logit_ref_both_se = logit_ref_se)]
sex_ratios_male <- dataset_agesexsplitting[sex == "Male" & widest_age == 1, .(unique_label, location_id, date_covid_start, date_covid_end, logit_covid_male = logit_covid, logit_covid_male_se = logit_covid_se, logit_ref_male = logit_ref, logit_ref_male_se = logit_ref_se)]
sex_ratios_female <- dataset_agesexsplitting[sex == "Female" & widest_age == 1, .(unique_label, location_id, date_covid_start, date_covid_end, logit_covid_female = logit_covid, logit_covid_female_se = logit_covid_se, logit_ref_female = logit_ref, logit_ref_female_se = logit_ref_se)]
sex_ratios <- merge(sex_ratios_both, sex_ratios_male, by = c("unique_label", "location_id", "date_covid_start", "date_covid_end"))
sex_ratios <- merge(sex_ratios, sex_ratios_female, by = c("unique_label", "location_id", "date_covid_start", "date_covid_end"))
sex_ratios[, `:=` (covid_male_crosswalk = logit_covid_both - logit_covid_male, covid_male_crosswalk_se = sqrt(logit_covid_both_se^2 + logit_covid_male_se^2),
                   covid_female_crosswalk = logit_covid_both - logit_covid_female, covid_female_crosswalk_se = sqrt(logit_covid_both_se^2 + logit_covid_female_se^2))]
sex_ratios[, `:=` (ref_male_crosswalk = logit_ref_both - logit_ref_male, ref_male_crosswalk_se = sqrt(logit_ref_both_se^2 + logit_ref_male_se^2),
                   ref_female_crosswalk = logit_ref_both - logit_ref_female, ref_female_crosswalk_se = sqrt(logit_ref_both_se^2 + logit_ref_female_se^2))]
sex_ratios <- sex_ratios[,.(unique_label, location_id, date_covid_start, date_covid_end, covid_male_crosswalk, covid_male_crosswalk_se, covid_female_crosswalk, covid_female_crosswalk_se, ref_male_crosswalk, ref_male_crosswalk_se, ref_female_crosswalk, ref_female_crosswalk_se)]

dataset_age_male <- merge(dataset_agesexsplitting[widest_age == 0,], sex_ratios[,.(unique_label, location_id, date_covid_start, date_covid_end, covid_male_crosswalk, covid_male_crosswalk_se, ref_male_crosswalk, ref_male_crosswalk_se)], by = c("unique_label", "location_id", "date_covid_start", "date_covid_end"))
dataset_age_male[,`:=`(sex = "Male", percent_female = 0, logit_covid = logit_covid - covid_male_crosswalk, logit_covid_se = sqrt(logit_covid_se^2 + covid_male_crosswalk_se^2),
                       logit_ref= logit_ref - ref_male_crosswalk, logit_ref_se = sqrt(logit_ref_se^2 + ref_male_crosswalk_se^2))]
dataset_age_female <- merge(dataset_agesexsplitting[widest_age == 0,], sex_ratios[,.(unique_label, location_id, date_covid_start, date_covid_end, covid_female_crosswalk, covid_female_crosswalk_se, ref_female_crosswalk, ref_female_crosswalk_se)], by = c("unique_label", "location_id", "date_covid_start", "date_covid_end"))
dataset_age_female[,`:=`(sex = "Female", percent_female = 1, logit_covid = logit_covid - covid_female_crosswalk, logit_covid_se = sqrt(logit_covid_se^2 + covid_female_crosswalk_se^2),
                         logit_ref= logit_ref - ref_female_crosswalk, logit_ref_se = sqrt(logit_ref_se^2 + ref_female_crosswalk_se^2))]

dataset_final <- rbind(dataset_raw, dataset_age_male, dataset_age_female, fill = T)
  
dataset_final[, `:=` (mid_age = (age_start + age_end) / 2, follow_up_days = ((date_covid_start+date_covid_end)/2) -  ((date_ref_start+date_ref_end)/2))]
dataset_final[, `:=` (logit_dif = logit_ref-logit_covid, logit_dif_se = sqrt(logit_ref_se^2+logit_covid_se^2))]
dataset_final[is.na(percent_female) & sex == "Both", percent_female := 0.5]

location_meta <- get_location_metadata(location_set_id = 35, gbd_round_id = 7, decomp_step = 'iterative')
dataset_final <- merge(dataset_final, location_meta[,.(location_id = as.character(location_id), parent_id)], by = "location_id", all.x = T)

# Re-incorporate recall of questionnaire ----------------------------------
dataset_final[case_diagnostics %in% c("BOTH (dass-d and dass-a)", "CES-D", "DASS-A", "DASS-D", "HAD - Anxiety",  "HADS - Anxiety",  "HADS - Depression"), recall_days := 7]
dataset_final[case_diagnostics %in% c("GAD-7", "MFQ", "PHQ-9", "PHQ-8", "SMFQ", "WHO-5"), recall_days := 14]
dataset_final[case_diagnostics %in% c("GHQ-12"), recall_days := 21]
dataset_final[case_diagnostics %in% c("CIDI", "K6", "K10", "MHI-5"), recall_days := 30]

dataset_final[case_diagnostics == "MINI" & case_definition == "Major depressive disorder", recall_days := 14]

# Ambiguous recalls capped at past month:
dataset_final[case_diagnostics == "MINI" & case_definition == "Anxiety disorder", recall_days := 90] 
dataset_final[case_diagnostics == "HBQ", recall_days := 30] 
dataset_final[case_diagnostics == "SCARED - Generalised Anxiety subscale", recall_days := 90] 

dataset_final[, date_covid_start := date_covid_start - recall_days]

#  IHME infectionator estimates ----------------------------------------

## Deaths
ihme_deaths <- fread("/FILEPATH/daily_deaths.csv")

dataset_final[location_id %in% ihme_deaths$location_id, location_id_deaths := location_id]
dataset_final[!(location_id %in% ihme_deaths$location_id), location_id_deaths := parent_id]

ihme_deaths <- ihme_deaths[location_id %in% dataset_final$location_id_deaths]
ihme_deaths <- melt.data.table(ihme_deaths, id.vars = names(ihme_deaths)[!(names(ihme_deaths) %like% "draw")], value.name="deaths", variable.name="draw")
ihme_deaths[is.na(deaths), deaths := 0]
ihme_deaths[, deaths := mean(deaths), by = c("location_id", "date")]
ihme_deaths <- unique(ihme_deaths[,.(location_id, date, deaths)])

ihme_deaths[, date := as.numeric(as.Date(paste0(date)) - as.Date(0, origin="1899-12-30", tz='UTC')), by = "date"]

population <- get_population(age_group_id = 22, location_id = unique(dataset_final$location_id_deaths), year_id = 2020, gbd_round_id = 7, decomp_step = 'iterative')
ihme_deaths <- merge(ihme_deaths, population, by = "location_id", all.x = T)
ihme_deaths[, deaths := deaths / population]

## Infections
ihme_infections <- fread("/FILEPATH/daily_infections.csv")

ihme_infections <- ihme_infections[location_id %in% dataset_final$location_id_deaths]
ihme_infections <- melt.data.table(ihme_infections, id.vars = names(ihme_infections)[!(names(ihme_infections) %like% "draw")], value.name="infections", variable.name="draw")
ihme_infections[is.na(infections), infections := 0]
ihme_infections[, infections := mean(infections), by = c("location_id", "date")]
ihme_infections <- unique(ihme_infections[,.(location_id, date, infections)])

ihme_infections[, date := as.numeric(as.Date(paste0(date)) - as.Date(0, origin="1899-12-30", tz='UTC')), by = "date"]

ihme_infections <- merge(ihme_infections, population, by = "location_id", all.x = T)
ihme_infections[, infections := infections / population]

# IHME mobility estimates -------------------------------------------------
mobility_data <- fread("/FILEPATH/mobility.csv") 

mobility_data[, date := as.numeric(as.Date(paste0(date)) - as.Date(0, origin="1899-12-30", tz='UTC')), by = "date"]
mobility_data[, mobility_reference := mean]

dataset_final[location_id %in% mobility_data$location_id, location_id_mob := location_id]
dataset_final[!(location_id %in% mobility_data$location_id), location_id_mob := parent_id]

# Match date-specific covariates to dataset -------------------------------
dataset_final[location_id != location_id_mob & location_id_mob == location_id_deaths, location_id := parent_id]

population <- get_population(age_group_id = 22, location_id = unique(dataset_final$location_id), gbd_round_id = 7, year_id = 2020, decomp_step = 'iterative')
  
dataset_final[, row_id := seq_len(.N)] # re-do row IDs
for(r in dataset_final$row_id){
  date_range <- dataset_final[row_id == r, seq(date_covid_start,date_covid_end)]
  date_range_lag_1m <- dataset_final[row_id == r, seq(date_covid_start,date_covid_end)]-30
  date_range_lag_2m <- dataset_final[row_id == r, seq(date_covid_start,date_covid_end)]-60
  location <- as.numeric(dataset_final[row_id == r, location_id])
  location_infectionator <- as.numeric(dataset_final[row_id == r, location_id_deaths])
  # death data
  deaths <- ihme_deaths[date %in% date_range & location_id == location_infectionator, deaths]
  if(length(deaths) == 0){print(paste0("No deaths data for ", dataset_final[row_id == r, location_name]))}
  daily_rate_deaths <- mean(deaths)
  ## Infections
  infections <- ihme_infections[date %in% date_range & location_id == location_infectionator, infections]
  if(length(infections) == 0){print(paste0("No infections data for ", dataset_final[row_id == r, location_name]))}
  daily_rate_infections <- mean(infections)
  # mobility data
  mobility <- mean(mobility_data[date %in% date_range & location_id == location,mobility_reference])/100

  if(length(mobility) == 0){print(paste0("No mobility data for ", dataset_final[row_id == r, location_name]))}
  dataset_final[row_id == r, `:=` (covid_death_daily = daily_rate_deaths, covid_infections_daily = daily_rate_infections, social_mobility = mobility)]

  if(r %in% seq(1, dataset_final[,max(row_id)], round(dataset_final[,max(row_id)]/10))){print(paste0(round(r/dataset_final[,max(row_id)], 2)*100, "% complete"))}
}
 
dataset_final[, covid_death_log := log(covid_death_daily)]
dataset_final[, covid_death_sqrt := sqrt(covid_death_daily)]

dataset_final[, covid_infections_log := log(covid_infections_daily)]
dataset_final[, covid_infections_sqrt := sqrt(covid_infections_daily)]

# Create interaction terms ------------------------------------------------

mean_mid_age <- dataset_final[, mean(mid_age)]
mean_mid_age
dataset_final[, sd(mid_age)]
dataset_final[, m_mid_age := mid_age - mean_mid_age]
dataset_final[, m_percent_female := percent_female -0.5]

## Fix online panel * cross=sectional interaction
dataset_final[, cv_cross_panel := 0]
dataset_final[cv_cross_sectional == 1 & cv_online_panel == 1, cv_cross_panel := 1]
dataset_final[cv_cross_panel == 1, `:=` (cv_cross_sectional = 0, cv_online_panel = 0)]

## test making k-6 etc. measures 1 on both dep and anx CVs
dataset_final[cv_anx == 0 & cv_dep == 0, `:=` (cv_anx = 1, cv_dep = 1)]

## Social mobility covariates
dataset_final[, `:=` (int_soc_mob_sex = social_mobility*m_percent_female,
                      int_soc_mob_age = social_mobility*m_mid_age,
                      int_soc_cross = social_mobility*cv_cross_sectional,
                      int_soc_panel = social_mobility*cv_online_panel)]

## Mortality covariates
dataset_final[, `:=` (int_death_sqrt_sex = covid_death_sqrt*m_percent_female,
                      int_death_sqrt_age = covid_death_sqrt*m_mid_age,
                      int_death_sqrt_cross = covid_death_sqrt*cv_cross_sectional,
                      int_death_sqrt_panel = covid_death_sqrt*cv_online_panel)]

## Infections covariates
dataset_final[, `:=` (int_infections_sqrt_sex = covid_infections_sqrt*m_percent_female,
                      int_infections_sqrt_age = covid_infections_sqrt*m_mid_age,
                      int_infections_sqrt_cross = covid_infections_sqrt*cv_cross_sectional,
                      int_infections_sqrt_panel = covid_infections_sqrt*cv_online_panel)]

## covariates for disorder-specific models
dataset_final[, `:=` (cv_both_md = ifelse(cv_anx == 1 & cv_dep == 1,1, 0))]
dataset_final[, `:=` (int_soc_both_md = social_mobility*cv_both_md,
                      int_death_daily_both_md = covid_death_daily*cv_both_md,
                      int_death_log_both_md = covid_death_log*cv_both_md,
                      int_death_sqrt_both_md = covid_death_sqrt*cv_both_md,
                      int_infections_daily_both_md = covid_infections_daily*cv_both_md,
                      int_infections_log_both_md = covid_infections_log*cv_both_md,
                      int_infections_sqrt_both_md = covid_infections_sqrt*cv_both_md)]

# Run MR-BRT --------------------------------------------------------------

## Histagrams showing skew of primary indicators ensuring each unique observation within each cohort is counted once
hist(unique(dataset_final[,.(social_mobility, cohort)])$social_mobility)
hist(unique(dataset_final[,.(covid_infections_daily, cohort)])$covid_infections_daily, main = "", ylab = "Frequency of unique estimates by cohort", xlab = "Estimated daily COVID-19 infection rate")
hist(unique(dataset_final[,.(covid_infections_sqrt, cohort)])$covid_infections_sqrt, main = "", ylab = "Frequency of unique estimates by cohort", xlab = "Square root of estimated daily COVID-19 infection rate")
hist(unique(dataset_final[,.(covid_infections_log, cohort)])$covid_infections_log, main = "", ylab = "Frequency of unique estimates by cohort", xlab = "Log of estimated daily COVID-19 infection rate")
hist(unique(dataset_final[,.(covid_death_daily, cohort)])$covid_death_daily, main = "", ylab = "Frequency of unique estimates by cohort", xlab = "Estimated daily excess mortality rate")
hist(unique(dataset_final[,.(covid_death_sqrt, cohort)])$covid_death_sqrt, main = "", ylab = "Frequency of unique estimates by cohort", xlab = "Square root of estimated daily excess mortality rate")
hist(unique(dataset_final[,.(covid_death_log, cohort)])$covid_death_log, main = "", ylab = "Frequency of unique estimates by cohort", xlab = "Log of estimated daily excess mortality rate")

## Correlations between indicator x demographic interactions
cor(dataset_final[,.(int_death_sqrt_age, int_soc_mob_age)])
cor(dataset_final[,.(int_death_sqrt_age, int_infections_sqrt_age)])
cor(dataset_final[,.(int_infections_sqrt_age, int_soc_mob_age)])

cor(dataset_final[,.(int_death_sqrt_sex, int_soc_mob_sex)])
cor(dataset_final[,.(int_death_sqrt_sex, int_infections_sqrt_sex)])
cor(dataset_final[,.(int_infections_sqrt_sex, int_soc_mob_sex)])

## Correlations between indicators
cor(dataset_final[,.(covid_death_sqrt, covid_infections_sqrt)])
cor(dataset_final[,.(covid_death_daily, social_mobility)])
cor(dataset_final[,.(covid_infections_sqrt, social_mobility)])

# Functions for draws -----------------------------------------------------
get_beta_vcov <- function(model){
  model_specs <- mrbrt001::core$other_sampling$extract_simple_lme_specs(model)
  beta_hessian <- mrbrt001::core$other_sampling$extract_simple_lme_hessian(model_specs)
  solve(beta_hessian)
}

get_beta_sd <- function(model){
  beta_sd <- sqrt(diag(get_beta_vcov(model)))
  names(beta_sd) <- model$cov_names
  return(beta_sd)
}

get_gamma_sd <- function(model){
  gamma <- model$gamma_soln
  gamma_fisher <- model$lt$get_gamma_fisher(gamma)
  return(sqrt(diag(solve(gamma_fisher))))
}

get_beta_with_alpha <- function(model){
  beta_sd <- sqrt(diag(get_beta_vcov(model)))
  alpha <- sd((model$data$obs - model$predict(model$data, predict_for_study=T))/model$data$obs_se)
  return(beta_sd <- beta_sd * alpha)
}


# Specify disorder-specific datasets --------------------------------------
dataset_final[, social_mobility := -social_mobility] # to make it represent 'decreasing human mobility'

dataset_mdd <- dataset_final[!(cv_anx == 1 & cv_dep == 0),] # For MDD
dataset_anx <- dataset_final[!(cv_anx == 0 & cv_dep == 1),] # For amx

mean_mid_age_mdd <- dataset_mdd[, mean(mid_age)]
mean_mid_age_anx <- dataset_anx[, mean(mid_age)]

dataset_mdd[, m_mid_age := mid_age - mean_mid_age_mdd]
dataset_anx[, m_mid_age := mid_age - mean_mid_age_anx]

# Creation of impact indicator --------------------------------------------
trimming <- 0.95

indicator_variables <- c("social_mobility", "covid_infections_sqrt", "covid_death_sqrt")

indicator_model_table_mdd <- data.table(disorder = "mdd", cov = indicator_variables)
indicator_model_table_anx <- data.table(disorder = "anx", cov = indicator_variables)

for(d in c("mdd", "anx")){
    potential_covs <- indicator_variables

    mr_dataset <- MRData()
    if(d == "mdd"){
      mr_dataset$load_df(
        data = dataset_mdd,
        col_obs = "logit_dif", col_obs_se = "logit_dif_se",
        col_covs = as.list(c(potential_covs, "cv_cross_panel", "row_id")), col_study_id = "cohort" )
    } else if(d == "anx"){
      mr_dataset$load_df(
        data = dataset_anx,
        col_obs = "logit_dif", col_obs_se = "logit_dif_se",
        col_covs = as.list(c(potential_covs, "cv_cross_panel", "row_id")), col_study_id = "cohort" )
    }

    selected_covs <- potential_covs

    cov_list <- list(LinearCovModel('intercept', use_re = T,  prior_beta_uniform = c(0, 0)))
    
    for(c in selected_covs){
      if(c == "covid_death_sqrt"){
      cov_list <- c(cov_list, list(LinearCovModel(c, use_re = T, prior_beta_uniform = c(-Inf, 0))))
      } else {
      cov_list <- c(cov_list, list(LinearCovModel(c, use_re = T)))
      }
    }

    model <- MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)

    model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)

    betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
    betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
    betas <- betas[,.(cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]

    used_data <- data.table(cbind(model$data$to_df(), data.frame(w = model$w_soln)))

    log_lik <- -model$get_objective()
    
    aic_current <- -2 * log_lik + 2 * (length(model$beta_soln) + length(model$gamma_soln))
  
    worst_p <- betas[cov != "intercept", max(p) ]
    betas[p == worst_p & cov != "intercept", cov]
    betas[, disorder := d]

    assign(paste0('indicator_model_table_', d), merge(get(paste0("indicator_model_table_", d)), betas[,.(disorder, cov, beta_initial = beta, lower_initial = lower, upper_initial = upper, p_initial = p)], all.x = T, by = c("disorder", "cov")))
    old_betas <- as.numeric(model$beta_soln)[2:length(model$beta_soln)]

    py_save_object(object = model, filename = paste0("/FILEPATH/", d, "_", "indicator_model.pkl"), pickle = "dill")

    assign(paste0('indicator_model_', d), model)
    assign(paste0('indicator_model_table_', d), merge(get(paste0("indicator_model_table_", d)), betas[,.(disorder, cov, beta_final = beta, lower_final = lower, upper_final = upper, p_final = p)], all.x = T, by = c("disorder", "cov")))
}

indicator_model_table <- rbind(indicator_model_table_mdd, indicator_model_table_anx)
indicator_model_table[, `:=` (beta_initial = round(beta_initial, 1), lower_initial = round(lower_initial, 1), upper_initial = round(upper_initial, 1), p_initial = round(p_initial, 4), beta_final = round(beta_final, 1),  lower_final = round(lower_final, 1), upper_final = round(upper_final, 1), p_final = round(p_final , 4))]

indicator_model_table

write.csv(indicator_model_table, "/FILEPATH/indicator_model_table.csv", row.names = F, na = "-")

dataset_mdd[, indicator := social_mobility * as.numeric(indicator_model_mdd$fe_soln["social_mobility"]) + covid_infections_sqrt * as.numeric(indicator_model_mdd$fe_soln["covid_infections_sqrt"])]
dataset_anx[, indicator := social_mobility * as.numeric(indicator_model_anx$fe_soln["social_mobility"]) + covid_infections_sqrt * as.numeric(indicator_model_anx$fe_soln["covid_infections_sqrt"])]

dataset_mdd[, `:=` (age_int = m_mid_age * indicator, sex_int = m_percent_female * indicator, int_both_md = cv_both_md * indicator, int_cross_sectional = cv_cross_sectional * indicator, int_online_panel = cv_online_panel * indicator)]
dataset_anx[, `:=` (age_int = m_mid_age * indicator, sex_int = m_percent_female * indicator, int_both_md = cv_both_md * indicator, int_cross_sectional = cv_cross_sectional * indicator, int_online_panel = cv_online_panel * indicator)]

write.csv(dataset_mdd, "/FILEPATH/dataset_final_prepped_mdd.csv", row.names = F)
write.csv(dataset_anx, "/FILEPATH/dataset_final_prepped_anx.csv", row.names = F)

# Run final model ---------------------------------------------------------

potential_covs <- c('indicator', 'age_int', 'sex_int', "int_both_md", "int_cross_sectional", "int_online_panel", "cv_cross_panel")

trimming <- 0.95

final_model_table_mdd <- data.table(disorder = "mdd", cov = potential_covs)
final_model_table_mdd[, row_id := seq_len(.N)]

final_model_table_anx <- data.table(disorder = "anx", cov = potential_covs)
final_model_table_anx[, row_id := seq_len(.N)]

for(d in c("mdd", "anx")){
  mr_dataset <- MRData()
  if(d == "mdd"){
    mr_dataset$load_df(
      data = dataset_mdd,
      col_obs = "logit_dif", col_obs_se = "logit_dif_se",
      col_covs = as.list(c(potential_covs, "row_id")), col_study_id = "cohort" )
  } else if(d == "anx"){
    mr_dataset$load_df(
      data = dataset_anx,
      col_obs = "logit_dif", col_obs_se = "logit_dif_se",
      col_covs = as.list(c(potential_covs, "row_id")), col_study_id = "cohort" )
  }

    selected_covs <- potential_covs[potential_covs != 'indicator']

    cov_list <- list(LinearCovModel('intercept', use_re = T,  prior_beta_uniform = c(0, 0)), LinearCovModel('indicator', use_re = T))
    for(c in selected_covs){
      cov_list <- c(cov_list, list(LinearCovModel(c, use_re = F)))
    }

    model <- MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)

    model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
    
    betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
    betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]

    betas <- betas[,.(cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]
    
    log_lik <- -model$get_objective()
    aic_current <- -2 * log_lik + 2 * (length(model$beta_soln) + length(model$gamma_soln))
   
    worst_p <- betas[cov != "intercept", max(p) ]
    betas[p == worst_p & cov != "intercept", cov]
    betas[, disorder := d]

    assign(paste0('final_model_table_', d), merge(get(paste0("final_model_table_", d)), betas[,.(disorder, cov, beta_initial = round(beta,3), lower_initial = round(lower, 3), upper_initial = round(upper, 3), p_initial = round(p,4), aic_initial = round(aic_current, 1))], all.x = T, by = c("disorder", "cov")))

    stop <- F
    while(stop == F & worst_p > 0){
      previous_model <- model
      aic_previous <- aic_current
      drop_cv <- betas[p == worst_p & cov != "intercept", cov]
      print(paste0('dropping ', drop_cv))
      selected_covs <- selected_covs[selected_covs != drop_cv]

      cov_list <- list(LinearCovModel('intercept', use_re = T,  prior_beta_uniform = c(0, 0)), LinearCovModel('indicator', use_re = T))
      for(c in selected_covs){
        cov_list <- c(cov_list, list(LinearCovModel(c, use_re = F)))
      }

      model <- MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming)

      model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)

      log_lik <- -model$get_objective()
      aic_current <- -2 * log_lik + 2 * (length(model$beta_soln) + length(model$gamma_soln))
     
      if(aic_previous <= aic_current){
        print(selected_covs)
        stop <- T
        model <- previous_model
        aic_current <- aic_previous
      }

      betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
      betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]

      betas <- betas[,.(cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p )]

      betas
      worst_p <- betas[cov != "intercept", max(p)]

      if(worst_p < 0.05){stop <- T}
    }

    assign(paste0('final_model_', d), model)
    assign(paste0('final_model_table_', d), merge(get(paste0("final_model_table_", d)), betas[,.(disorder = d, cov, beta_final = round(beta,3), lower_final = round(lower,3), upper_final = round(upper,3), p_final = round(p, 4), aic_final = round(aic_current,1))], all.x = T, by = c("disorder", "cov")))

    infections_variable <- indicator_model_anx$cov_model_names[grepl("infections", indicator_model_anx$cov_model_names)]
    deaths_variable <- indicator_model_anx$cov_model_names[grepl("death", indicator_model_anx$cov_model_names)]

    if(d == 'mdd'){
      indicator_draws <- rnorm(1000, as.numeric(final_model_mdd$fe_soln["indicator"]), get_beta_sd(final_model_mdd)["indicator"])
      social_mobility_beta <- round(as.numeric(indicator_model_mdd$fe_soln["social_mobility"]) * as.numeric(final_model_mdd$fe_soln["indicator"]),1)
      social_mobility_draws <- rnorm(1000, as.numeric(indicator_model_mdd$fe_soln["social_mobility"]), get_beta_sd(indicator_model_mdd)["social_mobility"])
      social_mobility_draws <- social_mobility_draws * indicator_draws
      social_mobility_lower <- quantile(social_mobility_draws, 0.025)
      social_mobility_upper <- quantile(social_mobility_draws, 0.975)
      covid_infections_beta <- round(as.numeric(indicator_model_mdd$fe_soln[infections_variable]) * as.numeric(final_model_mdd$fe_soln["indicator"]),1)
      covid_infections_draws <- rnorm(1000, as.numeric(indicator_model_mdd$fe_soln[infections_variable]), get_beta_sd(indicator_model_mdd)[infections_variable])
      covid_infections_draws <- covid_infections_draws * indicator_draws
      covid_infections_lower <- quantile(covid_infections_draws, 0.025)
      covid_infections_upper <- quantile(covid_infections_draws, 0.975)
     } else if(d == 'anx'){
      indicator_draws <- rnorm(1000, as.numeric(final_model_anx$fe_soln["indicator"]), get_beta_sd(final_model_anx)["indicator"])
      social_mobility_beta <- round(as.numeric(indicator_model_anx$fe_soln["social_mobility"]) * as.numeric(final_model_anx$fe_soln["indicator"]),1)
      social_mobility_draws <- rnorm(1000, as.numeric(indicator_model_anx$fe_soln["social_mobility"]), get_beta_sd(indicator_model_anx)["social_mobility"])
      social_mobility_draws <- social_mobility_draws * indicator_draws
      social_mobility_lower <- quantile(social_mobility_draws, 0.025)
      social_mobility_upper <- quantile(social_mobility_draws, 0.975)
      covid_infections_beta <- round(as.numeric(indicator_model_anx$fe_soln[infections_variable]) * as.numeric(final_model_anx$fe_soln["indicator"]),1)
      covid_infections_draws <- rnorm(1000, as.numeric(indicator_model_anx$fe_soln[infections_variable]), get_beta_sd(indicator_model_anx)[infections_variable])
      covid_infections_draws <- covid_infections_draws * indicator_draws
      covid_infections_lower <- quantile(covid_infections_draws, 0.025)
      covid_infections_upper <- quantile(covid_infections_draws, 0.975)
    }

    indicator_row <- get(paste0('final_model_table_', d))[cov == 'indicator', ]
    social_mobility_row <- indicator_row[,.(disorder, cov = 'social_mobility', row_id = 1.1, beta_initial = "-", lower_initial = "-", upper_initial ="-", p_initial = "-", aic_initial, beta_final = social_mobility_beta, lower_final = social_mobility_lower, upper_final = social_mobility_upper, p_final = NA, aic_final)]
    covid_infections_row <- indicator_row[,.(disorder, cov = 'covid_infections_daily', row_id = 1.2, beta_initial = "-", lower_initial = "-", upper_initial ="-", p_initial = "-", aic_initial, beta_final = covid_infections_beta, lower_final = covid_infections_lower, upper_final = covid_infections_upper, p_final = NA, aic_final)]
    assign(paste0('final_model_table_', d), rbind(get(paste0('final_model_table_', d)), social_mobility_row, covid_infections_row))
    assign(paste0('final_model_table_', d), get(paste0('final_model_table_', d))[order(row_id)])

  py_save_object(object = model, filename = paste0("/FILEPATH/adj_model_", d, ".pkl"), pickle = "dill")
}
print('finished')

final_model_table_mdd
final_model_table_anx



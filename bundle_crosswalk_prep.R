###########################################################################################
#### Make changes to bundle, save new bundle version, and update bundle meta-data file ####
###########################################################################################

rm(list=ls())

library(data.table)
library(openxlsx)
library(msm)

release <- 9 # 9 = GBD 2021

acause <- 'mental_unipolar_mdd'

cause_meta_data <- data.table(acause_label = 'mental_unipolar_mdd', bundle_id = 10051)

bundle_id <- cause_meta_data[acause_label == acause, bundle_id]

## Load required functions
source("/FILEPATH/get_bundle_data.R")
source("/FILEPATH/upload_bundle_data.R")
source("/FILEPATH/save_bundle_version.R")
source("/FILEPATH/get_ids.R")
source("/FILEPATH/get_draws.R")
source("/FILEPATH/get_bundle_version.R")
source("/FILEPATH/get_crosswalk_version.R")
source("/FILEPATH/save_crosswalk_version.R")
source("/FILEPATH/get_covariate_estimates.R")
source("/FILEPATH/get_outputs.R")
source("/FILEPATH/get_population.R")

source("/FILEpATH/custom_functions.R")

rlogit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}
get_beta_vcov <- function(model){
  model_specs <- mr$core$other_sampling$extract_simple_lme_specs(model)
  beta_hessian <- mr$core$other_sampling$extract_simple_lme_hessian(model_specs)
  solve(beta_hessian)
}
get_beta_sd <- function(model){
  beta_sd <- sqrt(diag(get_beta_vcov(model)))
  names(beta_sd) <- model$cov_names
  return(beta_sd)
}

v_id <- 42865 

# Pull bundle version and set aside non-sex-split data for later crosswalk --------
## Pull dataset
review_sheet <- get_bundle_version(v_id, fetch = 'all')

length(unique(review_sheet$nid))
length(unique(review_sheet$group))
unique(review_sheet$location_name)
length(review_sheet$nid)
length(review_sheet[group_review == 1 & sex == "Both", unique(group)])
length(review_sheet[group_review == 1 & sex != "Both", unique(group)])
length(review_sheet[group_review == 1 & cv_antidep == 1, unique(group)])
length(review_sheet[group_review == 1 & cv_any_mental == 1, unique(group)])
length(review_sheet[group_review == 1 & cv_antidep == 0 & cv_any_mental == 0, unique(group)])

## Remove excluded estimates ##
review_sheet[is.na(group_review), group_review := 1]
review_sheet <- review_sheet[group_review == 1, ]
review_sheet <- review_sheet[year_start > 1979,]

review_sheet[, study_covariate := "ref"]

## Save mat vs any mental pairs for later
review_sheet[, pair_label := paste0(location_name, "_", sex, "_", age_start, "_", age_end, "_", year_start, "_", year_end, "_recall", cv_recall_1yr)]

covariates <- names(review_sheet)[names(review_sheet) %like% "cv_"]
covariates <- covariates[!(covariates %in% c("cv_mat", "cv_recall_1yr", "cv_wmhs", "cv_whs", "cv_psychiatrist_only", "cv_psychologist_only", "cv_mat_lenient", "cv_antidep_psychiatrist"))]

## Check added value of covariates
studies_cv_antidep <- unique(review_sheet[cv_antidep == 1, group])
studies_cv_any_mental <- unique(review_sheet[cv_any_mental== 1, group])

length(studies_cv_antidep[!(studies_cv_antidep %in% review_sheet[cv_antidep == 0, group])]) # cv_antidep provides 5 studies
length(studies_cv_any_mental[!(studies_cv_any_mental %in% review_sheet[cv_any_mental == 0, group])]) # cv_any_mental provides 10 studies

review_sheet[, covariates_applied := ""]
for(c in covariates){
  c_simple <- substring(c, 4)
  review_sheet[get(c) == 1, covariates_applied := paste0(covariates_applied, "-", c_simple)]
}
review_sheet[, covariates_applied := substring(covariates_applied,2)] # get rid of leading - symbol

x_walk_pairs <- data.table()
for(c in covariates){
  pairs <- merge(review_sheet[get(c) == 0, .(pair_label, location_id, mid_year = round((year_start + year_end) / 2), mean_ref = mean, se_ref = standard_error, ref = covariates_applied)],
                 review_sheet[get(c) == 1, .(pair_label, mean_alt = mean, se_alt = standard_error, alt = covariates_applied)], by = "pair_label")
  x_walk_pairs <- rbind(x_walk_pairs, pairs)
}

for(c in covariates){
  c_simple <- substring(c, 4)
  x_walk_pairs[grepl(c_simple, alt) & grepl(c_simple, ref), `:=` (alt = gsub(c_simple, "", alt), ref = gsub(c_simple, "", ref))]
}
x_walk_pairs[, `:=` (ref = gsub('\\--', '-', ref), alt = gsub('\\--', '-', alt))]
x_walk_pairs[, `:=` (ref = gsub('^\\-|\\-$', '', ref), alt = gsub('^\\-|\\-$', '', alt))]

x_walk_pairs[mean_ref > mean_alt, `:=` (ref = alt, mean_ref = mean_alt, se_ref = se_alt, alt = ref, mean_alt = mean_ref, se_alt = se_ref)]
x_walk_pairs <- unique(x_walk_pairs[mean_ref != 0 & mean_alt != 0])

x_walk_pairs[, row_id := seq_len(.N)]
x_walk_pairs[, `:=` (log_r = log(mean_ref), log_r_se = deltamethod(~log(x1), mean_ref, se_ref^2),
                     log_a = log(mean_alt), log_a_se = deltamethod(~log(x1), mean_alt, se_alt^2),
                     logit_r = logit(mean_ref), logit_r_se = deltamethod(~log(x1/(1-x1)), mean_ref, se_ref^2),
                     logit_a = logit(mean_alt), logit_a_se = deltamethod(~log(x1/(1-x1)), mean_alt, se_alt^2)), by = "row_id"]

x_walk_pairs[, `:=` (log_diff = log_a - log_r, log_diff_se = sqrt(log_r_se^2 + log_a_se^2),
                     logit_diff = logit_a - logit_r, logit_diff_se = sqrt(logit_r_se^2 + logit_a_se^2))]

haqi <- get_covariate_estimates(covariate_id = 1099, location_id = unique(x_walk_pairs$location_id), year_id = unique(x_walk_pairs$mid_year), release_id = release)
x_walk_pairs <- merge(x_walk_pairs, haqi[,.(location_id, haqi = mean_value, mid_year = year_id)], by = c("location_id", "mid_year"))
x_walk_pairs[ref == "", ref := "mat"]
x_walk_pairs[alt == "", alt := "mat"]
x_walk_pairs[, pair_type := ifelse(ref < alt, paste(ref, alt), paste(alt, ref))]
pair_count <- data.table(table(x_walk_pairs$pair_type))
for(c in covariates){
  c_simple <- gsub("cv_", "", c)
  x_walk_pairs[, paste0(c) := 0]
  x_walk_pairs[alt == c_simple, paste0(c) := 1]
  x_walk_pairs[ref == c_simple, paste0(c) := -1]
  pair_count[grepl(c_simple, V1), paste0(c_simple) := 1 * N]
  pair_count[!grepl(c_simple, V1), paste0(c_simple) := 0]
}

pair_count
plot(rbind(x_walk_pairs[ref == "mat" & alt == "any_mental",.(haqi, val = exp(log_diff))], x_walk_pairs[ref == "any_mental" & alt == "mat",.(haqi, val = 1/exp(log_diff))]))

plot(rbind(x_walk_pairs[ref == "mat" & alt == "any_mental",.(haqi, val = log_diff)], x_walk_pairs[ref == "any_mental" & alt == "mat",.(haqi, val = -log_diff)]))

plot(rbind(x_walk_pairs[ref == "mat" & alt == "any_mental",.((haqi)^3, val = log_diff)], x_walk_pairs[ref == "any_mental" & alt == "mat",.((haqi)^3, val = -log_diff)]))



plot(rbind(x_walk_pairs[ref == "mat" & alt == "antidep",.(haqi, val = exp(log_diff))], x_walk_pairs[ref == "antidep" & alt == "mat",.(haqi, val = 1/exp(log_diff))]))

## Inspecting haqi variation to determine whether to test haqi interaction

x_walk_pairs[, haqi_c_100 := haqi - 100]
x_walk_pairs[, `:=` (haqi_x_any_mental = cv_any_mental * haqi_c_100, haqi_x_antidep = cv_antidep * haqi_c_100)]

x_walk_pairs[, `:=` (haqi_c_log = log(haqi) - log(100), haqi_c_sqrt = sqrt(haqi) - sqrt(100), haqi_c_cr = (haqi)^(1/3) - 100^(1/3), haqi_c_squared = haqi^2 - 100^2, haqi_c_cubed = haqi^3 - 100^3)]

x_walk_pairs[, `:=` (haqi_x_any_mental_log = cv_any_mental * haqi_c_log, haqi_x_antidep_log = cv_antidep * haqi_c_log)]
x_walk_pairs[, `:=` (haqi_x_any_mental_sqrt = cv_any_mental * haqi_c_sqrt, haqi_x_antidep_sqrt = cv_antidep * haqi_c_sqrt)]
x_walk_pairs[, `:=` (haqi_x_any_mental_cr = cv_any_mental * haqi_c_cr, haqi_x_antidep_cr = cv_antidep * haqi_c_cr)]
x_walk_pairs[, `:=` (haqi_x_any_mental_squared = cv_any_mental * haqi_c_squared, haqi_x_antidep_squared = cv_antidep * haqi_c_squared)]
x_walk_pairs[, `:=` (haqi_x_any_mental_cubed = cv_any_mental * haqi_c_cubed, haqi_x_antidep_cubed = cv_antidep * haqi_c_cubed)]

# Write study-characteristics table ---------------------------------------
study_characteristics_table <- unique(review_sheet[, .(Citation = field_citation_value, Location = location_name, Sex = sex, `Female cases (%)` = paste0(round(prop_female*100, 1)), Years = paste0(year_start, " to ", year_end), Ages = paste0(age_start, " to ", age_end), Estimate = round(mean*100, 1), `Lower CI` = round(lower*100, 1), `Upper CI` = round(upper*100, 1), `Total cases` = sample_size, `Treated cases` = cases, Recall = paste0(recall_type_value, " ", gsub('Period: ', "", recall_type)), `Service type` = case_name)])
study_characteristics_table[, Citation := gsub("<i>", "", Citation)]
study_characteristics_table[, Citation := gsub("</i>", "", Citation)]
study_characteristics_table[`Service type` == "mat", `Service type` := "Minimally adequate treatment"]
study_characteristics_table[`Female cases (%)` == "NA", `Female cases (%)` := "Not reported"]
write_table <- F
if(write_table){
  write.xlsx(study_characteristics_table, "/FILEPATH/study_characteristics.xlsx")
}

# Run models --------------------------------------------------------------
library(reticulate)
reticulate::use_python("/FILEPATH/python")
mr <- import("mrtool")

trimming <- 0.95

#fractional_polys_test <- c("haqi_x_any_mental_log", "haqi_x_antidep_log",
#                           "haqi_x_any_mental_sqrt", "haqi_x_antidep_sqrt",
#                           "haqi_x_any_mental_cr", "haqi_x_antidep_cr", 
#                           "haqi_x_any_mental_squared", "haqi_x_antidep_squared", 
#                           "haqi_x_any_mental_cubed", "haqi_x_antidep_cubed")

## Model 1:
mr_dataset <- mr$MRData()
mr_dataset$load_df(
  data = x_walk_pairs, col_obs = "logit_diff", col_obs_se = "logit_diff_se", # no trim p = 0.1035, 10% trim p = 0.0514
  #col_covs = as.list(c(covariates, "haqi_x_any_mental", "haqi_x_antidep", fractional_polys_test)), col_study_id = "pair_label" )
  col_covs = as.list(c(covariates, "haqi_x_any_mental", "haqi_x_antidep")), col_study_id = "pair_label" )

cov_list <- list()
for(c in c(covariates, "haqi_x_any_mental", "haqi_x_antidep")){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
#for(c in c(covariates, "haqi_x_any_mental_log", "haqi_x_antidep_log")){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
#for(c in c(covariates, "haqi_x_any_mental_sqrt", "haqi_x_antidep_sqrt")){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
#for(c in c(covariates, "haqi_x_any_mental_cr", "haqi_x_antidep_cr")){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
#for(c in c(covariates, "haqi_x_any_mental_squared", "haqi_x_antidep_squared")){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
#for(c in c(covariates, "haqi_x_any_mental_cubed", "haqi_x_antidep_cubed")){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}

model_1 <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming) 

model_1$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)

aic_m1 <- mrbrt_aic(model_1)
bic_m1 <- mrbrt_bic(model_1)

betas <- data.table(cov = model_1$cov_names, coef = as.numeric(model_1$beta_soln), se = get_beta_sd(model_1))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas_1 <- betas[,.(Step = "Step 1", cov, odds = round(exp(coef),3), odds_lower = round(exp(lower),3), odds_upper = round(exp(upper),3),beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = signif(p,4), aic = aic_m1, bic = bic_m1)]

betas_1

  ## Model 2:
  cov_list <- list()
  for(c in c(covariates, "haqi_x_any_mental")){cov_list <- c(cov_list, list(mr$LinearCovModel(c, use_re = F)))}
  
  model_2 <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming) # Betas stabalise at 5% trimming so 10% not necessary. Inspection of outliers supports their exclusion.
  
  model_2$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
  
  aic_m2 <- mrbrt_aic(model_2)
  bic_m2 <- mrbrt_bic(model_2)
  
  betas <- data.table(cov = model_2$cov_names, coef = as.numeric(model_2$beta_soln), se = get_beta_sd(model_2))
  betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
  betas_2 <- betas[,.(Step = "Step 2", cov, odds = round(exp(coef),3), odds_lower = round(exp(lower),3), odds_upper = round(exp(upper),3),beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = signif(p,4), aic = aic_m2, bic = bic_m2)]
  betas_2
  
  aic_m1
  aic_m2
  
  bic_m1
  bic_m2

  model <- model_2


model_results <- rbind(betas_1, betas_2)

write_model <- F
if(write_model){
  write.csv(model_results, "/FILEPATH/model_results.csv", row.names = F)
}
  
used_data <- data.table(cbind(model$data$to_df(), data.frame(w = model$w_soln)))
used_data[w<1]

# Age-sex-split data where available prior to crosswalking ----------------
# Remove overall data from studies with sex-specific data that were only retained to inform crosswalks
review_sheet <- review_sheet[!(group == 15 & sex == "Both" & age_start == 30 & age_end == 99),]
review_sheet <- review_sheet[!((nid %in% c(120059, 120058) & age_start == 16 & age_end > 84)),]

review_sheet[, unique_study_label := paste0(group, location_id, year_start, year_end, cv_antidep, cv_any_mental)]

review_sheet[, `:=` (min_age = min(age_start), max_age = max(age_end)), by = "unique_study_label"]
review_sheet[, `:=` (widest_age = ifelse(min_age == age_start & max_age == age_end, 1, 0))]
review_sheet_bysex_allage <- review_sheet[sex != "Both" & widest_age == 1,]
review_sheet_bothsex_byage <- review_sheet[sex == "Both" & widest_age == 0 & (unique_study_label %in% review_sheet_bysex_allage$unique_study_label),]

age_sex_split <- function(x){
  cases_male <- review_sheet_bysex_allage[unique_study_label == x & sex == "Male", cases] / (review_sheet_bysex_allage[unique_study_label == x & sex != "Both", sum(cases)])
  cases_female <- review_sheet_bysex_allage[unique_study_label == x & sex == "Female", cases] / (review_sheet_bysex_allage[unique_study_label == x & sex != "Both", sum(cases)])
  cases_male_se <- sqrt(cases_male*(1-cases_male)/(review_sheet_bysex_allage[unique_study_label == x & sex != "Both", sum(cases)]))
  cases_female_se <- sqrt(cases_female*(1-cases_female)/(review_sheet_bysex_allage[unique_study_label == x & sex != "Both", sum(cases)]))
  sample_male <- review_sheet_bysex_allage[unique_study_label == x & sex == "Male", sample_size] / (review_sheet_bysex_allage[unique_study_label == x & sex != "Both", sum(sample_size)])
  sample_female <- review_sheet_bysex_allage[unique_study_label == x & sex == "Female", sample_size] / (review_sheet_bysex_allage[unique_study_label == x & sex != "Both", sum(sample_size)])
  sample_male_se <- sqrt(sample_male*(1-sample_male)/(review_sheet_bysex_allage[unique_study_label == x & sex != "Both", sum(sample_size)]))
  sample_female_se <- sqrt(sample_female*(1-sample_female)/(review_sheet_bysex_allage[unique_study_label == x & sex != "Both", sum(sample_size)]))
  ratio_male <- cases_male/sample_male
  ratio_female <- cases_female/sample_female
  ratio_male_se <- sqrt((cases_male^2/sample_male^2)*(cases_male_se^2/cases_male^2 + sample_male_se^2/sample_male^2))
  ratio_female_se <- sqrt((cases_female^2/sample_female^2)*(cases_female_se^2/cases_female^2 + sample_female_se^2/sample_female^2))
  male_by_age <- copy(review_sheet_bothsex_byage[unique_study_label == x])
  female_by_age <- copy(review_sheet_bothsex_byage[unique_study_label == x])
  male_by_age[,`:=`(mean=mean*ratio_male, standard_error=sqrt(standard_error^2*ratio_male_se^2 + standard_error^2*ratio_male^2 + ratio_male_se^2*mean^2), sex = "Male", prop_female = 0)]
  female_by_age[,`:=`(mean=mean*ratio_female, standard_error=sqrt(standard_error^2*ratio_female_se^2 + standard_error^2*ratio_female^2 + ratio_female_se^2*mean^2), sex = "Female", prop_female = 1)]
  by_age <- rbind(male_by_age, female_by_age)
  rm(male_by_age, female_by_age)
  by_age[, `:=` (study_covariate = "sex", crosswalk_parent_seq = seq, seq = NA, sample_size = NA, cases = NA)]
  return(by_age)
}

age_sex_split_list <- lapply(unique(review_sheet_bothsex_byage$unique_study_label), age_sex_split)
age_sex_split_list <- Reduce(function(x, y){rbind(x, y)}, age_sex_split_list)

review_sheet <- rbind(review_sheet[!(unique_study_label %in% age_sex_split_list$unique_study_label),], age_sex_split_list, fill = T)

# Crosswalk data ----------------------------------------------------------
review_sheet[, `:=` (mid_year = round((year_start + year_end) / 2), mean_logit = logit(mean), se_logit = deltamethod(~log(x1/(1-x1)), mean, standard_error^2)), by = c("mean", "standard_error")]

if(test_haqi){
  haqi <- get_covariate_estimates(covariate_id = 1099, location_id = unique(review_sheet$location_id), year_id = unique(review_sheet$mid_year), release_id = release)
  
  review_sheet <- merge(review_sheet, haqi[,.(location_id, mid_year = year_id, haqi = mean_value)], by = c("location_id", "mid_year"), all.x = T)
  review_sheet[, haqi_c_100 := haqi - 100]
  review_sheet[, `:=` (haqi_x_any_mental = cv_any_mental * haqi_c_100)]

  predict_matrix <- unique(review_sheet[,(c(covariates, "haqi_x_any_mental")), with = F])
  dat_pred1 <- mr$MRData()
  dat_pred1$load_df(
    data = predict_matrix,
    col_covs=as.list(c(covariates, "haqi_x_any_mental"))
  )
} else {
  predict_matrix <- unique(review_sheet[,(c(covariates)), with = F])
  dat_pred1 <- mr$MRData()
  dat_pred1$load_df(
    data = predict_matrix,
    col_covs=as.list(c(covariates)))
}

beta_samples <- mr$core$other_sampling$sample_simple_lme_beta(1000L, model)
gamma_outer_samples <- matrix(rep(model$gamma_soln, each = 1000L), nrow = 1000L)
draws <- model$create_draws(dat_pred1,
                            beta_samples = beta_samples,
                            gamma_samples = gamma_outer_samples,
                            random_study = T)

draws <- data.table(draws)
names(draws) <- paste0("draw_", as.numeric(gsub("V", "", names(draws)))-1)
predict_matrix <- cbind(predict_matrix, draws)
predict_matrix[, row_id := seq_len(.N)]
predict_matrix <- melt.data.table(predict_matrix, id.vars = names(predict_matrix)[!(names(predict_matrix) %like% "draw")], value.name="crosswalk", variable.name="draw")
predict_matrix[, `:=` (crosswalk = mean(crosswalk), crosswalk_se = sd(crosswalk)), by = 'row_id']
predict_matrix[, draw := NULL]
predict_matrix <- unique(predict_matrix)

if(test_haqi){
  review_sheet <- merge(review_sheet, predict_matrix, by = c(covariates, "haqi_x_any_mental"), all.x = T)
} else {
  review_sheet <- merge(review_sheet, predict_matrix, by = c(covariates), all.x = T)
}


review_sheet[, unique_label := paste0(group, "_", location_name, "_Both_", age_start, "_", age_end, "_", year_start, "_", year_end)] # Set sex as "Both" in label to exclude duplicative estimates flagged in x_walk_pairs

review_sheet[crosswalk != 0, `:=` (mean_logit = mean_logit - crosswalk, se_logit = sqrt(se_logit^2 + crosswalk_se^2), sample_size = NA, cases = NA)]
for(c in covariates){review_sheet[get(c) == 1, study_covariate := gsub("cv_", "", c)]}

# Estimate sex-ratios -----------------------------------------------------

## Create paired dataset where each row is a sex pair ##
review_sheet[, mid_year := round((year_start + year_end)/2)]
review_sheet[, mid_year_for_prev := ifelse(mid_year < 1990, 1990, mid_year)]

prevalence <- data.table()
for(y in sort(unique(review_sheet[is.na(prop_female), mid_year_for_prev]))){
  prev_temp <- get_draws(year_id = y, gbd_id_type = "cause_id", gbd_id = 568, source = "como", release_id = 9,  metric_id = 3, sex_id = c(1, 2), measure_id = 5, age_group_id = 22, location_id = unique(review_sheet[mid_year_for_prev== y & is.na(prop_female), location_id]))
  
  prevalence <- rbind(prevalence, prev_temp)
  rm(prev_temp)
  print(paste0("Finished year ", y, " of ", max(review_sheet[is.na(prop_female), mid_year_for_prev])))
}
prevalence <- melt.data.table(prevalence, id.vars = names(prevalence)[!(names(prevalence) %like% "draw")], value.name="prev", variable.name="draw")
population <- get_population(age_group_id = 22, sex_id = c(1, 2), location_id = unique(prevalence$location_id), year_id = unique(prevalence$year_id), release_id = release)
prevalence <- merge(prevalence, population[,.(location_id, year_id, sex_id, population)], all.x = T, by = c("location_id", "year_id", "sex_id"))
prevalence[, cases := prev * population]

prevalence <- merge(prevalence[sex_id == 1, .(location_id, year_id, draw, cases_m = cases)],
                    prevalence[sex_id == 2, .(location_id, year_id, draw, cases_f = cases)], by = c("location_id", "year_id", "draw"))
prevalence[, prop_female := cases_f / (cases_m + cases_f)]
prevalence[, prop_female_mean := mean(prop_female), by = c("location_id", "year_id")]
prevalence <- unique(prevalence[,.(location_id, mid_year_for_prev = year_id, prop_female_impute = prop_female_mean)])

review_sheet <- merge(review_sheet, all.x = T, prevalence, by = c("location_id", "mid_year_for_prev"))
review_sheet[sex == "Male", prop_female := 0]
review_sheet[sex == "Female", prop_female := 1]
review_sheet[is.na(prop_female), prop_female := prop_female_impute]
review_sheet[, prop_female_impute := NULL]
review_sheet[, unique_label := paste0(group, "_", location_name, "_", year_start, "_", year_end)]

# Remove non-MAT estimates that have an equivalent MAT estimate -----------
review_sheet[, mat_unique_label := paste0(group, "_", location_name, "_", sex, "_", age_start, "_", age_end, "_", year_start, "_", year_end, "_recall", cv_recall_1yr)]

mat_unique_labels <- review_sheet[cv_antidep == 0 & cv_any_mental == 0 & mean > 0, unique(mat_unique_label)]  # don't exclude alternative data for 0 MAT estimates
mat_unique_labels_0_mat <- review_sheet[cv_antidep == 0 & cv_any_mental == 0 & mean == 0, unique(mat_unique_label)] # don't exclude alternative data for 0 MAT estimates

mat_unique_labels

# Exception made for Nigerian estimate where MAT estimate is 0 but there is an available non-zero but crosswalked any mental health service estimate. 
review_sheet <- review_sheet[!(!( cv_antidep == 0 & cv_any_mental == 0) & (mat_unique_label %in% mat_unique_labels)),]

review_sheet <- review_sheet[!(mat_unique_label == "36_Nigeria_Both_18_100_2002_2003_recall1" & cv_any_mental != 1),] # only non-zero estimate in Nigeria is a crosswalked any-mental health service

# Run MR-BRT sex model ----------------------------------------------------
mr_dataset <- mr$MRData()
mr_dataset$load_df(
  data = review_sheet[group_review == 1 & mean_logit != -Inf,], col_obs = "mean_logit", col_obs_se = "se_logit",
  col_covs = list("prop_female"), col_study_id = "unique_label" )
cov_list <- c(list(mr$LinearCovModel('intercept', use_re = T)), list(mr$LinearCovModel("prop_female", use_re = F)))

model <- mr$MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =trimming) 

model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)

betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas <- betas[,.(cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), odds = round(exp(coef), 2), odds_lower = round(exp(lower), 2), odds_upper = round(exp(upper), 2), z = round(z, 3), p = signif(p,4))]
betas

write_model <- F
if(write_model){
  py_save_object(object = model, filename = "/FILEPATH/sex_model.pkl", pickle = "dill")
  write.csv(betas, "/FILEPATH/sex_betas.csv", row.names = F)
}

## Inspect data and model fit
used_data <- data.table(cbind(model$data$to_df(), data.frame(w = model$w_soln)))
df_pred1 <- data.frame(cv_mat_lenient = 1, prop_female = seq(0, 1, by = 0.01))
dat_pred1 <- mr$MRData()
dat_pred1$load_df(
  data = df_pred1,
  col_covs=c(list('prop_female'))
)
pred1 <- model$predict(data = dat_pred1)
df_pred1$pred1 <- pred1
used_data[, trimmed := as.factor(ifelse(w < 1, "yes", "no"))]
with(used_data[w == 1], plot(prop_female, obs))
with(used_data, plot(prop_female, obs, col = c("black", "red")[used_data$trimmed]))
with(df_pred1, lines(prop_female, pred1))

predict_matrix <- data.table(prop_female = unique(c(0, 1, review_sheet$prop_female)))
dat_pred1 <- mr$MRData()
dat_pred1$load_df(
  data = predict_matrix,
  col_covs=c(list('prop_female'))
)

beta_samples <- mr$core$other_sampling$sample_simple_lme_beta(1000L, model)
gamma_outer_samples <- matrix(rep(model$gamma_soln, each = 1000L), nrow = 1000L)
draws <- model$create_draws(dat_pred1,
                            beta_samples = beta_samples,
                            gamma_samples = gamma_outer_samples,
                            random_study = T)

draws <- data.table(draws)
names(draws) <- paste0("draw_", as.numeric(gsub("V", "", names(draws)))-1)
predict_matrix <- cbind(predict_matrix, draws)

predict_matrix <- melt.data.table(predict_matrix, id.vars = names(predict_matrix)[!(names(predict_matrix) %like% "draw")], value.name="service_cov", variable.name="draw")
predict_matrix <- merge(predict_matrix, predict_matrix[prop_female == 0, .(draw, service_cov_male = service_cov)], by = c("draw"))
predict_matrix <- merge(predict_matrix, predict_matrix[prop_female == 1, .(draw, service_cov_female = service_cov)], by = c("draw"))
predict_matrix[, `:=` (dif_male = service_cov - service_cov_male, dif_female = service_cov - service_cov_female)]

predict_matrix <- predict_matrix[,.(dif_male = mean(dif_male), dif_male_se = sd(dif_male), dif_female = mean(dif_female), dif_female_se = sd(dif_female)), by = c("prop_female")]

review_sheet <- merge(review_sheet, predict_matrix, by =c("prop_female"))

review_sheet_male <- review_sheet[sex == "Both",]
review_sheet_female <- review_sheet[sex == "Both",]
review_sheet <- review_sheet[sex != "Both"]
review_sheet_male[, `:=` (mean_logit = mean_logit - dif_male, se_logit = sqrt(se_logit^2 + dif_male_se^2), sex = "Male", study_covariate = ifelse(study_covariate == "ref", "sex", paste0(study_covariate, ", sex")), crosswalk_parent_seq = seq, seq = NA, sample_size = NA, cases = NA, lower = NA, upper = NA, uncertainty_type_value = NA)]
review_sheet_female[, `:=` (mean_logit = mean_logit - dif_female, se_logit = sqrt(se_logit^2 + dif_female_se^2), sex = "Female", study_covariate = ifelse(study_covariate == "ref", "sex", paste0(study_covariate, ", sex")), crosswalk_parent_seq = seq, seq = NA, sample_size = NA, cases = NA, lower = NA, upper = NA, uncertainty_type_value = NA)]
review_sheet <- rbind(review_sheet, review_sheet_male, review_sheet_female)

review_sheet[study_covariate != "ref" & mean_logit != -Inf, `:=` (crosswalk_parent_seq = seq, mean = rlogit(mean_logit), standard_error = deltamethod(~exp(x1)/(1+exp(x1)), mean_logit, se_logit^2), lower = NA, upper = NA, uncertainty_type_value = NA), by = c("mean_logit", "se_logit")]


# Prepaire for upload -----------------------------------------------------
review_sheet[, (c("unique_study_label", "mat_unique_label", "unique_label", "dif_male", "dif_male_se", "dif_female", "dif_female_se", "crosswalk", "crosswalk_se", "row_id", "haqi", "haqi_c_100", "mean_logit", "se_logit", "min_age", "max_age", "widest_age", "pair_label", "covariates_applied")) := NULL]
duplicate_columns <- data.table(table(tolower(names(review_sheet))))[N>1, V1]
if(length(duplicate_columns) > 0){review_sheet[,(duplicate_columns) := NULL]} # duplicate colums can exist in database for some bundles and causes upload issues

# Get rid of special characters that seem to appear sometimes from the database
non_proper_char <- c("Ã", "ƒ", "Æ", "’", "‚", "Â", "¢", "â", "¬", "Å", "¡", "¾", "†", "€", "™", "„", "š", "ž", "¦", "…", "œ")
for(n in non_proper_char){
  review_sheet[, `:=` (site_memo = gsub(n, "", site_memo))]
  for(note in names(review_sheet)[names(review_sheet) %like% "note_"]){
    review_sheet[, paste0(note) := gsub(n, "", get(note))]
    review_sheet[nchar(get(note)) > 1999, paste0(note) := substring(get(note), 1, 1999)]
  }
}

# Exclude 20% for cross-validation testing if doing s ---------------------------
cross_validation <- F
if(cross_validation){
  source("/FILEPATH/get_location_metadata.R")
  locatins <- get_location_metadata(9, release_id = 9)
  
  set.seed(212015452) ## mashed numpad with fist 3 times 16/9/2024 - Damian
  
  review_sheet[, country := location_name]
  
  review_sheet[location_name == "Beijing", country := "China"]
  review_sheet[location_name == "Shanghai", country := "China"]
  review_sheet[location_name == "New York", country := "United States of America"]
  review_sheet[location_name == "Maryland", country := "United States of America"]
  review_sheet[location_name == "São Paulo", country := "Brazil"]
  review_sheet[location_name == "Minas Gerais", country := "Brazil"]
  
  countries_to_keep <- sample(unique(review_sheet$country), 31*0.8, replace=F)
  
  review_sheet <- review_sheet[country %in% countries_to_keep,]
}



crosswalk_save_folder <- paste0("FILEPATH/", acause, "/", bundle_id, "/FILEPATH/")
dir.create(file.path(crosswalk_save_folder), showWarnings = FALSE)
crosswalk_save_file <- paste0(crosswalk_save_folder, "crosswalk_mat_", gsub(":", "-", gsub(" ", "-", Sys.time())), ".xlsx")
write.xlsx(review_sheet, crosswalk_save_file, sheetName = "extraction")

## Upload crosswalked dataset to database
save_results <- save_crosswalk_version(v_id, crosswalk_save_file, description = paste0("LABEL", Sys.time()))


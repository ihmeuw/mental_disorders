
source("/FILEPATH/get_ids.R")
source("/FILEPATH/get_draws.R")
source("/FILEPATH/interpolate.R")
source("/FILEPATH/get_location_metadata.R")
source("/FILEPATH/get_population.R")
library(mrbrt001, lib.loc = '/FILEPATH/')
library(openxlsx)
library(data.table)

get_beta_sd_old <- function(model){
  model_specs <- mrbrt001::core$other_sampling$extract_simple_lme_specs(model)
  beta_hessian <- mrbrt001::core$other_sampling$extract_simple_lme_hessian(model_specs)
  beta_sd <- 1/sqrt(diag(beta_hessian))
  names(beta_sd) <- model$cov_names
  return(beta_sd)
}

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



rlogit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}

epi_age_groups <- c(2:3, 388:389, 238, 34,  6:20, 30:32, 235)

model_mdd <- py_load_object(filename = "/FILEPATH/adj_model_mdd.pkl", pickle = "dill")
model_anx <- py_load_object(filename = "/FILEPATH/adj_model_anx.pkl", pickle = "dill")

# Load datasets ------------------------------------------------------------

dataset_mdd <- fread("/FILEPATH/dataset_final_prepped_mdd.csv")
dataset_anx <- fread("/FILEPATH/dataset_final_prepped_anx.csv")

#### Put subnational locations to countries
table(unique(dataset_mdd[,.(location_name, cohort)])$location_name)
table(dataset_mdd$location_name)

table(unique(dataset_anx[,.(location_name, cohort)])$location_name)

dataset_mdd[location_id == 35436, `:=` (location_id = 67, location_name = 'Japan|JPN')]
dataset_mdd[location_id == 354, `:=` (location_id = 6, location_name = 'China|CHN')]
dataset_mdd[location_id == 491, `:=` (location_id = 6, location_name = 'China|CHN')]
dataset_mdd[location_id == 491, `:=` (location_id = 6, location_name = 'China|CHN')]
dataset_mdd[location_id == 4749, `:=` (location_id = 95, location_name = 'United Kingdom|GBR')]

dataset_anx[location_id == 35436, `:=` (location_id = 67, location_name = 'Japan|JPN')]
dataset_anx[location_id == 354, `:=` (location_id = 6, location_name = 'China|CHN')]
dataset_anx[location_id == 491, `:=` (location_id = 6, location_name = 'China|CHN')]
dataset_anx[location_id == 491, `:=` (location_id = 6, location_name = 'China|CHN')]
dataset_anx[location_id == 4749, `:=` (location_id = 95, location_name = 'United Kingdom|GBR')]

dataset_mdd[, `:=` (linear_time = (date_covid_end + date_covid_start)/2 - (date_ref_end + date_ref_start)/2)]
dataset_mdd[, linear_time := linear_time / 365.25]

dataset_anx[, `:=` (linear_time = (date_covid_end + date_covid_start)/2 - (date_ref_end + date_ref_start)/2)]
dataset_anx[, linear_time := linear_time / 365.25]

used_data_mdd <- data.table(cbind(model_mdd$data$to_df(), data.frame(w = model_mdd$w_soln)))
used_data_anx <- data.table(cbind(model_anx$data$to_df(), data.frame(w = model_anx$w_soln)))

dataset_mdd <- merge(dataset_mdd, used_data_mdd[,.(w, row_id)])
dataset_anx <- merge(dataset_anx, used_data_anx[,.(w, row_id)])

# Run sensitivity analyses ------------------------------------------------
rm(validity_table)

for(d in c("mdd", "anx")){
  for(l in unique(get(paste0("dataset_", d))[, location_id])){
    if(d == "mdd"){
      dataset_subset <- dataset_mdd[w == 1]
      dataset_validity <- dataset_mdd[location_id == l,]
      selected_covs <- model_mdd$cov_model_names
    } else if(d == "anx"){
      dataset_subset <- dataset_anx[w == 1]
      dataset_validity <- dataset_anx[location_id == l,]
      selected_covs <- model_anx$cov_model_names
    }

    selected_covs <- selected_covs[!(selected_covs %in% c("intercept", "indicator"))]

    mr_dataset <- MRData()
    mr_dataset$load_df(
      data = dataset_subset[location_id != l,],
      col_obs = "logit_dif", col_obs_se = "logit_dif_se",
      col_covs = as.list(c("indicator", selected_covs, "row_id")), col_study_id = "cohort" )

    ##### leave-one-one original model
    cov_list <- list(LinearCovModel('intercept', use_re = T,  prior_beta_uniform = c(0, 0)), LinearCovModel('indicator', use_re = T))

    for(c in selected_covs){
      cov_list <- c(cov_list, list(LinearCovModel(c, use_re = F)))
    }

    model <- MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =1)

    model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)

    betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
    betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]

    betas <- betas[,.(cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]

    predict_matrix <- unique(dataset_subset[,.(location_id, indicator, age_int, sex_int, int_both_md, int_cross_sectional, int_online_panel, cv_cross_panel)])

    predict_data <- MRData()
    predict_data$load_df(data = predict_matrix, col_covs=as.list(names(predict_matrix)))

    beta_samples <- mrbrt001::core$other_sampling$sample_simple_lme_beta(1000L, model)
    gamma_outer_samples <- matrix(rep(model$gamma_soln, each = 1000L), nrow = 1000L)
    gamma_outer_samples[,1] <- 0 # remove intercept gamma
    draws <- model$create_draws(predict_data,
                                    beta_samples = beta_samples,
                                    gamma_samples = gamma_outer_samples,
                                    random_study = F)
    draws <- data.table(draws)
    names(draws) <- paste0("draw_", as.numeric(gsub("V", "", names(draws)))-1)
    predict_matrix <- cbind(predict_matrix, draws)

    predict_matrix <- melt.data.table(predict_matrix, id.vars = names(predict_matrix)[!(names(predict_matrix) %like% "draw")], value.name="adjustment", variable.name="draw")
    predict_matrix[, `:=` (adjustment = mean(adjustment), adjustment_sd = sd(adjustment, 0.025)), by = c("indicator", "age_int", "sex_int", "int_both_md", "int_cross_sectional", "int_online_panel", "cv_cross_panel")]
    predict_matrix[, draw := NULL]
    predict_matrix <- unique(predict_matrix)

    ##### leave-one-one linear time model

    if(d == "mdd"){
      i_selected_covs <- c("cv_cross_panel")
    } else {
      i_selected_covs <- c("cv_both_md", "cv_cross_panel")
    }

    mr_dataset_lt <- MRData()
    mr_dataset_lt$load_df(
      data = dataset_subset[location_id != l,],
      col_obs = "logit_dif", col_obs_se = "logit_dif_se",
      col_covs = as.list(c(i_selected_covs, "linear_time", "row_id")), col_study_id = "cohort" )

    cov_list <- list(LinearCovModel('intercept', use_re = T), LinearCovModel('linear_time', use_re = T))

    for(c in i_selected_covs){
      cov_list <- c(cov_list, list(LinearCovModel(c, use_re = F)))
    }

    model_lt <- MRBRT(data = mr_dataset_lt, cov_models =cov_list, inlier_pct =1)

    model_lt$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)

    betas <- data.table(cov = model_lt$cov_names, coef = as.numeric(model_lt$beta_soln), se = get_beta_sd(model_lt))
    betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]

    betas <- betas[,.(cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]

    predict_matrix_lt <- unique(dataset_subset[,.(location_id, linear_time, m_mid_age, m_percent_female, cv_both_md, cv_cross_panel)])

    predict_data_lt <- MRData()
    predict_data_lt$load_df(data = predict_matrix_lt, col_covs=as.list(names(predict_matrix_lt)))

    beta_samples <- mrbrt001::core$other_sampling$sample_simple_lme_beta(1000L, model_lt)
    gamma_outer_samples <- matrix(rep(model_lt$gamma_soln, each = 1000L), nrow = 1000L)
    gamma_outer_samples[,1] <- 0 # remove intercept gamma
    draws <- model_lt$create_draws(predict_data_lt,
                                  beta_samples = beta_samples,
                                  gamma_samples = gamma_outer_samples,
                                  random_study = F)
    draws <- data.table(draws)
    names(draws) <- paste0("draw_", as.numeric(gsub("V", "", names(draws)))-1)
    predict_matrix_lt <- cbind(predict_matrix_lt, draws)

    predict_matrix_lt <- melt.data.table(predict_matrix_lt, id.vars = names(predict_matrix_lt)[!(names(predict_matrix_lt) %like% "draw")], value.name="adjustment_lt", variable.name="draw")
    predict_matrix_lt[, `:=` (adjustment_lt = mean(adjustment_lt), adjustment_lt_sd = sd(adjustment_lt, 0.025)), by = c("linear_time", "m_mid_age", "m_percent_female", "cv_both_md", "cv_cross_panel")]
    predict_matrix_lt[, draw := NULL]
    predict_matrix_lt <- unique(predict_matrix_lt)

    ##### Merge estimates

    dataset_subset <- merge(dataset_subset, predict_matrix, all.x = T, by = c("indicator", "age_int", "sex_int", "int_both_md", "int_cross_sectional", "int_online_panel", "cv_cross_panel", "location_id"))
    dataset_subset <- merge(dataset_subset, predict_matrix_lt, all.x = T, by = c("linear_time", "m_mid_age", "m_percent_female", "cv_both_md", "cv_cross_panel", "location_id"))
    dataset_subset[, `:=` (resid = logit_dif - adjustment, resid_lt = logit_dif - adjustment_lt)]
    resid_sd <- sd(dataset_subset[location_id != l,resid])
    
    alpha <- sd(dataset_subset[location_id != l,resid / logit_dif_se])
    dataset_validity <- dataset_subset[location_id == l,]
    dataset_validity[, `:=` (adjustment_low = adjustment - adjustment_sd * alpha * qnorm(0.975, 0, 1), adjustment_high = adjustment + adjustment_sd * alpha * qnorm(0.975, 0, 1))]
    dataset_validity[, adjustment_sd := sqrt((adjustment_sd*alpha)^2 + resid_sd^2)]
    dataset_validity[, `:=` (adjust_stochastic_low = adjustment - adjustment_sd * qnorm(0.975, 0, 1), adjust_stochastic_high = adjustment + adjustment_sd * qnorm(0.975, 0, 1))]
    dataset_validity[, `:=` (logit_dif_low = logit_dif-logit_dif_se*qnorm(0.975, 0, 1), logit_dif_high = logit_dif+logit_dif_se*qnorm(0.975, 0, 1))]
    dataset_validity[, `:=` (sig_dif_stochastic_noui = ifelse((logit_dif > adjust_stochastic_low & logit_dif < adjust_stochastic_high), "No", "Yes"))]
    dataset_validity[, `:=` (resid_sq = sum(resid^2), resid_sq_lt = sum(resid_lt^2))]
    dataset_validity[, `:=` (rmse = sqrt(mean(resid^2)), rmse_lt = sqrt(mean(resid_lt^2)))]


    
    dataset_validity <- dataset_validity[,.(disorder = d, location_name, field_citation_value, cohort, cv_cross_sectional, cv_online_panel, cv_cross_panel, cv_both_md, mid_age, percent_female, obs_change = paste0(round(logit_dif, 2), " (", round(logit_dif_low, 2), " to ", round(logit_dif_high, 2), ")"),
                        pred_change = paste0(round(adjustment, 2), " (", round(adjust_stochastic_low, 2), " to ", round(adjust_stochastic_high, 2), ")"),
                        sig_dif_stochastic_noui,
                        rmse = round(rmse,2), rmse_lt = round(rmse_lt,2))]
    dataset_validity[, rmse_lt_better := ifelse(rmse_lt < rmse, "Yes", "No")]
    

    if(exists("validity_table") == F){
      validity_table <- dataset_validity
    } else{
      validity_table <- rbind(validity_table, dataset_validity)
    }
    print(paste0("Finished ", d, " location ", l))
  }
}

View(validity_table)

giant_locations <- c("United States of America|USA", "United Kingdom|GBR")

# MDD ---------------------------------------------------------------------
### Inspect RMSE
table(unique(validity_table[disorder == "mdd",.(location_name, rmse_lt_better)]))

## RMSE
mean(unique(validity_table[disorder == 'mdd',.(location_name, rmse)])$rmse) # 0.53
## RMSE LT
mean(unique(validity_table[disorder == 'mdd',.(location_name, rmse_lt)])$rmse_lt) # 0.55

## remove giants

# RMSE
mean(unique(validity_table[disorder == 'mdd' & !(location_name %in% giant_locations),.(location_name, rmse)])$rmse) 
mean(unique(validity_table[disorder == 'mdd' & !(location_name %in% giant_locations),.(location_name, rmse_lt)])$rmse_lt)

# Anxiety disorders -------------------------------------------------------

table(unique(validity_table[disorder == "anx",.(location_name, rmse_lt_better)])) # Best

## RMSE
mean(unique(validity_table[disorder == 'anx',.(location_name, rmse)])$rmse) # 0.54
## RMSE LT
mean(unique(validity_table[disorder == 'anx',.(location_name, rmse_lt)])$rmse_lt) # 0.61

## remove giants
# RMSE
mean(unique(validity_table[disorder == 'anx' & !(location_name %in% giant_locations),.(location_name, rmse)])$rmse) # 0.48
mean(unique(validity_table[disorder == 'anx' & !(location_name %in% giant_locations),.(location_name, rmse_lt)])$rmse_lt) # 0.62


# -------------------------------------------------------------------------
validity_table_final <- validity_table[,.(Disorder = ifelse(disorder == "mdd", "MDD", "Anxiety disorders"), Location = location_name, Study = field_citation_value, cv_cross_panel, cv_both_md,
                                    `Mid age` = round(mid_age), `% female` = round(percent_female*100), `Observation (CI)` = obs_change,
                                    `Prediction (PI)` = pred_change, `Outside PI` = sig_dif_stochastic_noui, `Model RMSE` = rmse,
                                    `Benchmark RMSE` = rmse_lt)]
validity_table_final[Location == "United States of America|USA", Location := "USA"]
validity_table_final[Location == "China|CHN", Location := "China"]
validity_table_final[Location == "Japan|JPN", Location := "Japan"]
validity_table_final[Location == "Czechia|CZE", Location := "Czechia"]
validity_table_final[Location == "United Kingdom|GBR", Location := "UK"]
validity_table_final[Location == "Trøndelag|NOR_53432", Location := "Norway"]
validity_table_final[Location == "Australia|AUS", Location := "Australia"]
validity_table_final[Location == "New Zealand|NZL", Location := "New Zealand"]
validity_table_final[Location == "Austria|AUT", Location := "Austria"]
validity_table_final[Location == "Denmark|DNK"  , Location := "Denmark"]
validity_table_final[Location == "Germany|DEU"  , Location := "Germany"]
validity_table_final[Location == "France|FRA", Location := "France"]
validity_table_final[Location == "Ireland|IRL", Location := "Ireland"]
validity_table_final[Location == "Netherlands|NLD", Location := "Netherlands"]
validity_table_final[Location == "Spain|ESP" , Location := "Spain"]

validity_table_final[grepl('Marroquína', Study), Study := "Marroquin et al"]
validity_table_final[grepl('Wanberg', Study), Study := "Wanberg et al"]
validity_table_final[grepl('Katz', Study), Study := "Katz et al"]
validity_table_final[Study == "McGinty, E. E., Presskreischer, R., Han, H., Barry, C. L., (2020). Psychological Distress and Loneliness Reported by US Adults in 2018 and April 2020. JAMA, 324(1):93–94. doi:10.1001/jama.2020.9740", Study := "McGinty et al (a)"]
validity_table_final[grepl('Twenge', Study), Study := "Twenge et al"]
validity_table_final[Study == "McGinty EE, Presskreischer R, Anderson KE, Han H, Barry CL. Psychological Distress and COVID-19–Related Stressors Reported in a Longitudinal Cohort of US Adults in April and July 2020. JAMA. Published online November 23, 2020. doi:10.1001/jama.2020.21231", Study := "McGinty et al (b)"]
validity_table_final[grepl('Bryan', Study), Study := "Bryan et al"]
validity_table_final[grepl('Kantor', Study), Study := "Kantor et al"]
validity_table_final[grepl('Wilson', Study), Study := "Wilson et al"]
validity_table_final[grepl('Ettman', Study), Study := "Ettman et al"]
validity_table_final[grepl('Killgore', Study), Study := "Killgore et al"]
validity_table_final[grepl('Maxfield', Study), Study := "Maxfield et al"]
validity_table_final[Study == "Zhou Y, MacGeorge EL, Myrick JG. Mental Health and Its Predictors during the Early Months of the COVID-19 Pandemic Experience in the United States. International Journal of Environmental Research and Public Health. 2020; 17(17):6315. https://doi.org/10.3390/ijerph17176315", Study := "Zhou et al"]
validity_table_final[grepl('Vieira', Study), Study := "Vieira et al"]
validity_table_final[grepl('Zhang', Study), Study := "Zhang et al"]
validity_table_final[grepl('Choi', Study), Study := "Choi et al"]
validity_table_final[grepl('Kikuchi', Study), Study := "Kikuchi et al"]
validity_table_final[grepl('Yamamoto', Study), Study := "Yamamoto et al"]
validity_table_final[grepl('Kiuchi', Study), Study := "Kiuchi et al"]
validity_table_final[grepl('Ueda', Study), Study := "Ueda et al"]
validity_table_final[grepl('Fukase', Study), Study := "Fukase et al"]
validity_table_final[grepl('Winkler', Study), Study := "Winkler et al"]
validity_table_final[Study == "Shevlin, M., McBride, O., Murphy, J., Gibson Miller, J., Hartman, T. K., Levita, L., … Bentall, R. (2020, April 18). Anxiety, Depression, Traumatic Stress, and COVID-19 Related Anxiety in the UK General Population During the COVID-19 Pandemic. https://doi.org/10.31234/osf.io/hb6nq", Study := "Shevlin et al"]
validity_table_final[grepl("O'Connor", Study), Study := "O'Connor et al"]
validity_table_final[grepl("Groarke", Study), Study := "Groarke et al"]
validity_table_final[Study == "Pieh, Christoph, Budimir, Sanja, Delgadillo, Jaime, Barkham, Michael, Fontaine, Johnny & Probst, Thomas. (2020). Mental health during COVID-19 lockdown in the United Kingdom. Psychosomatic Medicine, Advance on-line publication. Retrieved from http://ovidsp.ovid.com/ovidweb.cgi?T=JS&PAGE=reference&D=ovftw&NEWS=N&AN=00006842-900000000-98497. https://doi.org/10.1097/PSY.0000000000000871 ", Study := "Pieh et al (c)"]
validity_table_final[grepl("Widnall", Study), Study := "Widnall et al"]
validity_table_final[grepl("Vizard", Study), Study := "Vizard et al"]
validity_table_final[grepl("Kwong", Study), Study := "Kwong et al"]
validity_table_final[Study == "Daly, M., Sutin, A., & Robinson, E. (2020, September 5). Longitudinal changes in mental health and the COVID-19 pandemic: Evidence from the UK Household Longitudinal Study. https://doi.org/10.31234/osf.io/qd5z7", Study := "Daly et al (a)"]
validity_table_final[grepl("Knudsen", Study), Study := "Knudsen et al"]
validity_table_final[Study == "Australian Bureau of Statistics (2020). Household Impacts of COVID-19 Survey (November 2020). ", Study := "ABS (2020)"]
validity_table_final[grepl("Biddle", Study), Study := "Biddle et al"]
validity_table_final[grepl("Sibley", Study), Study := "Sibley et al"]
validity_table_final[grepl("Every-Palmer", Study), Study := "Every-Palmer et al"]
validity_table_final[grepl("Bulbulia", Study), Study := "Bulbulia et al"]
validity_table_final[Study == "Pieh, C., Budimir, S., Probst, T. (2020). The effect of age, gender, income, work, and physical activity on mental health during coronavirus disease (COVID-19) lockdown in Austria. Journal of Psychosomatic Research, 136(110186). https://doi.org/10.1016/j.jpsychores.2020.110186", Study := "Pieh et al (b)"]
validity_table_final[Study == "Sønderskov KM, Dinesen PT, Santini ZI, and Østergaard SD. (2020) The depressive state of Denmark during the COVID-19 pandemic. Acta Neuropsychiatrica 32:226–228. doi: 10.1017/neu.2020.15", Study := "Sønderskov et al (b)"]
validity_table_final[Study == "Sønderskov KM, Dinesen PT, Santini ZI, and Østergaard SD. (2020) Increased psychological well-being after the apex of the COVID-19 pandemic. Acta Neuropsychiatrica 32:277–279. doi: 10.1017/neu.2020.26", Study := "Sønderskov et al (a)"]
validity_table_final[grepl("Peretti-Watel", Study), Study := "Peretti-Watel et al"]
validity_table_final[grepl("Peters", Study), Study := "Peters et al"]
validity_table_final[Study == "Daly, M., MacLachlan, M., Maquire, R., Power, J. M., Shevlin, M., Spikol, E., Vallièrese, F., Hyland, P. (in prep).  Changes in PTSD, depression, and generalized anxiety before and during the COVID-19 pandemic in the Republic of Ireland. ", Study := "Daly et al (b)"]
validity_table_final[grepl("Velden", Study), Study := "Van der Velden et al"]
validity_table_final[grepl("Ayuso-Mateos", Study), Study := "Ayuso-Mateos et al"]
validity_table_final[grepl("Valiente", Study), Study := "Valiente et al"]
validity_table_final[grepl("France métropolitaine", Study), Study := "Sante Publique France"]
validity_table_final[grepl("Ravens-Sieberer", Study), Study := "Ravens-Sieberer et al"]

write.xlsx(validity_table_final, "/FILEPATH/validity_table.xlsx")

table(unique(dataset_mdd[,.(location_name, cohort)])$location_name)
table(unique(dataset_anx[,.(location_name, cohort)])$location_name)


validity_table[disorder == 'mdd',.(location_name, mean(rmse - rmse_lt)), by = 'location_name']
validity_table[disorder == 'anx',.(location_name, mean(rmse - rmse_lt)), by = 'location_name']


dim(dataset_mdd[location_id == 102, ])[1]/dim(dataset_mdd)[1]
dim(dataset_mdd[location_id == 95, ])[1]/dim(dataset_mdd)[1]

dim(dataset_anx[location_id == 102, ])[1]/dim(dataset_anx)[1]
dim(dataset_anx[location_id == 95, ])[1]/dim(dataset_anx)[1]

## % estimates by location - MDD
round(100*table(dataset_mdd$location_name)/dim(dataset_mdd)[1], 2)

## % sources by location - MDD
round(100*table(unique(dataset_mdd, by =c('cohort', "location_name"))$location_name)/dim(unique(dataset_mdd, by =c('cohort', "location_name")))[1], 2)

## % estimates by location - anxiety
round(100*table(dataset_anx$location_name)/dim(dataset_anx)[1], 2)

## % sources by location - anxiety
round(100*table(unique(dataset_anx, by =c('cohort', "location_name"))$location_name)/dim(unique(dataset_anx, by =c('cohort', "location_name")))[1], 2)


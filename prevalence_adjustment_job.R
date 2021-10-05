
source("/FILEPATH/get_ids.R")
source("/FILEPATH/get_draws.R")
source("/FILEPATH/interpolate.R")
source("/FILEPATH/get_location_metadata.R")
source("/FILEPATH/get_population.R")
library(mrbrt001, lib.loc = '/FILEPATH/')

rlogit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}

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

args<-commandArgs(trailingOnly = TRUE)
location <- args[1]

epi_age_groups <- c(2:3, 388:389, 238, 34,  6:20, 30:32, 235)

parent <- get_location_metadata(location_set_id = 35, gbd_round_id = 7)[location_id == location, parent_id]
grand_parent <- get_location_metadata(location_set_id = 35, gbd_round_id = 7)[location_id == parent, parent_id]
# Load indicator models ---------------------------------------------------
model_indicator_mdd <- py_load_object(filename = "/FILEPATH/mdd_indicator_model.pkl", pickle = "dill")
model_indicator_anx <- py_load_object(filename = "/FILEPATH/anx_indicator_model.pkl", pickle = "dill")

# Load adjustment models --------------------------------------------------
model_mdd <- py_load_object(filename = "/FILEPATH/adj_model_mdd.pkl", pickle = "dill")
model_anxiety <- py_load_object(filename = "/FILEPATH/adj_model_anx.pkl", pickle = "dill")

# Load indicator data -----------------------------------------------------

  ########## Human mobility #########
  loc_mobility <- fread("/FILEPATH/mobility.csv")

  loc_mobility[, mobility_reference := (-mean)/100] # reverse direction to reflect the 'drop'
  loc_mobility[, `:=` (mean = NULL, upper= NULL, lower = NULL, modeled = NULL, observed = NULL)]

  loc_mobility <- loc_mobility[location_id == location,]

  loc_mobility[, date_n := as.numeric(as.Date(paste0(date)) - as.Date(0, origin="1899-12-30", tz='UTC')), by = "date"]
  if(location %in% loc_mobility$location_id){
    loc_mobility <- loc_mobility[location_id == location & date < 44197, ]
  } else if(parent %in% loc_mobility$location_id){
    loc_mobility <- loc_mobility[location_id == parent & date < 44197, ]
    loc_mobility[, location_id := location]
  } else {
    loc_mobility <- loc_mobility[location_id == grand_parent & date < 44197, ]
    loc_mobility[, location_id := location]
  }

  ########## Covid-19 daily infections #########
  ihme_infections <- fread("/FILEPATH/daily_infections.csv") 
  ihme_infections <- ihme_infections[grepl("2020", date),]
  ihme_infections <- melt.data.table(ihme_infections, id.vars = names(ihme_infections)[!(names(ihme_infections) %like% "draw")], value.name="infections", variable.name="draw")
  ihme_infections[is.na(infections), infections := 0]
  ihme_infections[, infections := mean(infections), by = c("location_id", "date")]
  ihme_infections <- unique(ihme_infections[,.(location_id, date, infections)])
  ihme_infections[, date := as.numeric(as.Date(paste0(date)) - as.Date(0, origin="1899-12-30", tz='UTC')), by = "date"]
  if(location %in% ihme_infections$location_id){
      ihme_infections <- ihme_infections[location_id == location, ]
      infections_location <- location
  } else if(parent %in% ihme_infections$location_id){
      ihme_infections <- ihme_infections[location_id == parent, ]
      ihme_infections[, location_id := location]
      infections_location <- parent
  } else {
      ihme_infections <- ihme_infections[location_id == grand_parent, ]
      ihme_infections[, location_id := location]
      infections_location <- grand_parent
  }
  population_infections <- get_population(age_group_id = 22, location_id = infections_location, year_id = 2020, gbd_round_id = 7, decomp_step = 'iterative')
  population_infections[, location_id := infections_location]
  ihme_infections[, infection_rate := infections / population_infections$population]
  ihme_infections[, covid_infections_sqrt := sqrt(infection_rate)]

  ########## Create indicators #########
  indicator <- merge(loc_mobility, ihme_infections, by = c("location_id", "date"), all = T)
  indicator[is.na(mobility_reference), mobility_reference := 0]

  indicator[, indicator_mdd := mobility_reference * as.numeric(model_indicator_mdd$fe_soln["social_mobility"]) + covid_infections_sqrt * as.numeric(model_indicator_mdd$fe_soln["covid_infections_sqrt"])]
  indicator[, indicator_anx := mobility_reference * as.numeric(model_indicator_anx$fe_soln["social_mobility"]) + covid_infections_sqrt * as.numeric(model_indicator_anx$fe_soln["covid_infections_sqrt"])]

  rm(loc_mobility, ihme_infections)

  mean_mid_age_mdd <- as.numeric(fread("/FILEPATH/mean_mid_age_mdd.csv"))
  mean_mid_age_anx <- as.numeric(fread("/FILEPATH/mean_mid_age_anx.csv"))

for(d in c(1981, 1989)){
  dismod_prev <- interpolate(gbd_id_type = "modelable_entity_id", gbd_id = d, source = "epi", age_group_id = epi_age_groups, measure_id = c(5), location_id = location, reporting_year_start=1990, reporting_year_end=2022, sex_id = c(1,2), status = "best", decomp_step='iterative', gbd_round_id = 7)
  dismod_inc <- interpolate(gbd_id_type = "modelable_entity_id", gbd_id = d, source = "epi", age_group_id = epi_age_groups, measure_id = c(6), location_id = location, reporting_year_start=1990, reporting_year_end=2022, sex_id = c(1,2), status = "best", decomp_step='iterative', gbd_round_id = 7)
  dismod_prev <- rbind(dismod_prev, dismod_inc)
  dismod_prev[, model_version_id := NULL]
  dismod_prev[, modelable_entity_id := NULL]
  if(d == 1981){
    dismod_prev[, modelable_entity_id  := 26756]
  } else {
    dismod_prev[, modelable_entity_id  := 26759]
  }

  for(y in dismod_prev$year_id[dismod_prev$year_id != 2020]){
    write.csv(dismod_prev[year_id == y,], paste0("/FILEPATH/", d, "/prev_", location, "_", y,".csv"), row.names=F)
  }

  dismod_prev <- melt.data.table(dismod_prev[year_id == 2020], id.vars = names(dismod_prev)[!(names(dismod_prev) %like% "draw")], value.name="raw_prev", variable.name="draw")

  dismod_inc <- dismod_prev[measure_id == 6,]
  dismod_prev <- dismod_prev[measure_id == 5,]

  ages <- get_ids('age_group')[age_group_id %in%  epi_age_groups,]
  ages[, `:=` (age_start = as.numeric(unlist(strsplit(age_group_name, " "))[1]), age_end = as.numeric(unlist(strsplit(age_group_name, " "))[3])), by = "age_group_id"]
  ages[age_group_id %in% c(2, 3, 238, 388, 389), `:=` (age_start = 0, age_end = 0)]
  ages[age_start == 95, age_end := 99]
  ages[, mid_age := (age_start+age_end)/2]

  dismod_prev <- merge(dismod_prev, ages[,.(age_group_id, mid_age)], all.x = T, by = "age_group_id")
  dismod_prev[, percent_female := sex_id - 1.5] # mean-centre at 50% female (sex_id, 1 = male, 2 = female)

  if(d == 1981){
    dismod_prev <- merge(dismod_prev, indicator[,.(location_id, date, indicator = indicator_mdd)], by = "location_id", all.x = T, allow.cartesian=TRUE)
    dismod_prev[, m_mid_age := mid_age - mean_mid_age_mdd]
  } else{
    dismod_prev <- merge(dismod_prev, indicator[,.(location_id, date, indicator = indicator_anx)], by = "location_id", all.x = T, allow.cartesian=TRUE)
    dismod_prev[, m_mid_age := mid_age - mean_mid_age_anx]
  }

  dismod_prev[is.na(indicator), indicator := 0]

  dismod_prev[, `:=` (age_int = m_mid_age * indicator, sex_int = percent_female * indicator)]

  # Create draws ------------------------------------------------------------

  if(d == 1981){
    predict_matrix <- unique(dismod_prev[,.(intercept = 0, indicator, age_int, sex_int, cv_cross_panel = 0)])
  } else{
    predict_matrix <- unique(dismod_prev[,.(intercept = 0, indicator, age_int, sex_int, int_both_md = 0, cv_cross_panel = 0)])
  }

  predict_data <- MRData()
  predict_data$load_df(data = predict_matrix, col_covs=as.list(names(predict_matrix)))

  if(d == 1981){
    beta_samples <- mrbrt001::core$other_sampling$sample_simple_lme_beta(1000L, model_mdd)
    
    gamma_outer_samples <- matrix(rep(model_mdd$gamma_soln, each = 1000L), nrow = 1000L)
    gamma_outer_samples[,1] <- 0 # remove intercept gamma
    draws <- model_mdd$create_draws(predict_data,
                                    beta_samples = beta_samples,
                                    gamma_samples = gamma_outer_samples,
                                    random_study = F)
  } else {
    beta_samples <- mrbrt001::core$other_sampling$sample_simple_lme_beta(1000L, model_anxiety)
    
    gamma_outer_samples <- matrix(rep(model_anxiety$gamma_soln, each = 1000L), nrow = 1000L)
    gamma_outer_samples[,1] <- 0 # remove intercept gamma
    draws <- model_anxiety$create_draws(predict_data,
                                        beta_samples = beta_samples,
                                        gamma_samples = gamma_outer_samples,
                                        random_study = F)
  }
  draws <- data.table(draws)
  names(draws) <- paste0("draw_", as.numeric(gsub("V", "", names(draws)))-1)
  predict_matrix <- cbind(predict_matrix, draws)

  predict_matrix <- melt.data.table(predict_matrix, id.vars = names(predict_matrix)[!(names(predict_matrix) %like% "draw")], value.name="adjustment", variable.name="draw")
  predict_matrix[, `:=` (cv_cross_panel = NULL, cv_cross_sectional = NULL, cv_both_md = NULL, cv_online_panel = NULL)]

  dismod_prev <- merge(dismod_prev, predict_matrix, all.x = T, by = c("indicator", "age_int", "sex_int", "draw"))
  
  # Apply adjustments -------------------------------------------------------

  dismod_prev[, `:=` (raw_prev = logit(raw_prev))]
  dismod_prev[, adj_prev := raw_prev - adjustment]
  dismod_prev[, `:=` (adj_prev = rlogit(adj_prev), raw_prev = rlogit(raw_prev))]

  ## Estimate annual point prevalence
  dismod_prev[, `:=` (annual_prev = mean(adj_prev)), by = c("draw", "age_group_id", "sex_id")]

  dismod_prev <- unique(dismod_prev[, .(modelable_entity_id, draw, location_id, age_group_id, measure_id, sex_id, year_id, metric_id, raw_prev, annual_prev)])

  ## Adjust incidence data
  dismod_prev[, ratio := annual_prev / raw_prev]
  dismod_prev[is.na(ratio), ratio := 1]

  dismod_inc <- merge(dismod_inc, dismod_prev[,.(age_group_id, sex_id, draw, ratio)], all = T, by = c("draw", "age_group_id", "sex_id"))
  dismod_inc[, annual_prev := raw_prev * ratio] # annual_prev is just to aid rbind and dcast functions below, but still represents incidence. raw_prev = raw incidence in incidence data frame.

  dismod_prev <- rbind(dismod_prev, dismod_inc)

  dismod_prev[, `:=` (raw_prev = NULL, ratio = NULL)]

  dismod_prev <- dcast(dismod_prev, modelable_entity_id + location_id + age_group_id + measure_id + sex_id + year_id + metric_id ~ draw, value.var="annual_prev")

  write.csv(dismod_prev, paste0("/FILEPATH/", d, "/prev_", location, "_2020.csv"), row.names=F)
}



### PREP THE DATA FROM THE GBD SURVEY IN ORDER TO CROSSWALK FROM SF-12 SCORES TO GBD DISABILITY WEIGHTS

library(data.table)
library(foreign)
library(nlme)
library(ggplot2)
library(mrbrt001, lib.loc = '/FILEPATH/')
rlogit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}

# Prep data for analysis -------------------------------------------------------

# open data file and keep just mcs and pcs info
sf_data <- fread("/FILEPATH/sf12_dw_survey_results.csv")
sf_data <- sf_data[,.(id=Member.LoginName, pcs=Score.PCS, mcs=Score.MCS, composite=Score.PCS+Score.MCS)]

# merge login ids to get lay decriptions associated with them
lay_descriptions <- data.table(read.dta("/FILEPATH/lay_descriptions.dta"))
sf_data <- merge(sf_data, lay_descriptions[,.(id, laydescriptions = conditiondescription)], by = "id", all.x = T)

# merge on disability weights
more_descriptions <- data.table(read.dta("/FILEPATH/desc_id_merge.dta"))
sf_data <- merge(sf_data, more_descriptions, by = "laydescriptions", all = T)

# fill missing hhseqid values
sf_data[grepl("has constant neck pain and arm pain, and difficulty turning the head, holding arms up, and lifting things", laydescriptions), `:=` (hhseqid = 49, group = "Muskuloskeletal")]
sf_data[grepl("has severe neck pain, and difficulty turning the head and lifting things", laydescriptions), `:=` (hhseqid = 48, group = "Muskuloskeletal")]

dw <- fread("/FILEPATH/dw_full.csv")
dw <- melt.data.table(dw, id.vars = names(dw)[!(names(dw) %like% "draw")], value.name = "dw", variable.name = "draw")
dw <- dw[!is.na(dw)]
dw[, `:=` (dw_mean = mean(dw), dw_lower = quantile(dw, 0.025), dw_upper = quantile(dw, 0.975)), by = "hhseqid"]
dw <- unique(dw[,.(hhseqid, healthstate, dw_mean, dw_lower, dw_upper, dw_se = (dw_upper - dw_lower)/3.92)])
sf_data <- merge(sf_data, dw, by = 'hhseqid', all = T)

sf_data[!is.na(composite) & is.na(dw_mean),unique(laydescriptions)]

sf_data <- sf_data[!is.na(composite) & !is.na(dw_mean),]
sf_data[, `:=` (dw_logit = logit(dw_mean), dw_lower_logit = logit(dw_lower), dw_upper_logit = logit(dw_upper))]
sf_data[, `:=` (sf12 = composite, dw_logit_se = (dw_upper_logit - dw_lower_logit)/(qnorm(0.975, 0, 1)*2))]

max_sf12 <- sf_data[, max(sf12)]
sf_data[, sf12_c := sf12 - max_sf12]
sf_data[, sf_outliers := 0]
sf_data[, row_id := seq_len(.N)]

# Functions for MR-BRT model --------------------------------------------------------
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

mrbrt_aic <- function(m){
  log_lik <- -m$get_objective()
  aic_val <- -2 * log_lik + 2 * (length(m$beta_soln) + length(m$gamma_soln))
  return(as.numeric(aic_val))
}

mrbrt_bic <- function(m){
  log_lik <- -m$get_objective()
  used_data <- data.table(cbind(m$data$to_df(), data.frame(w = m$w_soln)))
  bic_val <- -2 * log_lik + (length(m$beta_soln) + length(m$gamma_soln)) * log(used_data[,sum(w)])
  return(as.numeric(bic_val))
}

run_model <- function(trim, int_cov_list, sf12_cov_list, gamma, show_trimmed){
  cov_list <- c(int_cov_list,
                sf12_cov_list)
  model <- MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct = trim)
  model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)

  betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
  betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
  betas <- betas[,.(cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]

  used_data <- data.table(cbind(model$data$to_df(), data.frame(w = model$w_soln)))

  predict_matrix <- data.table(intercept = 1, sf12_c = seq(40, 120)- max_sf12)
  predict_data <- MRData()
  predict_data$load_df(data = predict_matrix, col_covs=as.list(model$cov_names))

  beta_samples <- mrbrt001::core$other_sampling$sample_simple_lme_beta(1000L, model)
  gamma_outer_samples <- matrix(rep(model$gamma_soln, each = 1000L), nrow = 1000L)
  draws <- model$create_draws(predict_data,
                              beta_samples = beta_samples,
                              gamma_samples = gamma_outer_samples,
                              random_study = gamma)

  predict_matrix$pred <- apply(draws, 1, function(x) mean(x))
  predict_matrix$pred_lo <- apply(draws, 1, function(x) quantile(x, 0.025))
  predict_matrix$pred_hi <- apply(draws, 1, function(x) quantile(x, 0.975))
  predict_matrix[, `:=` (sf12 = sf12_c + max_sf12, pred_lo = rlogit(pred_lo), pred_hi = rlogit(pred_hi), pred = rlogit(pred))]

  used_data <- data.table(cbind(model$data$to_df(), data.frame(w = model$w_soln)))
  used_data[, `:=` (dw = rlogit(obs), sf12 = sf12_c + max_sf12)]


  plot <- ggplot(data=predict_matrix, aes(x=sf12, y=pred), fill = "blue")+
    geom_ribbon(data= predict_matrix, aes(x=sf12, ymin=pred_lo, ymax=pred_hi),  fill="blue", alpha=.7) +
    geom_line(size=1) +
    ylab("Disability Weight") +
    xlab("SF-12") +
    theme_minimal() +
    scale_x_continuous(expand=c(0,0), limits = c(40, 121))+
    scale_y_continuous(expand=c(0,0), limits = c(0, 1))+
    theme(axis.line=element_line(colour="black")) +
    geom_point(data=used_data[w ==1,], aes(x=sf12, y=dw) , color="black", shape=16) +
    geom_point(data=used_data[w ==0,], aes(x=sf12, y=dw) , color="red", shape=16)

  if(show_trimmed == T){
    plot <- plot + geom_point(data=sf_data[sf_outliers == 1, .(sf12, dw = dw_mean)], aes(x=sf12, y=dw) , color="red", shape=16)
  }

  return(list(model, plot, betas))
  
}

# Attempt #1 - Person-level data, no outliering ---------------------------
mr_dataset <- MRData()
mr_dataset$load_df(data = sf_data, col_obs = "dw_logit", col_obs_se = "dw_logit_se", col_covs = as.list(c("sf12_c")), col_study_id = "row_id")
int_cov <- list(LinearCovModel('intercept', use_re = T))
sf12_cov <- list(LinearCovModel('sf12_c', use_re = F))

results <- run_model(trim = 1, int_cov_list = int_cov, sf12_cov_list = sf12_cov, gamma = F, show_trimmed = T)
results[2]
results[3]
mrbrt_aic(results[[1]])
mrbrt_bic(results[[1]])

# Attempt #2 - Person-level data, 10% quantile outliering by health state ---------------------------
sf_trimming <- 0.1 # trimming is cut in half on each side so 0.1 trims 5% on each side
z_trim <- qnorm(1 - sf_trimming/2, 0, 1)

sf_data[, sf_outliers := 0]
sf_data[, `:=` (sf_mean = mean(sf12), sf_sd = sd(sf12)), by = "hhseqid"]
sf_data[, `:=` (sf_lower = sf_mean - z_trim*sf_sd, sf_upper = sf_mean + z_trim*sf_sd)]
sf_data[, `:=` (sf_outliers = ifelse(sf12 < sf_lower | sf12 > sf_upper, 1, 0))]
sf_data[is.na(sf_outliers), sf_outliers := 0]

mr_dataset <- MRData()
mr_dataset$load_df(data = sf_data[sf_outliers == 0,], col_obs = "dw_logit", col_obs_se = "dw_logit_se", col_covs = as.list(c("sf12_c")), col_study_id = "row_id")
int_cov <- list(LinearCovModel('intercept', use_re = T))
sf12_cov <- list(LinearCovModel('sf12_c', use_re = F))

results <- run_model(trim = 1, int_cov_list = int_cov, sf12_cov_list = sf12_cov, gamma = F, show_trimmed = T)
results[2]
results[3]
aic_2 <- mrbrt_aic(results[[1]])
bic_2 <- mrbrt_bic(results[[1]])

# Attempt #3 - Person-level data, 10% quantile outliering by health state, spline ---------------------------
mr_dataset <- MRData()
mr_dataset$load_df(data = sf_data[sf_outliers == 0,], col_obs = "dw_logit", col_obs_se = "dw_logit_se", col_covs = as.list(c("sf12_c")), col_study_id = "row_id")
int_cov <- list(LinearCovModel('intercept', use_re = T))
sf12_cov <- LinearCovModel(
  alt_cov = "sf12_c",
  use_spline = T,
  spline_knots = array(seq(0, 1, by = 1/4)),
  spline_degree = 2L,
  spline_knots_type = 'frequency',
  spline_r_linear = T,
  spline_l_linear = T,
  use_re = F
)

results <- run_model(trim = 1, int_cov_list = int_cov, sf12_cov_list = sf12_cov, gamma = T, show_trimmed = T)
results[2]
results[3]

ggsave(results[2][[1]], filename="/FILEPATH/sf12_model_individuals.pdf", width = 8, height = 6)

aic_3 <- mrbrt_aic(results[[1]])
bic_3 <- mrbrt_bic(results[[1]])

# Attempt #4 - Aggregate-level data ---------------------------
sf_data[sf_outliers == 0, `:=` (sf12_c_mean = mean(sf12_c), sf12_c_sd = sd(sf12_c)), by = 'hhseqid']
sf_data_unique <- unique(sf_data[sf_outliers == 0,.(hhseqid, sf12_c = sf12_c_mean, sf12_c_sd, dw_logit = dw_logit, dw_logit_se = dw_logit_se)])

mr_dataset <- MRData()
mr_dataset$load_df(data = sf_data_unique, col_obs = "dw_logit", col_obs_se = "dw_logit_se", col_covs = as.list(c("sf12_c")), col_study_id = "hhseqid")
int_cov <- list(LinearCovModel('intercept', use_re = T))
sf12_cov <- list(LinearCovModel('sf12_c', use_re = F))

results <- run_model(trim = 1, int_cov_list = int_cov, sf12_cov_list = sf12_cov,  gamma = F, show_trimmed = F)
results[2]
results[3]

# Attempt #4 - Aggregate-level data, spline ---------------------------
mr_dataset <- MRData()
mr_dataset$load_df(data = sf_data_unique, col_obs = "dw_logit", col_obs_se = "dw_logit_se", col_covs = as.list(c("sf12_c")), col_study_id = "hhseqid")
int_cov <- list(LinearCovModel('intercept', use_re = T))
sf12_cov <- LinearCovModel(
  alt_cov = "sf12_c",
  use_spline = T,
  spline_knots = array(seq(0, 1, by = 1/4)),
  spline_degree = 2L,
  spline_knots_type = 'frequency',
  spline_r_linear = T,
  spline_l_linear = T,
  use_re = F
)

results <- run_model(trim = 1, int_cov_list = int_cov, sf12_cov_list = sf12_cov, gamma = F, show_trimmed = F)
results[2]
results[3]

ggsave(results[2][[1]], filename="/FILEPATH/sf12_model_aggregate.pdf", width = 8, height = 6)

model_agg <- results[[1]]
model_agg_betas <- results[[3]]

# Attempt #5 - person-level data, priors from aggregate spline model ---------------------------
mr_dataset <- MRData()
mr_dataset$load_df(data = sf_data[sf_outliers == 0], col_obs = "dw_logit", col_obs_se = "dw_logit_se", col_covs = as.list(c("sf12_c")), col_study_id = "group")

int_cov <- list(LinearCovModel('intercept', use_re = T, prior_beta_uniform = rbind(c(model_agg$beta_soln[1]), c(model_agg$beta_soln[1]))))

sf12_cov <- LinearCovModel(
  alt_cov = "sf12_c",
  use_spline = T,
  spline_knots = array(seq(0, 1, by = 1/4)),
  spline_degree = 2L,
  spline_knots_type = 'frequency',
  spline_r_linear = T,
  spline_l_linear = T,
  prior_beta_uniform = rbind(c(model_agg$beta_soln[2:4]), c(model_agg$beta_soln[2:4])), 
  use_re = F
)



results <- run_model(trim = 1, int_cov_list = int_cov, sf12_cov_list = sf12_cov, gamma = T, show_trimmed = T)
results[2]
results[3]

ggsave(results[2][[1]], filename="/FILEPATH/sf12_model_final.pdf", width = 8, height = 6)

aic_5 <- mrbrt_aic(results[[1]])
bic_5 <- mrbrt_bic(results[[1]])

# Output map --------------------------------------------------------------
model <- results[[1]]
py_save_object(object = model, filename = "/FILEPATH/sf12_to_dw.pkl", pickle = "dill")

predict_matrix <- data.table(intercept = 1, sf12_c = seq(40, 120)- max_sf12)
predict_data <- MRData()
predict_data$load_df(data = predict_matrix, col_covs=as.list(model$cov_names))

beta_samples <- mrbrt001::core$other_sampling$sample_simple_lme_beta(1000L, model)
gamma_outer_samples <- matrix(rep(model$gamma_soln, each = 1000L), nrow = 1000L)
draws <- model$create_draws(predict_data,
                            beta_samples = beta_samples,
                            gamma_samples = gamma_outer_samples,
                            random_study = T)

predict_matrix$dw_logit_mean <- apply(draws, 1, function(x) mean(x))
predict_matrix$dw_logit_lo <- apply(draws, 1, function(x) quantile(x, 0.025))
predict_matrix$dw_logit_hi <- apply(draws, 1, function(x) quantile(x, 0.975))
predict_matrix[, `:=` (sf12 = sf12_c + max_sf12)]

###############

draws <- data.table(draws)
names(draws) <- paste0("draw_", as.numeric(gsub("V", "", names(draws)))-1)
predict_draws <- cbind(predict_matrix, draws)

predict_draws <- melt.data.table(predict_draws, id.vars = names(predict_draws)[!(names(predict_draws) %like% "draw")], value.name="dw_logit", variable.name="draw")
predict_draws[, `:=` (sf12_c = NULL, intercept = NULL)]

write.csv(predict_draws, "/FILEPATH/sf12_to_dw_draws.csv", row.names = F)





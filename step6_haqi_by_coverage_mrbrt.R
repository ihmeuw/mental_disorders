source("/FILEPATH/get_covariate_estimates.R")
library(openxlsx)
library(msm)
library(ggplot2)
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

library(mrbrt001, lib.loc = '/FILEPATH/')

# Set up dataset ----------------------------------------------------------

dataset <- data.table(read.xlsx("/FILEPATH/mat_by_haqi.xlsx"))
dataset <- dataset[standard_error > 0,]
dataset[, location_id := loc_id]
dataset[, mid_year := round((year_start + year_end)/2)]

haqi <- get_covariate_estimates(covariate_id = 1099, location_id = c(unique(dataset$location_id), 492, 514), year_id = unique(dataset$mid_year), release_id = 7)

## China survey conducted in Beijing and Shanghai
china_haqi <- haqi[location_id %in% c(492, 514) & year_id == dataset[location_id == 6, unique(mid_year)],]
china_haqi[, `:=` (location_name = "China", location_id = 6, mean_value = mean(mean_value))]

haqi <- rbind(haqi[location_id != 6, .(location_id, mid_year = year_id, haqi = mean_value)], unique(china_haqi[, .(location_id, mid_year = year_id, haqi = mean_value)]))

dataset <- merge(dataset, haqi, by = c('location_id', "mid_year"))

dataset[, `:=` (mean_logit = logit(mean), se_logit = deltamethod(~log(x1/(1-x1)), mean, standard_error^2)), by = c("mean", "standard_error")]

dataset[, row_id := seq_len(.N)]

# Run logit model ---------------------------------------------------------

mr_dataset_logit <- MRData()
mr_dataset_logit$load_df(data = dataset, col_obs = "mean_logit", col_obs_se = "se_logit",
                   col_covs = as.list(c("haqi", "row_id")), col_study_id = "site")

cov_list <- list(LinearCovModel('intercept', use_re = T))
cov_list <- c(cov_list, list(LinearCovModel("haqi", use_re = F)))

model_logit <- MRBRT(data = mr_dataset_logit, cov_models =cov_list, inlier_pct= 1)

model_logit$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)

betas_logit <- data.table(cov = model_logit$cov_names, coef = as.numeric(model_logit$beta_soln), se = get_beta_sd(model_logit))
betas_logit[, `:=` (Relationship = "Logit", lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas_logit

predict_matrix <- data.table(intercept = 1, haqi = seq(0,100, 0.01))

predict_data <- MRData()
predict_data$load_df(data = predict_matrix, col_covs=as.list(c("haqi")))

beta_samples <- mrbrt001::core$other_sampling$sample_simple_lme_beta(1000L, model_logit)
gamma_outer_samples <- matrix(rep(model_logit$gamma_soln, each = 1000L), nrow = 1000L)
draws <- model_logit$create_draws(predict_data,
                            beta_samples = beta_samples,
                            gamma_samples = gamma_outer_samples,
                            random_study =F)

predict_matrix$pred_logit <- model_logit$predict(data = predict_data)
predict_matrix$pred_logit_lo <- apply(draws, 1, function(x) quantile(x, 0.025))
predict_matrix$pred_logit_hi <- apply(draws, 1, function(x) quantile(x, 0.975))
predict_matrix[, `:=` (pred_logit = rlogit(pred_logit), pred_logit_lo = rlogit(pred_logit_lo), pred_logit_hi = rlogit(pred_logit_hi))]

used_data <- data.table(cbind(model$data$to_df(), data.frame(w = model$w_soln)))
used_data <- merge(used_data, dataset[,.(mean, standard_error, row_id)], by = 'row_id')

residuals_data <- MRData()
residuals_data$load_df(data = used_data[, .(intercept = 1, haqi)], col_covs=as.list(c("haqi")))
used_data$pred <- model_logit$predict(data = residuals_data)

used_data[, resid := (rlogit(pred) - rlogit(obs))^2]
logit_rmse <- sqrt(sum(used_data$resid)/dim(used_data)[1])

# Run linear model --------------------------------------------------------

mr_dataset$load_df(data = dataset, col_obs = "mean", col_obs_se = "standard_error",
                   col_covs = as.list(c("haqi", "row_id")), col_study_id = "site")

model <- MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct= 1) 

model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)

betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (Relationship = "Linear", lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas

beta_samples <- mrbrt001::core$other_sampling$sample_simple_lme_beta(1000L, model)
gamma_outer_samples <- matrix(rep(model$gamma_soln, each = 1000L), nrow = 1000L)
draws <- model$create_draws(predict_data,
                                  beta_samples = beta_samples,
                                  gamma_samples = gamma_outer_samples,
                                  random_study =F)

predict_matrix$pred <- model$predict(data = predict_data)
predict_matrix$pred_lo <- apply(draws, 1, function(x) quantile(x, 0.025))
predict_matrix$pred_hi <- apply(draws, 1, function(x) quantile(x, 0.975))

used_data <- data.table(cbind(model$data$to_df(), data.frame(w = model$w_soln)))
used_data <- merge(used_data, dataset[,.(mean, standard_error, row_id)], by = 'row_id')

residuals_data <- MRData()
residuals_data$load_df(data = used_data[, .(intercept = 1, haqi)], col_covs=as.list(c("haqi")))
used_data$pred <- model$predict(data = residuals_data)

used_data[, resid := (pred - obs)^2]
linear_rmse <- sqrt(sum(used_data$resid)/dim(used_data)[1])

# Plot models and data ----------------------------------------------------

plot <- ggplot(data=predict_matrix, aes(x=haqi, y=pred), fill = "blue")+
  geom_ribbon(data= predict_matrix, aes(x=haqi, ymin=pred_lo, ymax=pred_hi),  fill="blue", alpha=.5) +
  geom_line(size=1, color = 'blue') +
  geom_line(data=predict_matrix, aes(x=haqi, y=pred_logit), color="springgreen4", size=1)+
  geom_ribbon(data= predict_matrix, aes(x=haqi, ymin=pred_logit_lo, ymax=pred_logit_hi),  fill="springgreen4", alpha=0.5) +
  geom_hline(yintercept=0, linetype="dashed", color="dark grey", size=1) +
  ylab("Teatment coverage") +
  xlab("HAQI") +
  theme_minimal() +
  scale_x_continuous(expand=c(0,0), breaks = seq(0, 100, 10), limits = c(0, 100))+
  scale_y_continuous(expand=c(0,0), breaks = seq(-0.2, 0.3, 0.1), limits = c(-0.2, 0.3))+
  theme(axis.line=element_line(colour="black")) +

  geom_point(data=used_data[w ==1,], aes(x=haqi, y=mean) , color="black", size=used_data[w ==1, 0.1/standard_error], shape=16) +
 geom_point(data=used_data[w == 0,], aes(x=haqi, y=mean) , color="red", size=used_data[w ==0, 0.1/standard_error], shape=16)
plot

ggsave(plot, filename="/FILEPATH/mat_by_haqi.pdf", width = 7, height = 6)

# finalise table ----------------------------------------------------------

table <- rbind(betas[,.(Relationship, RMSE = linear_rmse, Covariate = cov, Coefficient = round(coef, 3), `95% UI` = paste0(round(lower, 3), " to ", round(upper, 3)), p)],
               betas_logit[,.(Relationship, RMSE = logit_rmse, Covariate = cov, Coefficient = round(coef, 3), `95% UI` = paste0(round(lower, 3), " to ", round(upper, 3)), p)])
table

predict_matrix <- data.table(intercept = 1, haqi = seq(0,100, 0.01))

predict_matrix <- cbind(predict_matrix, draws)
predict_matrix <- melt.data.table(predict_matrix, id.vars = names(predict_matrix)[!(names(predict_matrix) %in% paste0("V", 1:1000))], variable.name = "draw", value.name = "val")
predict_matrix[, draw := as.numeric(gsub("V", "", draw))]
predict_matrix[, draw := paste0("draw_", draw - 1)]
predict_matrix[val < 0, min_haqi := max(haqi), by = 'draw']

predict_matrix <- unique(predict_matrix[!is.na(min_haqi), .(draw, min_haqi)])

hist(predict_matrix$min_haqi)
mean(predict_matrix$min_haqi)
quantile(predict_matrix$min_haqi, 0.025)
quantile(predict_matrix$min_haqi, 0.975)


write.csv(predict_matrix, "/FILEPATH/min_haqi_draws.csv", row.names = F)




or_2_rr <- function(or, s, p){
  rr = 1
  test_or = (rr*(1 - (s / (p*rr+1-p)))) / (1 - (rr*s)/(p*rr+1-p))
  if(!is.na(test_or) & !is.na(or)){
    if(test_or < or){
      while(test_or < or){
        rr = rr+0.0001
        test_or = (rr*(1 - (s / (p*rr+1-p)))) / (1 - (rr*s)/(p*rr+1-p))
      }
    }
    if(test_or > or){
      while(test_or > or){
        rr = rr-0.0001
        test_or = (rr*(1 - (s / (p*rr+1-p)))) / (1 - (rr*s)/(p*rr+1-p))
      }
    }
  } else if(is.na(test_or) | is.na(or)){
    rr = as.numeric(NA)
    warning("NA value given. Function requires OR, s, and p to estimate RRs")
  }
  return(rr)
}

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

get_gamma_sd <- function(model){
  gamma <- model$gamma_soln
  gamma_fisher <- model$lt$get_gamma_fisher(gamma)
  return(sqrt(diag(solve(gamma_fisher))))
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

rlogit <- function(x){exp(x)/(1+exp(x))}

logit <- function(x){log(x/(1-x))}

funnel_plot_mrbrt_old <- function(mod1){

    # assemble data
  dat1 <- data.table(cbind(model$data$to_df(), data.frame(w = model$w_soln)))

  obs_data <- dat1[,.(row_id, val = obs , se= obs_se, study = study_id, included= w)]
  obs_data[, lower:=val-1.96*se]
  obs_data[, upper:=val+1.96*se]
  obs_data <- obs_data[order(val)]
  obs_data[, data:= 1]
  obs_data[included > 0 & included < 1, included:=0.5]

  # create prediction matrix
  intercept_matrix <- dat1[1, (model$cov_names), with = F]
  intercept_matrix[, (model$cov_names[model$cov_names != 'intercept']) := 0]
  predict_data <- MRData()
  predict_data$load_df(data = intercept_matrix, col_covs=as.list(model$cov_names))
  beta_samples <- mr$core$other_sampling$sample_simple_lme_beta(1000L, model)
  gamma_outer_samples <- matrix(rep(model$gamma_soln, each = 1000L), nrow = 1000L)
  draws_int <- model$create_draws(predict_data, beta_samples = beta_samples, gamma_samples = gamma_outer_samples, random_study = F)
  draws_int_gamma <- model$create_draws(predict_data, beta_samples = beta_samples, gamma_samples = gamma_outer_samples, random_study = T)
  intercept_matrix$pred <- model$predict(data = predict_data)
  intercept_matrix$pred_lo <- apply(draws_int, 1, function(x) quantile(x, 0.025))
  intercept_matrix$pred_hi <- apply(draws_int, 1, function(x) quantile(x, 0.975))
  intercept_matrix$pred_lo_g <- apply(draws_int_gamma, 1, function(x) quantile(x, 0.025))
  intercept_matrix$pred_hi_g <- apply(draws_int_gamma, 1, function(x) quantile(x, 0.975))

  pred_matrix <- dat1[, (c("row_id", model$cov_names)), with = F]
  predict_data <- mr$MRData()
  predict_data$load_df(data = pred_matrix, col_covs=as.list(model$cov_names))
  beta_samples <- mr$core$other_sampling$sample_simple_lme_beta(1000L, model)
  gamma_outer_samples <- matrix(rep(model$gamma_soln, each = 1000L), nrow = 1000L)
  draws <- model$create_draws(predict_data, beta_samples = beta_samples, gamma_samples = gamma_outer_samples, random_study = F)
  draws_gamma <- model$create_draws(predict_data, beta_samples = beta_samples, gamma_samples = gamma_outer_samples, random_study = T)
  pred_matrix$pred <- model$predict(data = predict_data)
  pred_matrix$pred_sd <- apply(draws, 1, function(x) sd(x))
  x_walks <- t(apply(draws, 1, function(x) {x - draws_int}))
  pred_matrix$xwalk <- apply(x_walks, 1, function(x) mean(x))
  pred_matrix$xwalk_se <- apply(x_walks, 1, function(x) sd(x))

  obs_data <- merge(obs_data, pred_matrix[,.(row_id, xwalk, xwalk_se)], by = 'row_id')
  obs_data[, `:=` (adj = val - xwalk, adj_se = sqrt(se^2 + xwalk_se^2))]

  # add results
  results <- data.table("val" = intercept_matrix$pred, study = c("Result", "Result w/ gamma"), lower = c(intercept_matrix$pred_lo, intercept_matrix$pred_lo_g), upper = c(intercept_matrix$pred_hi, intercept_matrix$pred_hi_g) )
  results[,data:= 2]
  results[, included := 1]

  obs_data <- rbind(results, obs_data, fill = T)
  obs_data[, row := 1:nrow(obs_data)]

  cov_names <- mod1$cov_names[!mod1$cov_names=="intercept"]
  header <- paste0(as.character(mod1$data),"\ncovariates: ",paste0(cov_names, collapse=", ") ,"\nbeta: ", round(mod1$beta_soln[1], digits = 3),"     gamma: ", mod1$gamma_soln)

  # forest plot of the data
  p <- ggplot(data=obs_data[included == 1,],
              aes(y = se,x = val, size = 1/(se)))+
    geom_point(shape=1, color="black") +
    geom_vline(xintercept =0, linetype=2, color = "red")+
    geom_vline(xintercept =obs_data[study == "Result", val], linetype=1, color = "black")+
    geom_point(data=obs_data[included == 1,], aes(x=adj, y=adj_se) , color="black", shape=16, size = obs_data[included == 1, 1/adj_se]) + ## adj data
    geom_point(data=obs_data[included == 0,], aes(x=val, y=se) , color="red", shape=1, size = obs_data[included == 0, 1/se]) + ## adj data
    geom_point(data=obs_data[included == 0,], aes(x=adj, y=adj_se) , color="red", shape=16, size = obs_data[included == 0, 1/adj_se]) + ## adj data

    # add rectangles for the mr-brt results
    annotate("rect", xmin=obs_data[study == "Result", lower], xmax=obs_data[study == "Result", upper], ymin=0, ymax=obs_data[,max(adj_se, na.rm=T)], alpha=0.2, fill="purple")+
    annotate("rect", xmin=obs_data[study == "Result w/ gamma", lower], xmax=obs_data[study == "Result w/ gamma", upper], ymin=0, ymax=obs_data[,max(adj_se, na.rm=T)], alpha=0.2, fill="blue")+
    annotate("rect", xmin=obs_data[study == "Result", lower], xmax=obs_data[study == "Result", upper], ymin=0, ymax=-Inf, alpha=1, fill="white")+
    annotate("rect", xmin=obs_data[study == "Result", lower], xmax=obs_data[study == "Result", upper], ymin=Inf, ymax=obs_data[,max(adj_se, na.rm=T)], alpha=1, fill="white")+

    # add funnel plot lines mean +- 1.96*se
    geom_segment(aes(x = obs_data[study == "Result", val], y = 0, xend = obs_data[study == "Result", val]+1.96*obs_data[,max(adj_se, na.rm=T)], yend = obs_data[,max(adj_se, na.rm=T)], size = .01))+
    geom_segment(aes(x = obs_data[study == "Result", val], y = 0, xend = obs_data[study == "Result", val]-1.96*obs_data[,max(adj_se, na.rm=T)], yend = obs_data[,max(adj_se, na.rm=T)], size = .01))+

    ylab('standard error')+ xlab(paste0("test"))+
    labs(subtitle = header)+
      scale_size_continuous(guide = F)+
    scale_y_continuous(trans = "reverse")+theme_bw()+
    #scale_x_continuous(expand=c(0,0))+
    labs(title= "Funnel plot")

  print(p)

}

egger_mr_brt_pval <- function(m){
  detach(package:metafor,unload=TRUE)
  require(nlme)
  # Can include precision as a moderator in multi-level analysis to test for publication bias
  # https://stats.stackexchange.com/questions/155693/metafor-package-bias-and-sensitivity-diagnostics
  obs_data <- data.table(cbind(m$data$to_df(), data.frame(w = m$w_soln)))
  obs_data <- obs_data[w == 1,]

  pred_matrix <- obs_data[, (c("row_id", m$cov_names)), with = F]
  predict_data <- mr$MRData()
  predict_data$load_df(data = pred_matrix, col_covs=as.list(m$cov_names))
  pred_matrix$pred <- m$predict(data = predict_data)

  obs_data <- merge(obs_data, pred_matrix[,.(row_id, pred)], by = 'row_id')
  obs_data[, `:=` (resid = obs - pred)]
  obs_data[, `:=` (w_resid = resid/obs_se, precision = 1/obs_se)]
  require(metafor)

  egg_model <- lme(fixed=w_resid ~ precision, random = ~1|study_id, data=obs_data)
  egg_results <- data.table(summary(egg_model)$tTable)
  p <- egg_results[1,`p-value`]
  return(p)
}

funnel_plot_mrbrt_bydata <- function(m){
  

  obs_data <- data.table(cbind(m$data$to_df(), data.frame(w = m$w_soln)))

  pred_matrix <- obs_data[, (c("row_id", m$cov_names)), with = F]
  predict_data <- mr$MRData()
  predict_data$load_df(data = pred_matrix, col_covs=as.list(m$cov_names))
  pred_matrix$pred <- m$predict(data = predict_data)
  obs_data <- merge(obs_data, pred_matrix[,.(row_id, pred)], by = 'row_id')
  obs_data[, `:=` (resid = obs - pred)]

  residuals <- obs_data$resid
  vei <- obs_data$obs_se
  included <- obs_data$w
  funnelplotdata <- funnel(residuals, vei, xlab = "Residual value")
  my_colors <- c('red','black')[(obs_data$w == 1) + 1]
  with(funnelplotdata, points(x, y, col = my_colors, pch = 19))
}

funnel_plot_mrbrt_bystudy <- function(m){
  

  obs_data <- data.table(cbind(m$data$to_df(), data.frame(w = m$w_soln)))

  pred_matrix <- obs_data[, (c("row_id", m$cov_names)), with = F]
  predict_data <- mr$MRData()
  predict_data$load_df(data = pred_matrix, col_covs=as.list(m$cov_names))
  pred_matrix$pred <- m$predict(data = predict_data)
  obs_data <- merge(obs_data, pred_matrix[,.(row_id, pred)], by = 'row_id')
  obs_data[, `:=` (resid = obs - pred)]

  study_to_pool <- data.table(table(obs_data[w == 1, study_id]))[N > 1, V1]

  study_level_residuals <- obs_data[!(study_id %in% study_to_pool), .(study_id, resid, obs_se, w)]

  for(s in study_to_pool){
    mr_dataset_pool <- mr$MRData()
    mr_dataset_pool$load_df(
      data = obs_data[study_id == s & w == 1,],
      col_obs = "resid", col_obs_se = "obs_se",
      col_covs = as.list(c("row_id")), col_study_id = "study_id" )
    cov_list <- list(mr$LinearCovModel('intercept', use_re = F))
    pool_model <- mr$MRBRT(data = mr_dataset_pool, cov_models =cov_list, inlier_pct =1)
    pool_model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)
    betas <- data.table(study_id = s, resid = as.numeric(pool_model$beta_soln), obs_se = get_beta_sd(pool_model), w = 1)
    study_level_residuals <- rbind(study_level_residuals, betas)
  }

  residuals <- study_level_residuals$resid
  vei <- study_level_residuals$obs_se

  funnelplotdata <- funnel(residuals, vei, xlab = "Residual value")
  my_colors <- c('red','black')[(study_level_residuals$w == 1) + 1]
}



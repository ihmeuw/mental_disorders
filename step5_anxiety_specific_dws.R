
library(data.table)
library(readstata13)
library(openxlsx)
library(mrbrt001, lib.loc = '/FILEPATH/')
rlogit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}

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

# Load in microdata -------------------------------------------------------
microdata <- data.table(read.dta13("/FILEPATH/slim MHS Aus97 12 mnth diagnoses added.dta"))

rename_matrix <- data.table(old = c("c1a", paste0("c", 2:12), "danx12", "ddepa12", "ddepb12", "ddepc12", "ddys12", "dalcd12",
                                    "ddrgd12", "danx1", "ddepa1", "ddepb1", "ddepc1", "ddys1", "dalcd1", "ddrgd1"),
                            new = c("Iasthma", "Ibronchitis", "Ianaemia", "Iblood_press", "Iheart_troub", "Iarthritis",
                                    "Ikidney_dis", "Idiabetes", "Icancer", "Iulcer", "Iliver_gallbladder", "Ihernia_rupture",
                                    "Ianxiety12", "Xmild_depr12", "Xmod_depr12", "Xsevere_depr12", "Idysthymia12",
                                    "Ialcohol_depend12", "Idrug_depend12", "imental_anxiety", "Xmild_depr1", "Xmod_depr1",
                                    "Xsevere_depr1", "imental_unipolar_dys", "imental_alcohol", "Idrug_depend1"))
for(r in rename_matrix$old){
  setnames(microdata, r, rename_matrix[old == r, new])
  microdata[get(rename_matrix[old == r, new]) %in% c(0, 1), paste(rename_matrix[old == r, new]) := 0]
  microdata[get(rename_matrix[old == r, new]) %in% c(2, 5), paste(rename_matrix[old == r, new]) := 1]
}

microdata[, `:=` (Idepression12 = ifelse(Xmild_depr12 == 1 | Xmod_depr12 == 1 | Xsevere_depr12 == 1, 1, 0),
                  imental_unipolar_mdd = ifelse(Xmild_depr1 == 1 | Xmod_depr1 == 1 | Xsevere_depr1 == 1, 1, 0), I_NONE = 1)]
microdata[, agegr := (agegr+2)*5]
microdata[, agegr := paste(agegr, "to", agegr+4)]
microdata[agegr == "15 to 19", agegr := "18 to 19"]
microdata[agegr == "75 to 79", agegr := "75+"]
setnames(microdata, "a1", "sex")
microdata[, id := seq_len(.N)]

## Calculate SF-12 total score ##
microdata[, `:=` (sf12= pcs_12 + mcs_12)]

# Load in SF-12 to DW draws -------------------------------------------------
dw_map_model <- py_load_object(filename = "/FILEPATH/sf12_to_dw.pkl", pickle = "dill")

dw_map <- data.table(intercept = 1, sf12_c = unique(microdata$sf - 120))
predict_data <- MRData()
predict_data$load_df(data = dw_map, col_covs=as.list(dw_map_model$cov_names))

beta_samples <- mrbrt001::core$other_sampling$sample_simple_lme_beta(1000L, dw_map_model)
gamma_outer_samples <- matrix(rep(dw_map_model$gamma_soln, each = 1000L), nrow = 1000L)
draws <- dw_map_model$create_draws(predict_data, beta_samples = beta_samples, gamma_samples = gamma_outer_samples, random_study = T)

dw_map$dw_logit_mean <- apply(draws, 1, function(x) mean(x))
dw_map$dw_logit_lo <- apply(draws, 1, function(x) quantile(x, 0.025))
dw_map$dw_logit_hi <- apply(draws, 1, function(x) quantile(x, 0.975))
dw_map[, `:=` (sf12 = sf12_c + 120, dw_logit_se = (dw_logit_hi - dw_logit_lo)/3.92)]

microdata <- merge(microdata, dw_map[,.(sf12, dw_logit_mean, dw_logit_se)], by = 'sf12', all.x = T)

draws <- data.table(draws)
names(draws) <- paste0("draw_", as.numeric(gsub("V", "", names(draws)))-1)
dw_map <- cbind(dw_map, draws)
dw_map <- melt.data.table(dw_map, id.vars = names(dw_map)[!(names(dw_map) %like% "draw")], value.name="dw_logit", variable.name="draw")

# Prep comorbidity variables for analysis ---------------------------------
setnames(microdata,names(microdata),tolower(names(microdata)))
comos<-grep('^i',names(microdata),value=T)
comos <- comos[!(grepl("12", comos))]

# create a list of comos to predict for
comos_to_predict<-fread("/FILEPATH/ICD_healthstate_GBD2020.csv")
comos_to_predict<-comos_to_predict[!yld_cause=="",]

# create map from source to MEPS dummy variables
source<-unique(comos_to_predict[,source])
dummy_map<-data.table(source)
dummy_map[grep('AHS',source),"AHS":=1]
dummy_map[is.na(AHS),AHS:=0]

# merge dummy map onto comos_to_predict, only keep those that use MEPS as sources
comos_to_predict<-merge(comos_to_predict,dummy_map,by="source")
comos_to_predict<-comos_to_predict[AHS==1]
comos_to_predict[,comos:=paste0("i",yld_cause)]
comos_to_predict<-tolower(unique(comos_to_predict$comos[comos_to_predict$comos %in% names(microdata)]))

# convert condition indicators to numeric variables
microdata[, (comos):=lapply(.SD, as.numeric), .SD=comos]

# drop rows with NAs for conditions
microdata<-microdata[complete.cases(microdata[,comos,with=F]),]

## Retain meaningful comorbidities and drop uninformative ones
mr_dataset <- MRData()
mr_dataset$load_df(
  data = microdata, col_obs = "dw_logit_mean", col_obs_se = "dw_logit_se",
  col_covs = as.list(comos), col_study_id = "id" )
covfinder <- CovFinder(
  data = mr_dataset,
  covs = as.list(comos[!(comos %in% comos_to_predict)]),
  pre_selected_covs = as.list(c("intercept", comos_to_predict)),
  normalized_covs = T,
  #num_samples = 1000L,
  power_range = list(-4, 4),
  power_step_size = 1,
  laplace_threshold = 1e-5,
  inlier_pct = 1
)
covfinder$select_covs(verbose = TRUE)
comos <- covfinder$selected_covs

# Run MR-BRT regression ---------------------------------------------------
comos <- comos[!(comos %in% 'intercept')]
cov_list <- list(LinearCovModel('intercept', use_re = T))
for(c in comos){cov_list <- c(cov_list, list(LinearCovModel(c, use_re = F)))}

model <- MRBRT(data = mr_dataset, cov_models =cov_list, inlier_pct =1)

model$fit_model(inner_print_level = 5L, inner_max_iter = 1000L, inner_acceptable_tol=1e-3)

betas <- data.table(cov = model$cov_names, coef = as.numeric(model$beta_soln), se = get_beta_sd(model))
betas[, `:=` (lower = coef-(qnorm(0.975)*se), upper = coef+(qnorm(0.975)*se), z = abs(coef/se), p = (1 - pnorm(abs(coef/se)))*2)]
betas <- betas[,.(cov, beta = round(coef,3), lower = round(lower,3), upper = round(upper,3), z = round(z, 3), p = round(p,4))]

## Subset to only those with condition of interest
microdata <- microdata[imental_anxiety == 1, ]

predict_matrix <- microdata[, (c(comos, "id")), with = F]
predict_matrix[, `:=`(intercept = 1, imental_anxiety = 0)]
predict_data <- MRData()
predict_data$load_df(data = predict_matrix, col_covs=as.list(model$cov_names))

beta_samples <- mrbrt001::core$other_sampling$sample_simple_lme_beta(1000L, model)
gamma_outer_samples <- matrix(rep(model$gamma_soln, each = 1000L), nrow = 1000L)
draws <- model$create_draws(predict_data,
                            beta_samples = beta_samples,
                            gamma_samples = gamma_outer_samples,
                            random_study = F)
draws <- data.table(draws)
names(draws) <- paste0("draw_", as.numeric(gsub("V", "", names(draws)))-1)
predict_matrix <- cbind(predict_matrix, draws)

predict_matrix <- melt.data.table(predict_matrix, id.vars = names(predict_matrix)[!(names(predict_matrix) %like% "draw")], value.name="dw_cf_logit", variable.name="draw")

## Bring in observed DW draws
microdata <- merge(microdata, dw_map[,.(sf12, dw_logit, draw)], by = "sf12", all.x = T, allow.cartesian = T)

## Merge in counterfactual DW draws
microdata <- merge(microdata, predict_matrix[,.(id, draw, dw_cf_logit)], by = c("id", "draw"), all.x = T)

## Calculate disorder-specific DW
microdata[, `:=` (dw = rlogit(dw_logit), dw_cf = rlogit(dw_cf_logit))]
microdata[, dw_ds := 1 - (1-dw)/(1-dw_cf)]

write.csv(microdata, "/FILEPATH/anxiety_dw_distributions.csv", row.names = F)



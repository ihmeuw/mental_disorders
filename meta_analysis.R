## Clear memory
rm(list=ls()) 
citation()

## Load required packages
if(require(openxlsx) == F){install.packages('openxlsx')} 
if(require(data.table) == F){install.packages('data.table')}
if(require(metafor) == F){install.packages('metafor')}
if(require(msm) == F){install.packages('msm')}

library(metafor)
library(data.table)
library(openxlsx)
library(msm)

logit <- function(x){log(x/(1-x))}
rlogit <- function(x){exp(x)/(1+exp(x))}


# Read in dataset ---------------------------------------------------------
dataset <- as.data.table(read.xlsx("FILEPATH/dataset.xlsx", sheet = "extraction"))

dataset <- dataset[(group_review == 1 | overall == 1) & cv_symptom_scale == 0 & measure == "prevalence", ]
dataset <- dataset[!(recall_type %in% c("Lifetime", "n.s."))]

#dataset[mean != 0 & cases == 0, cases := NA]
dataset[, mean := as.numeric(mean)]
dataset[, age_start := as.numeric(age_start)]
dataset[, age_end := as.numeric(age_end)]
dataset[, age_mean := as.numeric(age_mean)]

dataset[!is.na(sample_size), cases := mean * sample_size]
dataset[!is.na(cases) & is.na(sample_size), sample_size := cases / mean]

dataset[standard_error == 0, standard_error := NA]
dataset[is.na(standard_error) & !is.na(lower), standard_error := (upper - lower) / (qnorm(0.975,0, 1)*2)]
dataset[is.na(standard_error), standard_error := sqrt(1/sample_size * mean * (1-mean) + 1/(4*sample_size^2)*qnorm(0.975,0,1)^2)]

## Set up age
dataset[!is.na(age_mean), age := age_mean]
dataset[is.na(age_mean), age := (age_start + age_end)/2]

## Set up recall
dataset[, `:=` (cv_3to6months = 0, cv_1year = 0)]
dataset[recall_type == "Period: months" & (recall_type_value %in% c(3, 6)), cv_3to6months := 1]
dataset[recall_type == "Period: months" & recall_type_value == 12, cv_1year := 1]
dataset[recall_type == "Period: years" & recall_type_value == 1, cv_1year := 1]

## Set up location_variables
table(dataset$GBD.Regions)
table(dataset$GBD.Super.Regions)
dataset[, GBD.Super.Regions := trimws(GBD.Super.Regions)] # Remove spaces in excel making R think there two separate High-Income regions
dataset[, GBD.Regions := trimws(GBD.Regions)] # As a just-in-case
## We will treat high-income as the reference for the model given this is where most of the data exists.
dataset[, cv_not_hi_location := ifelse(GBD.Super.Regions == "High-income", 0, 1)]

# Prep model variables ----------------------------------------------------
dataset[, yi := logit(mean)]
dataset[!is.na(lower) & lower > 0, vi := ((logit(upper) - logit(lower))/(qnorm(0.975,0, 1)*2))^2]
dataset[, row_id := seq_len(.N)]
dataset[is.na(vi), vi := deltamethod(~log(x1/(1-x1)), mean, standard_error^2)^2, by = "row_id"]

dataset[group != 1, unique_study := as.character(group)]
dataset[group == 1, unique_study := field_citation_value] 

dataset[, obs_number := seq_len(.N)]

dataset[, m_percent_female := as.numeric(percent_female) - 0.5]
length(dataset[is.na(m_percent_female), unique(unique_study)])
dataset[is.na(m_percent_female), m_percent_female := 0]

mean_age_value <- mean(dataset[group_review == 1, age])
mean_age_value

dataset[, m_age := age - mean_age_value]

dataset <- dataset[mean != 0,]
dataset <- dataset[disoder != "Symptoms",]




dataset[, study_author := unlist(strsplit(field_citation_value, ","))[1], by = 'field_citation_value']

dataset[, study_year := unlist(strsplit(field_citation_value, "\\("))[2], by = 'field_citation_value']
dataset[, study_year := substr(study_year, 1, 4)]

dataset[!grepl("Report", field_citation_value), study_label := paste0(study_author, " et al. (", study_year, ")")]
dataset[grepl("Report", field_citation_value), study_label := "NMHS Nepal (2020)"]



dataset_overall <- dataset[overall == 1,]
dataset <- dataset[group_review == 1,]

# Statistical analyses ----------------------------------------------------


# Any SD ------------------------------------------------------------------
data_any_sd <- dataset[disoder == "Any SD", ]


table(unique(data_any_sd[,.(field_citation_value, GBD.Super.Regions)])$GBD.Super.Regions)
table(unique(data_any_sd[,.(field_citation_value, GBD.Super.Regions)])$GBD.Regions)

table(data_any_sd[,.(cv_1year, cv_not_hi_location)])
plot(data_any_sd[,.(m_percent_female, yi)])
data_any_sd[, mean(mean), by = c('cv_1year', 'cv_not_hi_location')]


## Run multi-level regression model
model_anysd_uni_intercept <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_any_sd, slab=study_label)
model_anysd_uni_intercept

model_anysd_uni_sex <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_any_sd, slab=study_label, mods = ~ m_percent_female)
model_anysd_uni_sex 

(model_anysd_uni_intercept$sigma2[1] - model_anysd_uni_sex$sigma2[1])/model_anysd_uni_intercept$sigma2[1]
(model_anysd_uni_intercept$sigma2[2] - model_anysd_uni_sex$sigma2[2])/model_anysd_uni_intercept$sigma2[2]

model_anysd_uni_age <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_any_sd, slab=study_label, mods = ~ m_age)
model_anysd_uni_age 

((model_anysd_uni_intercept$sigma2[1] - model_anysd_uni_age$sigma2[1])/model_anysd_uni_intercept$sigma2[1])*100
((model_anysd_uni_intercept$sigma2[2] - model_anysd_uni_age$sigma2[2])/model_anysd_uni_intercept$sigma2[2])*100

model_anysd_uni_3m_6m_recall <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_any_sd, slab=study_label, mods = ~ cv_3to6months)
model_anysd_uni_3m_6m_recall 

(model_anysd_uni_intercept$sigma2[1] - model_anysd_uni_3m_6m_recall$sigma2[1])/model_anysd_uni_intercept$sigma2[1]
(model_anysd_uni_intercept$sigma2[2] - model_anysd_uni_3m_6m_recall$sigma2[2])/model_anysd_uni_intercept$sigma2[2]

model_anysd_uni_1y_recall <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_any_sd, slab=study_label, mods = ~ cv_1year)
model_anysd_uni_1y_recall 

((model_anysd_uni_intercept$sigma2[1] - model_anysd_uni_1y_recall$sigma2[1])/model_anysd_uni_intercept$sigma2[1])*100
((model_anysd_uni_intercept$sigma2[2] - model_anysd_uni_1y_recall$sigma2[2])/model_anysd_uni_intercept$sigma2[2])*100

model_anysd_uni_no_hi <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_any_sd, slab=study_label, mods = ~ cv_not_hi_location)
model_anysd_uni_no_hi 

((model_anysd_uni_intercept$sigma2[1] - model_anysd_uni_no_hi$sigma2[1])/model_anysd_uni_intercept$sigma2[1])*100
((model_anysd_uni_intercept$sigma2[2] - model_anysd_uni_no_hi$sigma2[2])/model_anysd_uni_intercept$sigma2[2])*100

model_anysd_uni_maybe_bdd <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_any_sd, slab=study_label, mods = ~ maybe_bdd)
model_anysd_uni_maybe_bdd 

((model_anysd_uni_intercept$sigma2[1] - model_anysd_uni_maybe_bdd$sigma2[1])/model_anysd_uni_intercept$sigma2[1])*100
((model_anysd_uni_intercept$sigma2[2] - model_anysd_uni_maybe_bdd$sigma2[2])/model_anysd_uni_intercept$sigma2[2])*100

model_anysd_uni_jbi <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_any_sd, slab=study_label, mods = ~ JBI_Score)
model_anysd_uni_jbi 

((model_anysd_uni_intercept$sigma2[1] - model_anysd_uni_jbi$sigma2[1])/model_anysd_uni_intercept$sigma2[1])*100
((model_anysd_uni_intercept$sigma2[2] - model_anysd_uni_jbi$sigma2[2])/model_anysd_uni_intercept$sigma2[2])*100

model_anysd_multi <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_any_sd, slab=study_label, mods = ~ m_percent_female)
model_anysd_multi

((model_anysd_uni_intercept$sigma2[1] - model_anysd_multi$sigma2[1])/model_anysd_uni_intercept$sigma2[1])*100
((model_anysd_uni_intercept$sigma2[2] - model_anysd_multi$sigma2[2])/model_anysd_uni_intercept$sigma2[2])*100

forest_rma_anysd <- rma(yi = yi, vi = vi, data = dataset_overall[disoder == "Any SD", ] , slab=study_label)
forest(forest_rma_anysd)


pdf(file='FILEPATH/funnelplot_anysd.pdf', width = 10, height = 9)
funnel(forest_rma_anysd)
dev.off() # Turn the PDF device off

pdf(file='FILEPATH/trimandfill_anysd.pdf', width = 10, height = 9)
funnel(trimfill(forest_rma_anysd))
dev.off() # Turn the PDF device off


pub_bias <- function(x, d){
  eggers <-data.table(disorder = d, z = regtest(x)$zval, p = regtest(x)$pval)
  prev_uni <- as.data.table(predict(x, digits=3, transf=transf.ilogit))
  prev_trim <- as.data.table(predict(trimfill(x), digits=3, transf=transf.ilogit))
  eggers[, `:=` (prevalence = paste0(round(prev_uni$pred*100, 1), "% (", round(prev_uni$ci.lb*100, 1), "-", round(prev_uni$ci.ub*100, 1), ")"),
                 prevalence_trimnfill = paste0(round(prev_trim$pred*100, 1), "% (", round(prev_trim$ci.lb*100, 1), "-", round(prev_trim$ci.ub*100, 1), ")"))]
  return(eggers)
}

eggers_table <- pub_bias(forest_rma_anysd, "Any SD")

results_bysex_total <- as.data.table(predict(model_anysd_uni_sex, newmods = rbind(c(0),c(-0.5), c(0.5)), digits=3, transf=transf.ilogit))
results_bysex_total

# SSD ---------------------------------------------------------------------
data_somatization <- dataset[disoder == "SSD", ]
table(unique(data_somatization[,.(field_citation_value, GBD.Super.Regions)])$GBD.Super.Regions)
table(unique(data_somatization[,.(field_citation_value, GBD.Super.Regions)])$GBD.Regions)
unique(data_somatization[,.(field_citation_value, JBI_Score, mean)])

model_sd_uni_intercept <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_somatization, slab=study_label)
model_sd_uni_intercept

model_sd_uni_sex <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_somatization, slab=study_label, mods = ~ m_percent_female)
model_sd_uni_sex 

((model_sd_uni_intercept$sigma2[1] - model_sd_uni_sex$sigma2[1])/model_sd_uni_intercept$sigma2[1])*100
((model_sd_uni_intercept$sigma2[2] - model_sd_uni_sex$sigma2[2])/model_sd_uni_intercept$sigma2[2])*100

model_sd_uni_age <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_somatization, slab=study_label, mods = ~ m_age)
model_sd_uni_age 

((model_sd_uni_intercept$sigma2[1] - model_sd_uni_age$sigma2[1])/model_sd_uni_intercept$sigma2[1])*100
((model_sd_uni_intercept$sigma2[2] - model_sd_uni_age$sigma2[2])/model_sd_uni_intercept$sigma2[2])*100

model_sd_uni_3m_6m_recall <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_somatization, slab=study_label, mods = ~ cv_3to6months)
model_sd_uni_3m_6m_recall 

((model_sd_uni_intercept$sigma2[1] - model_sd_uni_3m_6m_recall$sigma2[1])/model_sd_uni_intercept$sigma2[1])*100
((model_sd_uni_intercept$sigma2[2] - model_sd_uni_3m_6m_recall$sigma2[2])/model_sd_uni_intercept$sigma2[2])*100

model_sd_uni_1y_recall <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_somatization, slab=study_label, mods = ~ cv_1year)
model_sd_uni_1y_recall 

((model_sd_uni_intercept$sigma2[1] - model_sd_uni_1y_recall$sigma2[1])/model_sd_uni_intercept$sigma2[1])*100
((model_sd_uni_intercept$sigma2[2] - model_sd_uni_1y_recall$sigma2[2])/model_sd_uni_intercept$sigma2[2])*100

model_sd_uni_no_hi <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_somatization, slab=study_label, mods = ~ cv_not_hi_location)
model_sd_uni_no_hi 

((model_sd_uni_intercept$sigma2[1] - model_sd_uni_no_hi$sigma2[1])/model_sd_uni_intercept$sigma2[1])*100
((model_sd_uni_intercept$sigma2[2] - model_sd_uni_no_hi$sigma2[2])/model_sd_uni_intercept$sigma2[2])*100

model_sd_uni_jbi <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_somatization, slab=study_label, mods = ~ JBI_Score)
model_sd_uni_jbi 

((model_sd_uni_intercept$sigma2[1] - model_sd_uni_jbi$sigma2[1])/model_sd_uni_intercept$sigma2[1])*100
((model_sd_uni_intercept$sigma2[2] - model_sd_uni_jbi$sigma2[2])/model_sd_uni_intercept$sigma2[2])*100

model_anysd_multi <- model_sd_uni_intercept
model_anysd_multi

forest_rma_sd <- rma(yi = yi, vi = vi, data = dataset_overall[disoder == "SSD", ] , slab=study_label)

forest(forest_rma_sd)

pdf(file='FILEPATH/funnelplot_somatization.pdf', width = 10, height = 9)
funnel(forest_rma_sd)
dev.off() # Turn the PDF device off

pdf(file='FILEPATH/trimandfill_somatization.pdf', width = 10, height = 9)
funnel(trimfill(forest_rma_sd))
dev.off() # Turn the PDF device off

eggers_table <- rbind(eggers_table, pub_bias(forest_rma_sd, "SD"))

eggers_table

results_total <- as.data.table(predict(model_anysd_multi, digits=3, transf=transf.ilogit))
results_total


# IAD ---------------------------------------------------------------------
data_hyp <- dataset[disoder == "IAD", ]

plot(data_hyp[,.(m_age, yi)])

model_iad_uni_intercept <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_hyp, slab=study_label)
model_iad_uni_intercept

model_iad_uni_sex <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_hyp, slab=study_label, mods = ~ m_percent_female)
model_iad_uni_sex 

((model_iad_uni_intercept$sigma2[1] - model_iad_uni_sex$sigma2[1])/model_iad_uni_intercept$sigma2[1])*100
((model_iad_uni_intercept$sigma2[2] - model_iad_uni_sex$sigma2[2])/model_iad_uni_intercept$sigma2[2])*100

model_iad_uni_age <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_hyp, slab=study_label, mods = ~ m_age)
model_iad_uni_age 

((model_iad_uni_intercept$sigma2[1] - model_iad_uni_age$sigma2[1])/model_iad_uni_intercept$sigma2[1])*100
((model_iad_uni_intercept$sigma2[2] - model_iad_uni_age$sigma2[2])/model_iad_uni_intercept$sigma2[2])*100

model_iad_uni_1y_recall <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_hyp, slab=study_label, mods = ~ cv_1year)
model_iad_uni_1y_recall 

((model_iad_uni_intercept$sigma2[1] - model_iad_uni_1y_recall$sigma2[1])/model_iad_uni_intercept$sigma2[1])*100
((model_iad_uni_intercept$sigma2[2] - model_iad_uni_1y_recall$sigma2[2])/model_iad_uni_intercept$sigma2[2])*100

model_iad_uni_no_hi <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_hyp, slab=study_label, mods = ~ cv_not_hi_location)
model_iad_uni_no_hi 

((model_iad_uni_intercept$sigma2[1] - model_iad_uni_no_hi$sigma2[1])/model_iad_uni_intercept$sigma2[1])*100
((model_iad_uni_intercept$sigma2[2] - model_iad_uni_no_hi$sigma2[2])/model_iad_uni_intercept$sigma2[2])*100

model_iad_uni_jbi <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_hyp, slab=study_label, mods = ~ JBI_Score)
model_iad_uni_jbi 

((model_iad_uni_intercept$sigma2[1] - model_iad_uni_jbi$sigma2[1])/model_iad_uni_intercept$sigma2[1])*100
((model_iad_uni_intercept$sigma2[2] - model_iad_uni_jbi$sigma2[2])/model_iad_uni_intercept$sigma2[2])*100

model_iad_multi <- model_iad_uni_intercept
model_iad_multi

forest_rma_hyp <- rma(yi = yi, vi = vi, data = dataset_overall[disoder == "IAD", ] , slab=study_label)
forest(forest_rma_hyp)
unique(data_hyp$study_label)



pdf(file='FILEPATH/funnelplot_hyp.pdf', width = 10, height = 9)
funnel(forest_rma_hyp)
dev.off() # Turn the PDF device off

pdf(file='FILEPATH/trimandfill_hypn.pdf', width = 10, height = 9)
funnel(trimfill(forest_rma_hyp))
dev.off() # Turn the PDF device off

eggers_table <- rbind(eggers_table, pub_bias(forest_rma_hyp, "IAD"))

eggers_table

results_bysex_total <- as.data.table(predict(model_iad_multi, digits=3, transf=transf.ilogit))
results_bysex_total


# PD ----------------------------------------------------------------------

data_pd <- dataset[disoder == "PD", ]
table(unique(data_pd[,.(field_citation_value, GBD.Super.Regions)])$GBD.Super.Regions)
table(unique(data_pd[,.(field_citation_value, GBD.Super.Regions)])$GBD.Regions)

model_pd_uni_intercept <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_pd, slab=study_label)
model_pd_uni_intercept

model_pd_uni_sex <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_pd, slab=study_label, mods = ~ m_percent_female)
model_pd_uni_sex 

((model_pd_uni_intercept$sigma2[1] - model_pd_uni_sex$sigma2[1])/model_pd_uni_intercept$sigma2[1])*100
((model_pd_uni_intercept$sigma2[2] - model_pd_uni_sex$sigma2[2])/model_pd_uni_intercept$sigma2[2])*100

model_pd_uni_age <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_pd, slab=study_label, mods = ~ m_age)
model_pd_uni_age 

((model_pd_uni_intercept$sigma2[1] - model_pd_uni_age$sigma2[1])/model_pd_uni_intercept$sigma2[1])*100
((model_pd_uni_intercept$sigma2[2] - model_pd_uni_age$sigma2[2])/model_pd_uni_intercept$sigma2[2])*100

model_pd_uni_1y_recall <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_pd, slab=study_label, mods = ~ cv_1year)
model_pd_uni_1y_recall 

((model_pd_uni_intercept$sigma2[1] - model_pd_uni_1y_recall$sigma2[1])/model_pd_uni_intercept$sigma2[1])*100
((model_pd_uni_intercept$sigma2[2] - model_pd_uni_1y_recall$sigma2[2])/model_pd_uni_intercept$sigma2[2])*100

model_pd_uni_no_hi <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_pd, slab=study_label, mods = ~ cv_not_hi_location)
model_pd_uni_no_hi 

((model_pd_uni_intercept$sigma2[1] - model_pd_uni_no_hi$sigma2[1])/model_pd_uni_intercept$sigma2[1])*100
((model_pd_uni_intercept$sigma2[2] - model_pd_uni_no_hi$sigma2[2])/model_pd_uni_intercept$sigma2[2])*100


model_pd_uni_jbi <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_pd, slab=study_label, mods = ~ JBI_Score)
model_pd_uni_jbi 

((model_pd_uni_intercept$sigma2[1] - model_pd_uni_jbi$sigma2[1])/model_pd_uni_intercept$sigma2[1])*100
((model_pd_uni_intercept$sigma2[2] - model_pd_uni_jbi$sigma2[2])/model_pd_uni_intercept$sigma2[2])*100


model_pd_multi <- model_pd_uni_sex
model_pd_multi

((model_pd_uni_intercept$sigma2[1] - model_pd_multi$sigma2[1])/model_pd_uni_intercept$sigma2[1])*100
((model_pd_uni_intercept$sigma2[2] - model_pd_multi$sigma2[2])/model_pd_uni_intercept$sigma2[2])*100


forest_rma_pd <- rma(yi = yi, vi = vi, data = dataset_overall[disoder == "PD", ] , slab=study_label)
forest(forest_rma_pd)
unique(data_pd$study_label)

pdf(file='FILEPATH/funnelplot_pd.pdf', width = 10, height = 9)
funnel(forest_rma_pd)
dev.off() # Turn the PDF device off

pdf(file='FILEPATH/trimandfill_pd.pdf', width = 10, height = 9)
funnel(trimfill(forest_rma_pd))
dev.off() # Turn the PDF device off

eggers_table <- rbind(eggers_table, pub_bias(forest_rma_pd, "PD"))

eggers_table

results_bysex_total <- as.data.table(predict(model_pd_multi, newmods = rbind(c(0),c(-0.5), c(0.5)), digits=3, transf=transf.ilogit))
results_bysex_total

# CD ----------------------------------------------------------------------
data_cd <- dataset[disoder == "CD", ]

table(data_cd[,.(cv_1year, cv_not_hi_location)])

model_cd_uni_intercept <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_cd, slab=study_label)
model_cd_uni_intercept

model_cd_uni_sex <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_cd, slab=study_label, mods = ~ m_percent_female)
model_cd_uni_sex 

((model_cd_uni_intercept$sigma2[1] - model_cd_uni_sex$sigma2[1])/model_cd_uni_intercept$sigma2[1])*100
((model_cd_uni_intercept$sigma2[2] - model_cd_uni_sex$sigma2[2])/model_cd_uni_intercept$sigma2[2])*100


model_cd_uni_age <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_cd, slab=study_label, mods = ~ m_age)
model_cd_uni_age 

((model_cd_uni_intercept$sigma2[1] - model_cd_uni_age$sigma2[1])/model_cd_uni_intercept$sigma2[1])*100
((model_cd_uni_intercept$sigma2[2] - model_cd_uni_age$sigma2[2])/model_cd_uni_intercept$sigma2[2])*100

model_cd_uni_JBI <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_cd, slab=study_label, mods = ~ JBI_Score)
model_cd_uni_JBI 

((model_cd_uni_intercept$sigma2[1] - model_cd_uni_JBI$sigma2[1])/model_cd_uni_intercept$sigma2[1])*100
((model_cd_uni_intercept$sigma2[2] - model_cd_uni_JBI$sigma2[2])/model_cd_uni_intercept$sigma2[2])*100


model_cd_uni_intercept

forest_rma_cd <- rma(yi = yi, vi = vi, data = dataset_overall[disoder == "CD", ] , slab=study_label)
forest(forest_rma_cd)
unique(data_cd$study_label)

pdf(file='FILEPATH/funnelplot_cd.pdf', width = 10, height = 9)
funnel(forest_rma_cd)
dev.off() # Turn the PDF device off

pdf(file='FILEPATH/trimandfill_cd.pdf', width = 10, height = 9)
funnel(trimfill(forest_rma_cd))
dev.off() # Turn the PDF device off

eggers_table <- rbind(eggers_table, pub_bias(forest_rma_cd, "CD"))

eggers_table

results_total <- as.data.table(predict(model_cd_uni_intercept, digits=3, transf=transf.ilogit))
results_total


# USD simple meta-analysis ------------------------------------------------

data_usd <- dataset[disoder == "USD", ]
model_usd <- rma(yi = yi, vi = vi, data = data_usd, slab=unique_study)
model_usd  -- put as both initial and final model
bic_usd <- BIC(model_usd)
bic_usd  -- put as both initial and final model
results <- as.data.table(predict(model_usd, digits=3, transf=transf.ilogit))
results 

forest_rma_usd <- rma(yi = yi, vi = vi, data = dataset_overall[disoder == "USD", ] , slab=study_label)
forest(forest_rma_usd)
unique(data_usd$study_label)

pdf(file='FILEPATH/funnelplot_usd.pdf', width = 10, height = 9)
funnel(forest_rma_usd)
dev.off() # Turn the PDF device off


# BDD ---------------------------------------------------------------------
data_bdd <- dataset[disoder == "BDD", ]


## Run multi-level regression model
model_bdd_uni_intercept <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_bdd, slab=study_label)
model_bdd_uni_intercept

model_bdd_uni_sex <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_bdd, slab=study_label, mods = ~ m_percent_female)
model_bdd_uni_sex 

((model_bdd_uni_intercept$sigma2[1] - model_bdd_uni_sex$sigma2[1])/model_bdd_uni_intercept$sigma2[1])*100
((model_bdd_uni_intercept$sigma2[2] - model_bdd_uni_sex$sigma2[2])/model_bdd_uni_intercept$sigma2[2])*100

model_bdd_uni_age <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_bdd, slab=study_label, mods = ~ m_age)
model_bdd_uni_age 

((model_bdd_uni_intercept$sigma2[1] - model_bdd_uni_age$sigma2[1])/model_bdd_uni_intercept$sigma2[1])*100
((model_bdd_uni_intercept$sigma2[2] - model_bdd_uni_age$sigma2[2])/model_bdd_uni_intercept$sigma2[2])*100

model_bdd_uni_1y_recall <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_bdd, slab=study_label, mods = ~ cv_1year)
model_bdd_uni_1y_recall 

((model_bdd_uni_intercept$sigma2[1] - model_bdd_uni_1y_recall$sigma2[1])/model_bdd_uni_intercept$sigma2[1])*100
((model_bdd_uni_intercept$sigma2[2] - model_bdd_uni_1y_recall$sigma2[2])/model_bdd_uni_intercept$sigma2[2])*100

model_bdd_uni_no_hi <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_bdd, slab=study_label, mods = ~ cv_not_hi_location)
model_bdd_uni_no_hi 

model_bdd_multi <- rma.mv(random = ~ 1 | unique_study/obs_number,  yi = yi, V = vi, data = data_bdd, slab=study_label, mods = ~ m_age)
model_bdd_multi

forest_rma_bdd <- rma(yi = yi, vi = vi, data = dataset_overall[disoder == "BDD", ] , slab=study_label)
forest(forest_rma_bdd)
unique(data_bdd$study_label)

pdf(file='FILEPATH/funnelplot_bdd.pdf', width = 10, height = 9)
funnel(forest_rma_bdd)
dev.off() # Turn the PDF device off

eggers_table <- rbind(eggers_table, pub_bias(forest_rma_bdd, "BDD"))

eggers_table

results_total <- as.data.table(predict(model_bdd_uni_age, newmods = rbind(c(0)), digits=3, transf=transf.ilogit))
results_total

# Overall forest ----------------------------------------------------------

models <- list(forest_rma_anysd, forest_rma_sd, forest_rma_hyp, forest_rma_pd, forest_rma_cd, forest_rma_bdd)
group_labels  <- c('Somatoform disorders', 'Somatization disorder', 'Hypochondriasis', 'Pain disorder', 'Conversion disorder', 'Body dysmorphic disorder')

# Step 1: Extract study-level data and number of rows needed
study_data <- lapply(models, function(model) {
  data.frame(
    yi = model$yi,
    vi = model$vi,
    sei = sqrt(model$vi),
    slab = model$slab,
    row_type = "study"
  )
})

# Step 2: Combine and assign row numbers
buffer_rows <- 1  # space between groups
rows_needed <- sapply(study_data, nrow) + buffer_rows + 1  # +1 for summary row

total_rows <- sum(rows_needed)
row_indices <- rev(seq_len(total_rows))  # reverse for forest()

# Build one data frame
plot_data <- data.frame()
row_cursor <- total_rows + 1  

row_map <- list()
for (i in seq_along(models)) {
  this_data <- study_data[[i]]
  n <- nrow(this_data)
  
  # Position for this group's studies
  row_pos <- (row_cursor - 1):(row_cursor - n)
  row_cursor <- min(row_pos) - 1
  
  # Store mapping
  row_map[[i]] <- list(study_rows = row_pos, summary_row = row_cursor)
  
  # Add group info
  this_data$row <- row_pos
  this_data$group <- group_labels[i]
  
  plot_data <- rbind(plot_data, this_data)
  
  # Leave buffer row
  row_cursor <- row_cursor - buffer_rows
}

transf.percent <- function(x) {
  100 * transf.ilogit(x)
}


pdf(file='FILEPATH/forest.pdf', width = 10, height = 20)

# Step 3: Plot forest
forest(plot_data$yi, vi = plot_data$vi,
       transf = transf.percent,
       slab = plot_data$slab,
       xlab = "Prevalence (%)",
 
       alim = c(0, 17),

       rows = plot_data$row,
header = c("Study", "Prevalence [95% CI]"))


# Step 4: Add group-level summaries
for (i in seq_along(models)) {
  res <- models[[i]]
  rmap <- row_map[[i]]
  
    addpoly(res, row = rmap$summary_row-0.5, transf = transf.percent, clip = c(0, 0.17), mlab = "")
  
  # Add group label to the left of the summary
  text(-19, rmap$summary_row-0.5, group_labels[i], pos = 4, font = 2)
}

dev.off() # Turn the PDF device off


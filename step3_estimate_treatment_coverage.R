#############################################################################
### Author: Damian Santomauro
### Purpose: Estimates coverage adjusted treatment-effect for anxiety disorders
#############################################################################

rm(list=ls())
set.seed(65235232)

# Load libraries ----------------------------------------------------------

library(readstata13)
library(data.table)
library(survey)
library(foreign)
lookfor <- function(x){names(survey_data)[grep(tolower(x), tolower(names(survey_data)))]}
tab <- function(x, y=NA){if(is.na(y)){table(survey_data[, get(x)], exclude=NULL)} else {table(survey_data[, get(x)], survey_data[, get(y)], exclude=NULL)}}

# Load in 1997 NMHWS data -------------------------------------------------

survey_data <- data.table(read.dta("/FILEPATH/MHS97.dta"))
survey_data[, `:=` (survey_strata=1, weight = (WEIGHT2/10000))]
for(x in c(paste0("WTJK1", 1:9), paste0("WTJK", 10:30))){
  survey_data[, paste0(x) := get(x)/10000]
}

# Code therapies -------------------------------------------

survey_data[, `:=` (d_cbt = "", d_psychodynamic = "", d_supportive = "", d_meds = "")]
survey_data[R9_4 == "Yes" & (R7_7 == "6-10 days" | R7_7 == "More than 10 days" | R7_6 == "6-10 days" | R7_6 == "More than 10 days"), `:=` (d_cbt = "cbt")]
survey_data[R9_3 == "Yes" & (R7_7 == "6-10 days" | R7_7 == "More than 10 days" | R7_6 == "6-10 days" | R7_6 == "More than 10 days"), `:=` (d_psychodynamic = "psychodynamic")]
survey_data[R9_5 == "Yes" & (R7_7 == "6-10 days" | R7_7 == "More than 10 days" | R7_6 == "6-10 days" | R7_6 == "More than 10 days" | R7_10 == "More than 5 days"), `:=` (d_supportive = "supportive")]
survey_data[R9_2 == "Yes", `:=` (d_meds = "meds")]

# Apply therapy heirarchy 
survey_data[d_cbt == "cbt", `:=` (d_psychodynamic = "", d_supportive = "")]
survey_data[d_psychodynamic == "psychodynamic", `:=` (d_supportive = "")]

survey_data[, `:=` (all_treatment = paste0(d_cbt, d_psychodynamic, d_supportive, d_meds))]

# Code overall 30 day Anxiety disorders variable --------------------------

# GAD 
survey_data[!is.na(D300_02C), gad_month_dsm := "No"]
survey_data[D300_02C == "Within last 2 wks" | D300_02C == "2 wks - 1 mth ago" , gad_month_dsm := "Yes"]
# OCD 
survey_data[!is.na(D300_3C), ocd_month_dsm := "No"]
survey_data[D300_3C == "Within last 2 wks"| D300_3C == "2 wks - 1 mth ago" , ocd_month_dsm := "Yes"]
# Social phobia 
survey_data[!is.na(D300_23C), socialph_month_dsm := "No"]
survey_data[D300_23C == "Within last 2 wks" | D300_23C == "2 wks - 1 mth ago" , socialph_month_dsm := "Yes"]
# PTSD
survey_data[!is.na(D309_81C), ptsd_month_dsm := "No"]
survey_data[D309_81C == "Within last 2 wks" | D309_81C == "2 wks - 1 mth ago" , ptsd_month_dsm := "Yes"]
# panic with agora recency
survey_data[!is.na(D300_21C), panicWithagora_month_dsm := "No"]
survey_data[D300_21C == "Within last 2 wks" | D300_21C == "2 wks - 1 mth ago" , panicWithagora_month_dsm := "Yes"]
# panic w/out agora recency 
survey_data[!is.na(D300_01C ), panicW_out_agora_month_dsm := "No"]
survey_data[D300_01C  == "Within last 2 wks" | D300_01C  == "2 wks - 1 mth ago" , panicW_out_agora_month_dsm := "Yes"]
# Agoraphobia w/out history of panic disorder - recency
survey_data[!is.na(D300_22C), agoraW_opanic_month_dsm := "No"]
survey_data[D300_22C  == "Within last 2 wks" | D300_22C  == "2 wks - 1 mth ago" , agoraW_opanic_month_dsm := "Yes"]

# Estimating overall anxiety disorders 
survey_data[gad_month_dsm == "No" | ocd_month_dsm == "No" | socialph_month_dsm == "No" | ptsd_month_dsm == "No" | panicWithagora_month_dsm == "No" 
            | panicW_out_agora_month_dsm == "No"| agoraW_opanic_month_dsm == "No" , anx_30d := "No"]
survey_data[gad_month_dsm == "Yes" | ocd_month_dsm == "Yes" | socialph_month_dsm == "Yes" | ptsd_month_dsm == "Yes" | panicWithagora_month_dsm == "Yes" 
            | panicW_out_agora_month_dsm == "Yes"| agoraW_opanic_month_dsm == "Yes" , anx_30d := "Yes"]

# Obtain detailed medication data from 2007 NSMHWB ---------------------------

drug_classes <- fread('/FILEPATH/nsmhw2007_drugs.csv')
drug_classes[, Class := tolower(Class)]
antidepressants <- c("noradrenalin reuptake inhibitors", "tetracyclic antidepressants", "sari", "snri", "ssri", "monoamine oxidase inhibitors", "tricyclic antidepressants")
drug_classes[Class %in% antidepressants, Class := "antidepressants"]

survey_data_2007 <- merge(data.table(read.dta('/FILEPATH/mhw07bh.dta')), 
                          data.table(read.dta("/FILEPATH/mhw07bp.dta")) , by="ABSHID") 
survey_data_2007[, `:=` (survey_strata=1, weight = as.numeric(MHSFINWT))]

survey_data_2007[, all_meds := paste(MEDCDA, MEDCDB, MEDCDC, MEDCDD, MEDCDE)]

for(m in drug_classes[Interest == 1, unique(Medication)]){
  class_label <- gsub(" ", "_", drug_classes[Medication == m, Class])
  survey_data_2007[grepl(m, all_meds), paste0("d_", class_label) := class_label]
  survey_data_2007[is.na(get(paste0("d_", class_label))), paste0("d_", class_label) := ""]
}

setnames(survey_data_2007, "d_beta_blockers", "d_betablocker")

survey_data_2007[, `:=` (all_treatment = paste0(d_benzodiazepines, d_antipsychotics, d_antidepressants, d_betablocker, d_antihistamine, d_lithium))]

survey_data_2007[, `:=` (anx_30d = ifelse(D_ANXH30 == "Lifetime Anxiety Disorder with 30 day symptoms", 1, 0))]

# Set survey design #
design_2007 <- svydesign(id = ~ABSHID, type = "JK1", repweights =c(paste0("WPM010", 1:9), paste0("WPM01", 10:60)), weights = ~weight, check.strata=F, data = survey_data_2007)
design_2007_anx_meds <- subset(design_2007, anx_30d == 1 & all_treatment != "")

# Estimate treatment coverage #
treatment_coverage <- data.frame(svymean(~all_treatment, design_2007_anx_meds))
treatment_coverage <- data.table(cbind(treatment = row.names(treatment_coverage), treatment_coverage))
treatment_coverage[, treatment := gsub("all_treatment", "", treatment)]
treatment_coverage <- treatment_coverage[treatment != ""]

# Bring in treatment effects
treatment_effects <-fread("/FILEPATH/treatment_effects_predicted.csv")
FCOT <- treatment_effects[d_cognitive_therapy == 1 & d_behavioural_therapy == 1 & d_antidepressants == 1 & 
                            d_antipsychotics == 0 & d_azapirones == 0 & d_anticonvulsants == 0 & d_benzodiazepines == 0 & 
                            d_neurokinin == 0 & d_betablocker == 0 & d_lithium == 0 & d_antihistamine == 0 & 
                            d_cycloserine == 0 & d_psychodynamic == 0 & d_supportive == 0,]


# Remove unwanted treatments #
remove_treatments <- names(treatment_effects)[names(treatment_effects) %like% "d_"][!(names(treatment_effects)[names(treatment_effects) %like% "d_"] %in% names(survey_data_2007)[names(survey_data_2007) %like% "d_"])]
for(r in remove_treatments){
  treatment_effects <- treatment_effects[get(r) == 0, ]
  treatment_effects[,paste0(r) := NULL]
}

for(t in names(treatment_effects)[names(treatment_effects) %like% "d_"]){
  t_label <- gsub("d_", "", t)
  treatment_effects[, paste(t) := ifelse(get(t) == 1, t_label, "")]
}

treatment_effects[, `:=` (all_treatment = paste0(d_benzodiazepines, d_antipsychotics, d_antidepressants, d_betablocker, d_antihistamine, d_lithium))]
treatment_effects <- treatment_effects[,.(all_treatment, effect = Y_mean, effect_se = (Y_mean_hi - Y_mean_lo) / (qnorm(0.975, 0, 1)*2))]
treatment_effects <- treatment_effects[all_treatment != "", ]

## Create draws ##
coolbeta <- function(x, se){
  if(se == 0){
    rep(x, 1000)
  } else {
    a <- x*(x-x^2-se^2)/se^2
    b <- a*(1-x)/x
    rbeta(1000, a, b)
  }
}

for(t in treatment_coverage$treatment){
  draws <- coolbeta(treatment_coverage[treatment == t, mean], treatment_coverage[treatment == t, SE])
  draws <- data.table(treatment = t, draw = paste0("draw_", 1:1000), coverage = draws)
  if(t == treatment_coverage$treatment[1]){
    treatment_coverage_draws <- draws
  } else {
    treatment_coverage_draws <- rbind(treatment_coverage_draws, draws)
  }
}

for(t in treatment_coverage$treatment){
  draws <- rnorm(1000, treatment_effects[all_treatment == t, effect], treatment_effects[all_treatment == t, effect_se])
  draws <- data.table(treatment = t, draw = paste0("draw_", 1:1000), effect = draws)
  if(t == treatment_coverage$treatment[1]){
    treatment_effects_draws <- draws
  } else {
    treatment_effects_draws <- rbind(treatment_effects_draws, draws)
  }
}

final <- merge(treatment_coverage_draws, treatment_effects_draws, by=c("treatment", "draw"))
final[, weighted_effect := coverage*effect]
final[, total_effect := sum(weighted_effect), by="draw"]
final <- unique(final[,.(draw, total_effect)])
final[, `:=` (mean = mean(total_effect), lower = quantile(total_effect, 0.025), upper = quantile(total_effect, 0.975))]
final[, `:=` (se = (upper - lower) / (2*qnorm(0.975, 0, 1)))]
unique(final[,.(mean, lower, upper, se)])

# Create treatment coverage variables for 1997 NSMHWB ---------------------

# Bring in treatment effects
treatment_effects <-fread("/FILEPATH/treatment_effects_predicted.csv")
treatment_effects[, d_cbt := ifelse(d_cognitive_therapy == 1 & d_behavioural_therapy == 1, 1, 0)]
treatment_effects[d_cbt == 1, `:=` (d_cognitive_therapy = 0, d_behavioural_therapy = 0)]
treatment_effects[, `:=`(n_treatments = apply(.SD, 1, sum)), .SDcols=names(treatment_effects)[names(treatment_effects) %like% "d_"]]
treatment_effects <- treatment_effects[n_treatments == 1,.(d_cbt, d_psychodynamic, d_supportive, Y_mean, Y_mean_lo, Y_mean_hi)]
treatment_effects <- treatment_effects[d_cbt == 1 | d_psychodynamic == 1 | d_supportive == 1, .(d_cbt, d_psychodynamic, d_supportive, d_meds = 0, effect = Y_mean, effect_se = (Y_mean_hi - Y_mean_lo) / (qnorm(0.975, 0, 1)*2))]
treatment_effects <- rbind(treatment_effects, data.table(d_cbt = 0, d_psychodynamic = 0, d_supportive = 0, d_meds = 1, effect = unique(final[,mean]), effect_se = unique(final[,se])))
treatment_effects <- rbind(treatment_effects, treatment_effects[d_cbt == 1, .(d_cbt, d_psychodynamic, d_supportive, d_meds = 1, effect = effect + unique(final[,mean]), effect_se = sqrt(effect_se^2 + unique(final[,se^2])))])
treatment_effects <- rbind(treatment_effects, treatment_effects[d_psychodynamic == 1, .(d_cbt, d_psychodynamic, d_supportive, d_meds = 1, effect = effect + unique(final[,mean]), effect_se = sqrt(effect_se^2 + unique(final[,se^2])))])
treatment_effects <- rbind(treatment_effects, treatment_effects[d_supportive == 1, .(d_cbt, d_psychodynamic, d_supportive, d_meds = 1, effect = effect + unique(final[,mean]), effect_se = sqrt(effect_se^2 + unique(final[,se^2])))])

for(t in names(treatment_effects)[names(treatment_effects) %like% "d_"]){
  t_label <- gsub("d_", "", t)
  treatment_effects[, paste(t) := ifelse(get(t) == 1, t_label, "")]
}

treatment_effects[, `:=` (all_treatment = paste0(d_cbt, d_psychodynamic, d_supportive, d_meds))]
treatment_effects <- treatment_effects[,.(all_treatment, effect, effect_se)]
treatment_effects <- treatment_effects[all_treatment != "", ]

# Estimate treatment coverage for 1997 NSMHWB -----------------------------
survey_data[, `:=` (d_cbt_n = ifelse(d_cbt == "cbt", 1, 0), d_psychodynamic_n = ifelse(d_psychodynamic == "psychodynamic", 1, 0), 
                    d_supportive_n = ifelse(d_supportive == "supportive", 1, 0), d_meds_n = ifelse(d_meds == "meds", 1, 0))]

# Set survey design #
design_1997 <- svydesign(id = ~ABSPID, type = "JK1", repweights =c(paste0("WTJK1", 1:9), paste0("WTJK10", 10:60)), weights = ~weight, check.strata=F, data = survey_data)
design_1997_anx <- subset(design_1997, anx_30d == "Yes")

# Estimate treatment coverage #
treatment_coverage <- data.frame(svymean(~all_treatment, design_1997_anx))
treatment_coverage <- data.table(cbind(treatment = row.names(treatment_coverage), treatment_coverage))
treatment_coverage[, treatment := gsub("all_treatment", "", treatment)]
treatment_coverage <- treatment_coverage[treatment != ""]

# Adjust treatments effect for coverage -------------------------------------------------------------------------

for(t in treatment_coverage$treatment){
  draws <- coolbeta(treatment_coverage[treatment == t, mean], treatment_coverage[treatment == t, SE])
  draws <- data.table(treatment = t, draw = paste0("draw_", 1:1000), coverage = draws)
  if(t == treatment_coverage$treatment[1]){
    treatment_coverage_draws <- draws
  } else {
    treatment_coverage_draws <- rbind(treatment_coverage_draws, draws)
  }
}

for(t in treatment_coverage$treatment){
  draws <- rnorm(1000, treatment_effects[all_treatment == t, effect], treatment_effects[all_treatment == t, effect_se])
  draws <- data.table(treatment = t, draw = paste0("draw_", 1:1000), effect = draws)
  if(t == treatment_coverage$treatment[1]){
    treatment_effects_draws <- draws
  } else {
    treatment_effects_draws <- rbind(treatment_effects_draws, draws)
  }
}

final <- merge(treatment_coverage_draws, treatment_effects_draws, by=c("treatment", "draw"))
final[, weighted_effect := coverage*effect]
final[, total_effect := sum(weighted_effect), by="draw"]
final <- unique(final[,.(draw, total_effect)])
final[, `:=` (mean = mean(total_effect), lower = quantile(total_effect, 0.025), upper = quantile(total_effect, 0.975))]
final[, `:=` (se = (upper - lower) / (2*qnorm(0.975, 0, 1)))]
unique(final[,.(mean, lower, upper, se)])

write.csv(final, "/FILEPATH/treatment_coverage_draws.csv", row.names=F)

## write FCOT draws:

write.csv(treatment_effects_draws[treatment == 'cbtmeds',], "/FILEPATH/treatment_coverage_fcot_draws.csv", row.names=F)

treatment_effects_draws[treatment == 'cbtmeds',.(mean = mean(effect), lower = quantile(effect, 0.025), upper = quantile(effect, 0.975))]

## Summary estimates of service use ----------------------------------------

svyby(~d_cbt_n, ~survey_strata, design_1997_anx, svyciprop, vartype=c("se", "ci"))
svyby(~d_psychodynamic_n, ~survey_strata, design_1997_anx, svyciprop, vartype=c("se", "ci"))
svyby(~d_supportive_n, ~survey_strata, design_1997_anx, svyciprop, vartype=c("se", "ci"))
svyby(~d_meds_n, ~survey_strata, design_1997_anx, svyciprop, vartype=c("se", "ci"))

survey_data_2007[, `:=` (d_antidepressants_n = ifelse(d_antidepressants == "antidepressants", 1, 0), d_antipsychotics_n = ifelse(d_antipsychotics == "antipsychotics", 1, 0),
                         d_antihistamine_n = ifelse(d_antihistamine == "antihistamine", 1, 0), d_benzodiazepines_n = ifelse(d_benzodiazepines == "benzodiazepines", 1, 0), d_lithium_n = ifelse(d_lithium == "lithium", 1, 0))]
design_2007 <- svydesign(id = ~ABSHID, type = "JK1", repweights =c(paste0("WPM010", 1:9), paste0("WPM01", 10:60)), weights = ~weight, check.strata=F, data = survey_data_2007)
design_2007_anx <- subset(design_2007, anx_30d == 1)
svyby(~d_antidepressants_n, ~survey_strata, design_2007_anx, svyciprop, vartype=c("se", "ci"))
svyby(~d_antipsychotics_n, ~survey_strata, design_2007_anx, svyciprop, vartype=c("se", "ci"))
svyby(~d_antihistamine_n, ~survey_strata, design_2007_anx, svyciprop, vartype=c("se", "ci"))
svyby(~d_benzodiazepines_n, ~survey_strata, design_2007_anx, svyciprop, vartype=c("se", "ci"))
svyby(~d_lithium_n, ~survey_strata, design_2007_anx, svyciprop, vartype=c("se", "ci"))


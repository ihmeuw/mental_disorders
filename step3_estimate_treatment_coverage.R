
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

survey_data <- data.table(read.dta("/FILEPATH//MHS97.dta"))
survey_data[, `:=` (survey_strata=1, weight = (WEIGHT2/10000))]
repweights_1997 <- names(survey_data)[names(survey_data) %like% "WTJK"]

for(x in repweights_1997){
  survey_data[, paste0(x) := get(x)/10000]
}
repweights_1997 <- survey_data[,(repweights_1997), with = F]

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

# Bring in treatment effects
treatment_effects <-fread("/FILEPATH//treatment_effects_predicted.csv")
treatment_effects[, d_cbt := ifelse(d_cognitive_therapy == 1 & d_behavioural_therapy == 1, 1, 0)]
treatment_effects[d_cbt == 1, `:=` (d_cognitive_therapy = 0, d_behavioural_therapy = 0)]
treatment_effects[, `:=`(n_treatments = apply(.SD, 1, sum)), .SDcols=names(treatment_effects)[names(treatment_effects) %like% "d_"]]

for_plot <- treatment_effects[n_treatments == 1 & d_cognitive_therapy == 0 & d_behavioural_therapy == 0,]

treatment_effects_single <- treatment_effects[n_treatments == 1,.(d_cbt, d_psychodynamic, d_supportive, d_meds = d_antidepressants, Y_mean, Y_mean_lo, Y_mean_hi)]
treatment_effects_with_meds <- treatment_effects[d_antidepressants == 1 & n_treatments == 2,.(d_cbt, d_psychodynamic, d_supportive, d_meds = d_antidepressants, Y_mean, Y_mean_lo, Y_mean_hi)]
treatment_effects <- rbind(treatment_effects_single, treatment_effects_with_meds)

treatment_effects <- treatment_effects[, .(d_cbt, d_psychodynamic, d_supportive, d_meds, effect = Y_mean, effect_se = (Y_mean_hi - Y_mean_lo) / (qnorm(0.975, 0, 1)*2))]

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
design_1997 <- svrepdesign(data = survey_data, repweights =repweights_1997, type = 'JK1', weights = survey_data$weight)
design_1997_anx <- subset(design_1997, anx_30d == "Yes")

# Estimate treatment coverage #
treatment_coverage <- data.frame(svymean(~all_treatment, design_1997_anx))
treatment_coverage <- data.table(cbind(treatment = row.names(treatment_coverage), treatment_coverage))
treatment_coverage[, treatment := gsub("all_treatment", "", treatment)]
treatment_coverage <- treatment_coverage[treatment != ""]

# Adjust treatments effect for coverage -------------------------------------------------------------------------
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

write.csv(final, "/FILEPATH/treatment_coverage_draws.csv", row.names=F)

## write FCOT draws:
fcot_draws <- rnorm(1000, treatment_effects[all_treatment == "cbtmeds", effect], treatment_effects[all_treatment == "cbtmeds", effect_se])
fcot_draws <- data.table(treatment = 'cbtantidep', draw = paste0("draw_", 1:1000), effect = fcot_draws)
write.csv(fcot_draws, "/FILEPATH/treatment_coverage_fcot_draws.csv", row.names=F)

# make effect sizes figure ------------------------------------------------

cbt_cov <- data.table(svyby(~d_cbt_n, ~survey_strata, design_1997_anx, svyciprop, vartype=c("se", "ci")))
psycodynamic_cov <- data.table(svyby(~d_psychodynamic_n, ~survey_strata, design_1997_anx, svyciprop, vartype=c("se", "ci")))
supportive_cov <- data.table(svyby(~d_supportive_n, ~survey_strata, design_1997_anx, svyciprop, vartype=c("se", "ci")))
meds_cov <- data.table(svyby(~d_meds_n, ~survey_strata, design_1997_anx, svyciprop, vartype=c("se", "ci")))

reporting_cov <- function(x, l, u){return(paste0(sprintf("%.1f", x*100), " (", sprintf("%.1f", l*100), " to ", sprintf("%.1f", u*100), ")"))}

cov_for_plot <- rbind(cbt_cov[,.(Treatment = "CBT", cov = reporting_cov(d_cbt_n, ci_l, ci_u))],
                      psycodynamic_cov[,.(Treatment = "Psychodynamic therapy", cov = reporting_cov(d_psychodynamic_n, ci_l, ci_u))],
                      meds_cov[,.(Treatment = "Medication", cov = reporting_cov(d_meds_n, ci_l, ci_u))],
                      supportive_cov[,.(Treatment = "Supportive therapy", cov = reporting_cov(d_supportive_n, ci_l, ci_u))])
cov_for_plot[, cov := gsub("\\.", "·", cov)]
cov_for_plot[, position := 5 - seq_len(.N)]
cov_for_plot[Treatment == "Medication", cov := paste0(cov, "*")]
cov_for_plot[Treatment == "Psychodynamic therapy", cov := paste0(cov, "**")]
cov_for_plot[Treatment == "Supportive therapy", cov := paste0(cov, "***")]
                      
for_plot[, Treatment := ""]
for(t in names(for_plot)[names(for_plot) %like% "d_"]){
  t_label <- gsub("d_", "", t)
  for_plot[, Treatment := ifelse(get(t) == 1, t_label, Treatment)]
}

for_plot <- for_plot[order(Y_mean),]
for_plot[, Treatment := paste0(toupper(substr(Treatment, 1, 1)), substr(Treatment, 2, nchar(Treatment)))]
for_plot[Treatment == "Cbt", Treatment := "Cognitive behavioural therapy"]
for_plot[Treatment == "Psychodynamic", Treatment := "Psychodynamic therapy"]
for_plot[Treatment == "Supportive", Treatment := "Supportive therapy"]

for_plot[, ord := 5-seq_len(.N)]

library(ggplot2)
box_plot <- ggplot(for_plot, aes(reorder(Treatment, ord))) +
  geom_boxplot(aes(ymin = Y_mean_lo, lower = Y_mean_lo_fe, middle = Y_mean, upper = Y_mean_hi_fe, ymax = Y_mean_hi), stat = "identity") +
  geom_hline(yintercept = 0, linetype='dashed') +
  scale_y_continuous(breaks = seq(-1.2, 1, 0.2), labels = sprintf("%.1f", seq(-1.2, 1, 0.2)), limits = c(-1.2, 1))+
   geom_text(x=4.5, y=0.7, label="Coverage % (95% UI)", size = 5) +
  theme_classic() +
  ylab("Effect size") + 
  xlab("Treatment") +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold")) +
   coord_flip()
for(t in cov_for_plot$Treatment){
  box_plot <- box_plot + geom_text(x=cov_for_plot[Treatment == t, position], y=0.7, label=cov_for_plot[Treatment == t, cov], size = 4)
}

box_plot
ggsave(box_plot, filename="/FILEPATH/figure_2.pdf", width = 10, height = 5)

## FCOT summary
treatment_effects_with_meds[d_cbt == 1 & d_meds == 1, ]


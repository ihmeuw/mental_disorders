#####################################################################################
###
### Adjust anxiety disorder severity splits to account for treatemnt
###
#####################################################################################

rm(list = ls())

library(data.table)
library(msm)
library(readstata13)
library(ggplot2)
library(openxlsx)
source("/FILEPATH/mr_brt_functions.R")
source("/FILEPATH/get_covariate_estimates.R")
source("/FILEPATHr/get_draws.R")
source("/FILEPATH/get_location_metadata.R")
source("/FILEPATH/get_population.R")
library(mrbrt001, lib.loc = '/FILEPATH/')
rlogit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}


set.seed(645823) 

##### Bring in draws of treatment effects adjusted for coverage #####

treatment_effect <- fread('/FILEPATH/treatment_coverage_draws.csv')[,.(draw, pred_treat = total_effect)]
treatment_effect[, draw := paste0("draw_", as.numeric(gsub("draw_", "", draw))-1)]

FCOT <- fread('/FILEPATH/treatment_coverage_fcot_draws.csv')[,.(draw, fcot = effect)]
FCOT[, draw := paste0("draw_", as.numeric(gsub("draw_", "", draw))-1)]

##### Pull anxiety DWs from AMHS 1997 survey #####
anxiety_dws <- fread("/FILEPATH/anxiety_dw_distributions.csv")
anxiety_dws[, dw_t := dw_ds]

## conduct standard SF-12 to DW model and estimate an sf-12 score for each person attributable to anxiety disorders
dw_map_model <- py_load_object(filename = "/FILEPATH/sf12_to_dw.pkl", pickle = "dill")

dw_map <- data.table(intercept = 1, sf12_c = seq(40, 120, 0.01)-120)
predict_data <- MRData()
predict_data$load_df(data = dw_map, col_covs=as.list(dw_map_model$cov_names))

beta_samples <- mrbrt001::core$other_sampling$sample_simple_lme_beta(1000L, dw_map_model)
gamma_outer_samples <- matrix(rep(dw_map_model$gamma_soln, each = 1000L), nrow = 1000L)
draws <- dw_map_model$create_draws(predict_data, beta_samples = beta_samples, gamma_samples = gamma_outer_samples, random_study = T)

dw_map$dw_logit_mean <- apply(draws, 1, function(x) mean(x))
dw_map$dw_logit_lo <- apply(draws, 1, function(x) quantile(x, 0.025))
dw_map$dw_logit_hi <- apply(draws, 1, function(x) quantile(x, 0.975))
dw_map[, `:=` (sf12 = sf12_c + 120, dw_logit_se = (dw_logit_hi - dw_logit_lo)/3.92)]

draws <- data.table(draws)
names(draws) <- paste0("draw_", as.numeric(gsub("V", "", names(draws)))-1)
dw_map <- cbind(dw_map, draws)
dw_map <- melt.data.table(dw_map, id.vars = names(dw_map)[!(names(dw_map) %like% "draw")], value.name="dw_logit", variable.name="draw")
dw_map[, dw_t := rlogit(dw_logit)]

remap_dw <- function(x){
  draw_specific_anxiety_dws <- anxiety_dws[draw == x, ]
  s<-sort(dw_map[draw == x,dw_t])
  draw_specific_anxiety_dws[,closest_dw:=s[findInterval(dw_t,s,all.inside = T)]]
  draw_specific_anxiety_dws<-merge(draw_specific_anxiety_dws,dw_map[draw == x,.(sf = sf12, dw_t)],by.x="closest_dw",by.y="dw_t")
}

anxiety_dws_list <- lapply(unique(anxiety_dws$draw), remap_dw)
anxiety_dws <- Reduce(function(x, y){rbind(x, y)}, anxiety_dws_list)
rm(anxiety_dws_list) # to save memory

## Merge with treatment estimates ##
anxiety_dws <- merge(anxiety_dws, treatment_effect, by = "draw")
anxiety_dws <- merge(anxiety_dws, FCOT, by = "draw")

# calculate the SD of sf for each draw
anxiety_dws[,sf_sd:=sd(sf),by=c("draw")]

## Adjust sf-12 score by treatment effect
anxiety_dws[,sf_adj:=sf+(sf_sd*pred_treat)]
anxiety_dws[,sf_fcot:=sf_adj-(sf_sd*fcot)]

# Map back from SF to DW
remap_sf <- function(x){
  draw_specific_anxiety_dws <- anxiety_dws[draw == x, ]
  s<-sort(dw_map[draw == x,sf12])
  ## for adj_dw
  draw_specific_anxiety_dws[,closest_sf:=s[findInterval(sf_adj,s,all.inside = T)]]
  draw_specific_anxiety_dws<-merge(draw_specific_anxiety_dws,dw_map[draw == x,.(sf = sf12, adj_dw = dw_t)],by.x="closest_sf",by.y="sf")
  ## for FCOT
  draw_specific_anxiety_dws[,closest_sf:=s[findInterval(sf_fcot,s,all.inside = T)]]
  draw_specific_anxiety_dws<-merge(draw_specific_anxiety_dws,dw_map[draw == x,.(sf = sf12, fcot_dw = dw_t)],by.x="closest_sf",by.y="sf")
}

anxiety_dws_list <- lapply(unique(anxiety_dws$draw), remap_sf)
anxiety_dws <- Reduce(function(x, y){rbind(x, y)}, anxiety_dws_list)
rm(anxiety_dws_list) # to save memory

anxiety_dws[dw_t < 0, dw_t := 0]
anxiety_dws[dw_t > 1, dw_t := 1]
anxiety_dws[adj_dw < 0, adj_dw := 0]
anxiety_dws[adj_dw > 1, adj_dw := 1]
anxiety_dws[fcot_dw < 0, fcot_dw := 0]
anxiety_dws[fcot_dw > 1, fcot_dw := 1]

anxiety_dws[sf_adj > 120, adj_dw := 0]
anxiety_dws[sf_fcot > 120, fcot_dw := 0]

# Start plot
anxiety_dws[, `:=` (mean_dw = median(dw_t), mean_adj = median(adj_dw)), by = c("id")]
hist_data <- unique(rbind(anxiety_dws[,.(id, Scenario = "Observed", dw_val = mean_dw)], anxiety_dws[,.(id, Scenario = "No treatment", dw_val = mean_adj)]))
hist_data[dw_val < 0, dw_val := 0]
hist_data[dw_val > 1, dw_val := 1]


anxiety_dws[, `:=` (mean_dw = NULL, mean_adj = NULL)]
hist_data$Scenario <- factor(hist_data$Scenario, levels = c("Observed", "No treatment"))

# Average DWs before and after
hist_data[Scenario == "Observed", mean(dw_val)]
hist_data[Scenario == "Observed", sd(dw_val)]
hist_data[Scenario == "No treatment", mean(dw_val)]
hist_data[Scenario == "No treatment", sd(dw_val)]

### Conduct severity analysis ###
cause <- "mental_anxiety"

# load severity cutoffs
data <- fread("/FILEPATH/3b_meps_severity_cutoffs.csv")

## make 1000 distribution variables
data <- data[yld_cause==cause,]

data <- melt.data.table(data, id.vars = names(data)[!(names(data) %like% "MID")], value.name = "MID", variable.name = "draw")
data <- data[draw !="MID_calc",]
data[, `:=` (MID_type = draw)]
data[, `:=` (draw = substring(gsub("MID", "draw_", draw), 1, nchar(gsub("MID", "draw_", draw))-1))]

data <- merge(data[grepl("a", MID_type),.(yld_cause, healthstate_id, severity, draw, MID_lower = MID)],
              data[grepl("b", MID_type),.(yld_cause, healthstate_id, severity, draw, MID_upper = MID)],
              by = c("yld_cause", "healthstate_id", "severity", "draw"))

# finish plot
hist_data$Scenario <- relevel(hist_data$Scenario, ref = "No treatment")
plot_lines <- data.table(severity = c("Asymptomatic", "Mild", "Moderate"), cutoff = c(-0.02, data[severity == 1, mean(MID_upper)], data[severity == 2, mean(MID_upper)]))
plot <- ggplot(hist_data[, .(Scenario, dw_val = ifelse(dw_val == 0, -0.02, dw_val))], aes(x=dw_val, fill = Scenario, color = Scenario)) +
  geom_histogram(position="identity", alpha=0.3) +
  geom_vline(data=plot_lines, aes(xintercept=cutoff), linetype="dashed") + 
  scale_x_continuous(expand = c(0, 0), breaks = c(-0.02, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 170, 10))+
  annotate("text", x = plot_lines[severity == "Asymptomatic", cutoff-0.02], y = -8, label = "Asymp.", size = 3) + 
  annotate("text", x = mean(c(plot_lines[severity == "Asymptomatic", cutoff], plot_lines[severity == "Mild", cutoff])), y = -8, label = "Mild", size = 3) + 
  annotate("text", x = mean(c(plot_lines[severity == "Mild", cutoff], plot_lines[severity == "Moderate", cutoff])), y = -8, label = "Moderate", size = 3) + 
  annotate("text", x = mean(c(plot_lines[severity == "Moderate", cutoff], 1)), y = -8, label = "Severe", size = 3) + 
  coord_cartesian(ylim = c(0, 170), clip = "off")+
  ylab("Frequency") +
  xlab("\nDisability Weight") +
  theme(legend.position="top")
plot
ggsave(plot, filename="/FILEPATH/figure_2.pdf", width = 10, height = 6)

# continue analysis
anxiety_dws <- merge(anxiety_dws, data, by = "draw", allow.cartesian=TRUE)

anxiety_dws[severity == 0, case_severity := ifelse(dw_t <= 0, 1, 0)]
anxiety_dws[severity == 0, case_severity_adj := ifelse(adj_dw <= 0, 1, 0)]
anxiety_dws[severity == 0, case_severity_fcot := ifelse(fcot_dw <= 0, 1, 0)]

anxiety_dws[severity != 0, case_severity := ifelse(dw_t > MID_lower & dw_t <= MID_upper, 1, 0)]
anxiety_dws[severity != 0, case_severity_adj := ifelse(adj_dw > MID_lower & adj_dw <= MID_upper, 1, 0)]
anxiety_dws[severity != 0, case_severity_fcot := ifelse(fcot_dw > MID_lower & fcot_dw <= MID_upper, 1, 0)]

anxiety_dws[, total_cases := max(seq_len(.N)), by = c("draw", "severity")]
anxiety_dws[, sev_cases := sum(case_severity), by = c("draw", "severity")]
anxiety_dws[, sev_cases_adj := sum(case_severity_adj), by = c("draw", "severity")]
anxiety_dws[, sev_cases_fcot := sum(case_severity_fcot), by = c("draw", "severity")]
anxiety_dws[, sev_prop := sev_cases / total_cases]
anxiety_dws[, sev_prop_adj := sev_cases_adj / total_cases]
anxiety_dws[, sev_prop_fcot := sev_cases_fcot / total_cases]
anxiety_dws <- unique(anxiety_dws[,.(yld_cause, healthstate_id, severity, sev_prop, sev_prop_adj, sev_prop_fcot, draw)])

weights <- fread("/FILEPATH/dw_full.csv")
weights <- melt.data.table(weights, id.vars = c("hhseqid", "healthstate_id", "healthstate"), variable.name = "draw", value.name = "dw")
weights[, draw := gsub("draw", "draw_", draw)]

anxiety_dws <- merge(anxiety_dws, weights, by = c("healthstate_id", "draw"))
anxiety_dws[severity == 0, dw := 0]
anxiety_dws[, `:=` (dw_mean = sev_prop*dw, dw_mean_adj = sev_prop_adj * dw, dw_mean_fcot = sev_prop_fcot * dw)]
anxiety_dws[, `:=` (dw_mean = sum(dw_mean), dw_mean_adj = sum(dw_mean_adj), dw_mean_fcot = sum(dw_mean_fcot)), by = "draw"]

severity_table <- rbind(anxiety_dws[,.(Scenario = "Observed", mean= mean(sev_prop), lower = quantile(sev_prop, 0.025), upper = quantile(sev_prop, 0.975)), by = 'severity'],
                        anxiety_dws[,.(Scenario = "No treatment", mean= mean(sev_prop_adj), lower = quantile(sev_prop_adj, 0.025), upper = quantile(sev_prop_adj, 0.975)), by = 'severity'],
                        anxiety_dws[,.(Scenario = "FCOT", mean= mean(sev_prop_fcot), lower = quantile(sev_prop_fcot, 0.025), upper = quantile(sev_prop_fcot, 0.975)), by = 'severity'])

severity_table[, Sequela := ifelse(severity == 0, "Asymptomatic", ifelse(severity == 1, "Mild", ifelse(severity == 2, "Moderate", "Severe")))]
severity_table[, `:=` (Proportion = paste0(round(mean*100, 1), " (", round(lower*100, 1), " to ", round(upper*100, 1), ")"))]
severity_table
severity_table <- dcast(severity_table[,.(Scenario, Sequela, Proportion)], ...~Scenario, value.var="Proportion")
write.csv(severity_table[,.(Sequela, Observed, `No treatment`, FCOT)], "/FILEPATH/severity_table.csv", row.names = F)

## Just to show GBD DWs
anxiety_dws[,.(mean= mean(dw), lower = quantile(dw, 0.025), upper = quantile(dw, 0.975)), by = 'severity']

unique(anxiety_dws[, .(dw_mean = mean(dw_mean), lower = quantile(dw_mean, 0.025), upper = quantile(dw_mean, 0.975),
                       dw_mean_adj = mean(dw_mean_adj), lower_adj = quantile(dw_mean_adj, 0.025), upper_adj = quantile(dw_mean_adj, 0.975),
                       dw_mean_fcot = mean(dw_mean_fcot), lower_fcot = quantile(dw_mean_fcot, 0.025), upper_fcot = quantile(dw_mean_fcot, 0.975))])

anxiety_dws[, dw_adj_factor := dw_mean/dw_mean_adj]
anxiety_dws[, dw_fcot_factor := dw_mean/dw_mean_fcot]

adj_factor <- unique(anxiety_dws[,.(draw, dw_adj_factor)])
adj_factor[,.(dw_adj_factor = mean(dw_adj_factor), lower = quantile(dw_adj_factor, 0.025), upper = quantile(dw_adj_factor, 0.975))]

raw_dws <- unique(anxiety_dws[,.(raw_dw = dw_mean, draw)])
fcot_dws <- unique(anxiety_dws[,.(fcot_dw = dw_mean_fcot, draw)])
adj_dws <- unique(anxiety_dws[,.(adj_dw = dw_mean_adj, draw)])

raw_dws[,.(mean(raw_dw), quantile(raw_dw, 0.025), quantile(raw_dw, 0.975))]
adj_dws[,.(mean(adj_dw), quantile(adj_dw, 0.025), quantile(adj_dw, 0.975))]
fcot_dws[,.(mean(fcot_dw), quantile(fcot_dw, 0.025), quantile(fcot_dw, 0.975))]

percent_averted_in_nsmhwb <- merge(raw_dws, adj_dws, by = "draw")
percent_averted_in_nsmhwb[, change := (adj_dw-raw_dw)/adj_dw]
percent_averted_in_nsmhwb[,.(mean(change), quantile(change, 0.025), quantile(change, 0.975))]

##### Predict adjustment factor by HAQI  #####
min_haqi <- fread("/FILEPATH/min_haqi_draws.csv")
adj_factor <- merge(adj_factor, min_haqi[,.(draw, haqi = min_haqi)], by = 'draw')

haqi <- get_covariate_estimates(covariate_id = 1099, location_id = 71, year_id = c(1997), gbd_round_id = 6, decomp_step = "step4")

haqi <- data.table(draw = paste0("draw_", 0:999),haqi = rnorm(1000, haqi$mean_value, (haqi$upper_value - haqi$lower_value)/(qnorm(0.975, 0, 1)*2)))

observed <- merge(haqi, anxiety_dws[, .(draw, sev_prop, severity)], by = 'draw')

no_treatment <- merge(min_haqi[, .(draw, haqi = min_haqi)], anxiety_dws[, .(draw, sev_prop = sev_prop_adj, severity)], by = 'draw')

adj_factor <- rbind(observed, no_treatment)

pred_matrix <- data.table(expand.grid(haqi = seq(0, 100, 0.5), draw = paste0("draw_", 0:999)))
pred_matrix <- merge(pred_matrix, observed[,.(draw, severity)], by = 'draw', allow.cartesian = T)

for(d in paste0("draw_", 0:999)){
  mod_asymp <- lm(sev_prop~haqi,subset = (draw == d & severity == 0), data = adj_factor)
  mod_mild <- lm(sev_prop~haqi,subset = (draw == d & severity == 1), data = adj_factor)
  mod_mod <- lm(sev_prop~haqi,subset = (draw == d & severity == 2), data = adj_factor)
  mod_sev <- lm(sev_prop~haqi,subset = (draw == d & severity == 3), data = adj_factor)
  new <- data.frame(haqi = seq(0, 100, 0.5))
  pred_matrix[draw == d & severity == 0, sev_prop := predict(mod_asymp, newdata = new)]
  pred_matrix[draw == d & severity == 1, sev_prop := predict(mod_mild, newdata = new)]
  pred_matrix[draw == d & severity == 2, sev_prop := predict(mod_mod, newdata = new)]
  pred_matrix[draw == d & severity == 3, sev_prop := predict(mod_sev, newdata = new)]
  asymp_haqi_floor <- no_treatment[draw == d & severity == 0, sev_prop]
  mild_haqi_floor <- no_treatment[draw == d & severity == 1, sev_prop]
  mod_haqi_floor <- no_treatment[draw == d & severity == 2, sev_prop]
  sev_haqi_floor <- no_treatment[draw == d & severity == 3, sev_prop]
  pred_matrix[draw == d & severity == 0 & haqi <= min_haqi[draw == d, round(min_haqi/0.5)*0.5], sev_prop := asymp_haqi_floor]
  pred_matrix[draw == d & severity == 1 & haqi <= min_haqi[draw == d, round(min_haqi/0.5)*0.5], sev_prop := mild_haqi_floor]
  pred_matrix[draw == d & severity == 2 & haqi <= min_haqi[draw == d, round(min_haqi/0.5)*0.5], sev_prop := mod_haqi_floor]
  pred_matrix[draw == d & severity == 3 & haqi <= min_haqi[draw == d, round(min_haqi/0.5)*0.5], sev_prop := sev_haqi_floor]
  if(d %in% paste0("draw_", seq(99, 999, 100))){print(d)}
}

pred_matrix[, sum_prop := sum(sev_prop), by = c("draw", "haqi")]
pred_matrix[, sev_prop := sev_prop / sum_prop] # just to correct some that are 0.999999999999999 due to rounding

line_plot <- copy(pred_matrix)
line_plot[, `:=` (mean = mean(sev_prop), lower = quantile(sev_prop, 0.025), upper = quantile(sev_prop, 0.975)), by = c("haqi", "severity")]
line_plot <- unique(line_plot[,.(HAQI = haqi, Severity = ifelse(severity == 0, "Asymptomatic", ifelse(severity == 1, "Mild", ifelse(severity == 2, "Moderate", "Severe"))), Proportion = mean, lower, upper)])

plot <- ggplot(data=line_plot, aes(x=HAQI, y=Proportion, color = Severity))+
  geom_line(size=2) +
  geom_ribbon(data= line_plot, aes(x=HAQI, ymin=lower, ymax=upper, color = Severity), alpha=.1) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, 10), limits = c(0, 101)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0.0, 1, 0.1),  limits = c(0, 1))+
  ylab("Proportion of anxiety disorder cases") +
  xlab("Healthcare Access Quality Index") +
  theme(axis.line=element_line(colour="black"),
        legend.title=element_blank())
plot
ggsave(plot, filename="/FILEPATH/adjust_factor.pdf", width = 8, height = 6)


# Estimate burden ---------------------------------------------------------
library(dplyr)
location_metadata <- get_location_metadata(35, gbd_round_id = 6, decomp_step = 'step5')
countries <- location_metadata[level == 3, location_id]

prevalence <- get_draws(year_id = 2019, gbd_id_type = "cause_id", gbd_id = 571, source = "como", decomp_step = 'step5', gbd_round_id = 6, version_id = 470, metric_id = 3, sex_id = 3, measure_id = 5, age_group_id = 22, location_id = countries)
prevalence <- melt.data.table(prevalence, id.vars = names(prevalence)[!(names(prevalence) %like% "draw")], value.name="prev", variable.name="draw")

set.seed(645823) # reset seed

haqi_by_country <- get_covariate_estimates(covariate_id = 1099, location_id = countries, year_id = 2019, gbd_round_id = 6, decomp_step = "step4")
haqi_by_country[, (paste0("draw_", 0:999)) := as.list(rnorm(n = 1000, mean = mean_value, sd = (upper_value - lower_value)/(qnorm(0.975, 0, 1)*2))), by = "location_id"]
haqi_by_country <- melt.data.table(haqi_by_country, id.vars = names(haqi_by_country)[!(names(haqi_by_country) %like% "draw")], value.name="haqi", variable.name="draw")

prevalence <- merge(prevalence[,.(location_id, draw, prev)], haqi_by_country[,.(location_id, draw, haqi)], all.x = T, by = c("location_id", "draw"))

prevalence[, haqi := (5*round((haqi*10)/5))/10]

prevalence <- merge(prevalence, pred_matrix, all.x = T, by = c("haqi", "draw"))

prevalence <- merge(prevalence, unique(anxiety_dws[, .(draw, severity, dw)]), all.x = T, by = c("draw", "severity")) # merge with healthstate-specific DWs
prevalence[, adj_dw := sev_prop * dw]
prevalence[, adj_dw := sum(adj_dw), by = c("draw", "location_id")]
prevalence[, (c("sum_prop", "dw")) := NULL]
prevalence[, severity := paste0("sev_prop_", severity)]

prevalence <- dcast(prevalence, ...~severity, value.var="sev_prop")

prevalence <- merge(prevalence, raw_dws, all.x = T, by = "draw") # merge with raw DWs

prevalence <- merge(prevalence, fcot_dws, all.x = T, by = "draw") # merge with fcot DWs

prevalence <- merge(prevalence, adj_dws[,.(draw, worst_dw = adj_dw)], all.x = T, by = "draw") # merge with no-treatment DWs

population <- get_population(age_group_id = 22, sex_id = 3, location_id = countries, year_id = 2019, gbd_round_id = 6, decomp_step = 'step5')

prevalence <- merge(prevalence, population[,.(location_id, population)], all.x = T, by = "location_id")

prevalence <- merge(prevalence, location_metadata[,.(location_id, super_region_name, region_name)], all.x = T, by = "location_id")

best_haqi <- unique(haqi_by_country[mean_value == max(mean_value),location_id])

best_dw <- unique(prevalence[location_id == best_haqi, .(draw, best_dw = adj_dw)])

prevalence <- merge(prevalence, best_dw, all.x = T, by = "draw")

prevalence[, cases := population * prev]

prevalence[, `:=` (raw_yld = raw_dw * cases, adj_yld = adj_dw * cases, worst_yld = worst_dw * cases, best_yld = best_dw * cases, fcot_ylds = fcot_dw * cases)]

# Aggregate to region ------------------------------------------------------

prevalence[, `:=` (cases_region = sum(cases), population_region = sum(population), raw_yld_region = sum(raw_yld), worst_yld_region = sum(worst_yld), adj_yld_region = sum(adj_yld), best_yld_region = sum(best_yld), fcot_yld_region = sum(fcot_ylds)), by = c("draw", "region_name")]
prevalence[, `:=` (cases_super_region = sum(cases), population_super_region = sum(population), raw_yld_super_region = sum(raw_yld), worst_yld_super_region = sum(worst_yld), adj_yld_super_region = sum(adj_yld), best_yld_super_region = sum(best_yld), fcot_yld_super_region = sum(fcot_ylds)), by = c("draw", "super_region_name")]
prevalence[, `:=` (cases_global = sum(cases), population_global = sum(population), raw_yld_global = sum(raw_yld), worst_yld_global = sum(worst_yld), adj_yld_global = sum(adj_yld), best_yld_global = sum(best_yld), fcot_yld_global = sum(fcot_ylds)), by = c("draw")]

super_region_data <- unique(prevalence[,.(super_region_name, population = population_super_region, raw_yld = raw_yld_super_region, worst_yld = worst_yld_super_region, adj_yld = adj_yld_super_region, best_yld = best_yld_super_region, fcot_yld = fcot_yld_super_region, draw)])
super_region_data[, `:=` (averted_yld = (worst_yld - adj_yld)/worst_yld, brc_yld = (adj_yld - best_yld)/worst_yld, fcot_yld = (best_yld - fcot_yld)/worst_yld)]
super_region_data[averted_yld < 0, `:=` (averted_yld = 0)]
super_region_data[, `:=` (unavoidable_yld = 1 - averted_yld - brc_yld - fcot_yld)]

super_region_data[,.(averted_yld = round(mean(averted_yld)*100,1), averted_lower = round(quantile(averted_yld, 0.025)*100,1), averted_upper = round(quantile(averted_yld, 0.975)*100, 1)), by = "super_region_name"]
super_region_data[,.(brc_yld = mean(brc_yld), brc_lower = quantile(brc_yld, 0.025), brc_upper = quantile(brc_yld, 0.975)), by = "super_region_name"]
super_region_data[,.(fcot_yld = mean(fcot_yld), fcot_lower = quantile(fcot_yld), fcot_upper = quantile(fcot_yld, 0.975)), by = "super_region_name"]
super_region_data[,.(remaining_yld = mean(worst_yld), remaining_lower = quantile(worst_yld, 0.025), remaining_upper = quantile(worst_yld, 0.975)), by = "super_region_name"]

super_region_data[,.(total_BRC = mean(brc_yld+averted_yld), lower = quantile(brc_yld+averted_yld, 0.025), upper = quantile(brc_yld+averted_yld, 0.975)), by = "super_region_name"]
super_region_data[,.(total_FCOT = mean(fcot_yld+brc_yld+averted_yld), lower = quantile(fcot_yld+brc_yld+averted_yld, 0.025), upper = quantile(fcot_yld+brc_yld+averted_yld, 0.975)), by = "super_region_name"]

super_region_data_summary <- super_region_data[,.(averted_yld = mean(averted_yld), brc_yld = mean(brc_yld), fcot_yld = mean(fcot_yld), unavoidable_yld = mean(unavoidable_yld)), by = "super_region_name"]
super_region_data_summary <- melt.data.table(super_region_data_summary, id.vars = names(super_region_data_summary)[!(names(super_region_data_summary) %like% "yld")], value.name="value", variable.name="YLDs")

super_region_data_summary[, value := value * 100]

super_region_data_summary[, super_region_name := gsub(", ", ",\n", super_region_name)]
super_region_data_summary[, super_region_name := gsub(" and", "\nand", super_region_name)]
super_region_data_summary[, super_region_name := gsub("Saharan", "Saharan\n", super_region_name)]
super_region_data_summary[YLDs == "averted_yld", YLDs := "Averted"]
super_region_data_summary[YLDs == "brc_yld", YLDs := "Avoidable: BRC"]
super_region_data_summary[YLDs == "fcot_yld", YLDs := "Avoidable: FCOT"]
super_region_data_summary[YLDs == "unavoidable_yld", YLDs := "Remaining"]

gg<-ggplot(data=super_region_data_summary)+
  geom_bar(aes(y=value,x=super_region_name,fill=YLDs),stat="identity")+
  theme_bw()+
  scale_y_continuous(expand=c(0,0))+
  ylab("Percent of total burden")+
  xlab("GBD super region")+
  scale_fill_manual(values=c("darkorange", "#9BC53D","#00A5CE","#006594"))
gg
ggsave(gg, filename="/FILEPATH/figure_3.pdf", width = 10, height = 6)

global_data <- unique(prevalence[,.(population = population_global, raw_yld = raw_yld_global, worst_yld = worst_yld_global, adj_yld = adj_yld_global, best_yld = best_yld_global, fcot_yld = fcot_yld_global, draw)])
global_data[, `:=` (averted_yld = (worst_yld - adj_yld)/worst_yld, brc_yld = (adj_yld - best_yld)/worst_yld, fcot_yld = (best_yld - fcot_yld)/worst_yld)]
global_data[, `:=` (unavoidable_yld = 1 - averted_yld - brc_yld - fcot_yld)]

global_data[,.(averted_yld = round(mean(averted_yld)*100, 1), averted_lower = round(quantile(averted_yld, 0.025)*100, 1), averted_upper = round(quantile(averted_yld, 0.975)*100, 1))]
global_data[,.(brc_yld = mean(brc_yld), brc_lower = quantile(brc_yld, 0.025), brc_upper = quantile(brc_yld, 0.975))]
global_data[,.(fcot_yld = mean(fcot_yld), fcot_lower = quantile(fcot_yld), fcot_upper = quantile(fcot_yld, 0.975))]
global_data[,.(remaining_yld = mean(worst_yld), remaining_lower = quantile(worst_yld, 0.025), remaining_upper = quantile(worst_yld, 0.975))]

global_data[,.(total_BRC = round(mean(brc_yld+averted_yld)*100, 1), lower = round(quantile(brc_yld+averted_yld, 0.025)*100, 1), upper = round(quantile(brc_yld+averted_yld, 0.975)*100,1))]
global_data[,.(total_FCOT = round(mean(fcot_yld+brc_yld+averted_yld)*100, 1), lower = round(quantile(fcot_yld+brc_yld+averted_yld, 0.025)*100, 1), upper = round(quantile(fcot_yld+brc_yld+averted_yld, 0.975)*100, 1))]

# Create table with country-level estimates -------------------------------

prevalence[, `:=` (adj_dw_region = adj_yld_region / cases_region)]
prevalence[, `:=` (adj_dw_super_region = adj_yld_super_region / cases_super_region)]
prevalence[, `:=` (adj_dw_global = adj_yld_global / cases_global)]

prevalence[, `:=` (dw_change_percent = 100*(adj_dw - raw_dw)/raw_dw, dw_change_percent_region = 100*(adj_dw_region-raw_dw)/raw_dw, dw_change_percent_super_region = 100*(adj_dw_super_region-raw_dw)/raw_dw, dw_change_percent_global = 100*(adj_dw_global-raw_dw)/raw_dw)]

prevalence[, `:=` (haqi_region = sum(haqi * population)/population_region), by = c("draw", "region_name")]
prevalence[, `:=` (haqi_super_region = sum(haqi * population)/population_super_region), by = c("draw", "super_region_name")]
prevalence[, `:=` (haqi_global = sum(haqi * population)/population_global), by = "draw"]

pres_results_1 <- function(m){return(paste0(sprintf("%.1f", mean(m)), " (", sprintf("%.1f", quantile(m, 0.025)), "â", sprintf("%.1f", quantile(m, 0.975)), ")"))}
pres_results_3 <- function(m){return(paste0(sprintf("%.3f", mean(m)), " (", sprintf("%.3f", quantile(m, 0.025)), "â", sprintf("%.3f", quantile(m, 0.975)), ")"))}

region_table <- unique(prevalence[,.(location_name = region_name, haqi = pres_results_1(haqi_region), raw_dw = pres_results_3(raw_dw), adj_dw = pres_results_3(adj_dw_region), dw_change_percent = pres_results_1(dw_change_percent_region)), by = "region_name"])
super_region_table <- unique(prevalence[,.(location_name = super_region_name, haqi = pres_results_1(haqi_super_region), raw_dw = pres_results_3(raw_dw), adj_dw = pres_results_3(adj_dw_super_region), dw_change_percent = pres_results_1(dw_change_percent_super_region)), by = "super_region_name"])
global_table <- unique(prevalence[,.(location_name = "Global", haqi = pres_results_1(haqi_global), raw_dw = pres_results_3(raw_dw), adj_dw = pres_results_3(adj_dw_global), dw_change_percent = pres_results_1(dw_change_percent_global))])

region_table <- merge(region_table, location_metadata[,.(location_name, lancet_label, sort_order)], by = "location_name")
super_region_table <- merge(super_region_table, location_metadata[,.(location_name, lancet_label, sort_order)], by = "location_name")
global_table <- merge(global_table, location_metadata[,.(location_name, lancet_label, sort_order)], by = "location_name")

country_table <- unique(prevalence[,.(haqi = pres_results_1(haqi), raw_dw = pres_results_3(raw_dw), adj_dw = pres_results_3(adj_dw), dw_change_percent = pres_results_1(dw_change_percent)), by = "location_id"])
country_table <- merge(country_table, location_metadata[,.(location_id, lancet_label, sort_order)], by = "location_id")

final_table <- rbind(country_table, region_table, super_region_table, global_table, fill = T)
final_table <- final_table[order(sort_order), .(Location = lancet_label, HAQI = gsub("\\.", "Â·", haqi), `Original DW` = gsub("\\.", "Â·", raw_dw), `Adjusted DW` = gsub("\\.", "Â·", adj_dw), `DW % change` = gsub("\\.", "Â·", dw_change_percent))]
final_table[, duplicate := seq_len(.N), by = "Location"]
final_table <- final_table[duplicate == 1,]
final_table[, duplicate := NULL]

write.xlsx(final_table, "/FILEPATH/table_s3.xlsx")

prevalence[, dw_mean := mean(adj_dw), by = "location_id"]
location_metadata[location_id == unique(prevalence[dw_mean == min(dw_mean),location_id]), location_name]

final_table[, num := as.numeric(substring(`Adjusted DW`, 3, 6))]
final_table[num == min(num),]
final_table[num == max(num),]

# Create table with country-level severity splits -------------------------------
setnames(prevalence, paste0("sev_prop_", c(0:3)), c("adj_asymp_prop", "adj_mild_prop", "adj_mod_prop", "adj_sev_prop"))

prevalence[, `:=` (asymp_region = sum(adj_asymp_prop * population)/population_region,
                   mild_region = sum(adj_mild_prop * population)/population_region,
                   mod_region = sum(adj_mod_prop * population)/population_region,
                   sev_region = sum(adj_sev_prop*population)/population_region), by = c("draw", "region_name")]

prevalence[, `:=` (asymp_super_region = sum(adj_asymp_prop * population)/population_super_region,
                   mild_super_region = sum(adj_mild_prop * population)/population_super_region,
                   mod_super_region = sum(adj_mod_prop * population)/population_super_region,
                   sev_super_region = sum(adj_sev_prop*population)/population_super_region), by = c("draw", "super_region_name")]

prevalence[, `:=` (asymp_global = sum(adj_asymp_prop * population)/population_global,
                   mild_global = sum(adj_mild_prop * population)/population_global,
                   mod_global = sum(adj_mod_prop * population)/population_global,
                   sev_global = sum(adj_sev_prop*population)/population_global), by = c("draw")]

pres_results_prop <- function(m){return(paste0(sprintf("%.1f", mean(m)*100), " (", sprintf("%.1f", quantile(m, 0.025)*100), "â", sprintf("%.1f", quantile(m, 0.975)*100), ")"))}

region_table_sev <- unique(prevalence[,.(location_name = region_name, Asymptomatic = pres_results_prop(asymp_region), Mild = pres_results_prop(mild_region), Moderate = pres_results_prop(mod_region), Severe = pres_results_prop(sev_region)), by = "region_name"])
super_region_table_sev <- unique(prevalence[,.(location_name = super_region_name, Asymptomatic = pres_results_prop(asymp_super_region), Mild = pres_results_prop(mild_super_region), Moderate = pres_results_prop(mod_super_region), Severe = pres_results_prop(sev_super_region)), by = "super_region_name"])
global_table_sev <- unique(prevalence[,.(location_name = "Global", Asymptomatic = pres_results_prop(asymp_global), Mild = pres_results_prop(mild_global), Moderate = pres_results_prop(mod_global), Severe = pres_results_prop(sev_global))])

region_table_sev <- merge(region_table_sev, location_metadata[,.(location_name, lancet_label, sort_order)], by = "location_name")
super_region_table_sev <- merge(super_region_table_sev, location_metadata[,.(location_name, lancet_label, sort_order)], by = "location_name")
global_table_sev <- merge(global_table_sev, location_metadata[,.(location_name, lancet_label, sort_order)], by = "location_name")

country_table_sev <- unique(prevalence[,.(Asymptomatic = pres_results_prop(adj_asymp_prop), Mild = pres_results_prop(adj_mild_prop), Moderate = pres_results_prop(adj_mod_prop), Severe = pres_results_prop(adj_sev_prop)), by = "location_id"])
country_table_sev <- merge(country_table_sev, location_metadata[,.(location_id, lancet_label, sort_order)], by = "location_id")

final_table_sev <- rbind(country_table_sev, region_table_sev, super_region_table_sev, global_table_sev, fill = T)
final_table_sev <- final_table_sev[order(sort_order), .(Location = lancet_label, Asymptomatic, Mild, Moderate, Severe)]
final_table_sev[, duplicate := seq_len(.N), by = "Location"]
final_table_sev <- final_table_sev[duplicate == 1,]
final_table_sev[, duplicate := NULL]
final_table_sev[, `:=` (Asymptomatic = gsub("\\.", "Â·", Asymptomatic), Mild = gsub("\\.", "Â·", Mild), Moderate = gsub("\\.", "Â·", Moderate), Severe = gsub("\\.", "Â·", Severe))]

write.xlsx(final_table_sev, "/FILEPATH/table_s4.xlsx")



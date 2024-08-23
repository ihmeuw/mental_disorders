##########################################################
#### Estimation of PAFs of ASD to self-harm mortality ####
##########################################################

library(data.table)
library(ggplot2)
library(openxlsx)
source("/FILEPATH/get_draws.R")
source("/FILEPATH/interpolate.R")
source("/FILEPATH/get_outputs.R")
source("/FILEPATH/get_population.R")
source("/FILEPATH/get_location_metadata.R")
source("/FILEPATH/get_ids.R")
round_c <- function(x,y){sprintf(paste0("%.",y, "f"), x)}
mean_and_UIs <- function(v){ paste0(round_c(mean(v), 1), " (", round_c(quantile(v, 0.025), 1), "-", round_c(quantile(v, 0.975), 1), ")")}
rr_matrix <- fread("/FILEPATH/rr_matrix.csv")
source("/FILEPATH/gbd2021_map.R")

cod_correction_version <- 393
como_version <- 1471
compare_version <- 8016
epi_ages <- c(2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 31, 32, 34, 235, 238, 388, 389)

# Load age-metadata -------------------------------------------------------
age_ids <- get_ids('age_group')[age_group_id %in%  c(2:3, 388:389, 238, 34,  6:20, 30:32, 235),]
suppressWarnings(age_ids[, `:=` (age_start = as.numeric(unlist(strsplit(age_group_name, " "))[1]), age_end = as.numeric(unlist(strsplit(age_group_name, " "))[3])), by = "age_group_id"])
age_ids[age_group_id %in% c(2, 3, 388, 389), `:=` (age_start = 0, age_end = 0)]
age_ids[age_group_id %in% c(238), `:=` (age_start = 1, age_end = 1)]
age_ids[age_start == 95, age_end := 99]
age_ids[, mid_age := (age_start+age_end)/2]

locations <- get_location_metadata(35, release_id = 9)

country_overall_draws <- data.table()
region_draws <- data.table()
super_region_draws <- data.table()
for(super_region in locations[!is.na(super_region_id), unique(super_region_id)]){
  for(region in locations[super_region_id == super_region & !is.na(region_id), unique(region_id)]){
    ## Calculate PAFs
    asd_draws <- get_draws(gbd_id_type = "sequela_id", gbd_id = c(1327:1332), source = "como", measure_id = 5, metric_id = 3, location_id = locations[region_id == region & level == 3, unique(location_id)], year_id = 2021, sex_id = c(1, 2), version_id = como_version, release_id = 9)
    asd_draws <- melt.data.table(asd_draws, id.vars = names(asd_draws)[!(names(asd_draws) %like% "draw")], value.name = "prev", variable.name = "draw")
    asd_draws[, intel_dis := 0]
    asd_draws[sequela_id %in% c(1329:1332), intel_dis := 1]
    asd_draws[, prev := sum(prev), by = c("age_group_id", "location_id", "measure_id", "sex_id", "year_id", "metric_id", "draw", "intel_dis")]
    asd_draws <- unique(asd_draws[,.(age_group_id, location_id, sex_id, year_id, draw, prev, intel_dis)])
    asd_draws <- merge(asd_draws, rr_matrix, by = c("sex_id", "intel_dis", "draw"), all.x = T)
    asd_draws[, exp_rr := prev * rr]
    asd_draws[, sum_exp_rr := sum(exp_rr), by=c("age_group_id", "location_id", "sex_id", "draw")]
    asd_draws[, no_asd := 1-sum(prev), by=c("age_group_id", "location_id", "sex_id", "draw")]
    tmred <- 1
    asd_draws[, paf := (sum_exp_rr + no_asd - tmred) / (sum_exp_rr + no_asd)]
    asd_draws <- unique(asd_draws[,.(age_group_id, sex_id, draw, location_id, year_id, paf)])

    ## Calculate attributable self-harm deaths
    death_draws <- get_draws(gbd_id_type = "cause_id", gbd_id = 718, source = "codcorrect",  measure_id = c(1, 4), metric_id = 1, location_id = locations[region_id == region & level == 3, unique(location_id)], year_id = 2021, sex_id = c(1,2), version_id = cod_correction_version, release_id = 9)
    yll_draws <- melt.data.table(death_draws[measure_id == 4, ], id.vars = names(death_draws)[!(names(death_draws) %like% "draw")], value.name = "ylls", variable.name = "draw")
    death_draws <- melt.data.table(death_draws[measure_id == 1, ], id.vars = names(death_draws)[!(names(death_draws) %like% "draw")], value.name = "deaths", variable.name = "draw")
    asd_draws <- merge(asd_draws, death_draws, by = c("age_group_id", "location_id", "sex_id", "year_id", "draw"), all.x = T)
    asd_draws <- merge(asd_draws, yll_draws, by = c("age_group_id", "location_id", "sex_id", "year_id", "draw"), all.x = T)
    asd_draws[is.na(ylls), ylls := 0] # For age groups where self-harm deaths are not modelled as they are assumed to be 0 deaths
    asd_draws[is.na(deaths), deaths := 0] # For age groups where self-harm deaths are not modelled as they are assumed to be 0 deaths
    asd_draws <- asd_draws[,.(age_group_id, location_id, sex_id, year_id, draw, paf, deaths, ylls)]
    asd_draws[, `:=` (asd_deaths = paf * deaths, asd_ylls = paf * ylls)]

    population <- get_population(run_id = 359, location_id = locations[region_id == region & level == 3, unique(location_id)], age_group_id = epi_ages, year_id = 2021, sex_id = c(1,2), release_id = 9)

    asd_draws <- merge(asd_draws, population, by = c('age_group_id', 'location_id', 'year_id', 'sex_id'))

    asd_bothsex <- copy(asd_draws)
    asd_bothsex[, `:=` (deaths_country = sum(deaths), ylls_country = sum(ylls), asd_deaths_country = sum(asd_deaths), ylls_asd_country = sum(asd_ylls), pop_country = sum(population)), by = c('location_id', "draw")]
    asd_draws[, `:=` (deaths_country = sum(deaths), ylls_country = sum(ylls), asd_deaths_country = sum(asd_deaths), ylls_asd_country = sum(asd_ylls), pop_country = sum(population)), by = c('sex_id', 'location_id', "draw")]

    country_summaries_bothsex <- unique(asd_bothsex[,.(sex_id = 3, location_id, year_id, draw, death_count = deaths_country, yll_count = ylls_country, asd_death_count = asd_deaths_country, asd_yll_count = ylls_asd_country, pop = pop_country)])
    country_summaries_bysex <- unique(asd_draws[,.(sex_id, location_id, year_id, draw, death_count = deaths_country, yll_count = ylls_country, asd_death_count = asd_deaths_country, asd_yll_count = ylls_asd_country, pop = pop_country)])

    country_summaries <- rbind(country_summaries_bothsex, country_summaries_bysex)

    country_summaries[, `:=` (death_rate = death_count / pop, yll_rate = yll_count / pop, asd_death_rate = asd_death_count / pop, asd_yll_rate = asd_yll_count / pop)]
    country_overall_draws <- rbind(country_overall_draws, country_summaries)

    ## Region estimates
    population_region <- get_population(run_id = 359, location_id = region, age_group_id = epi_ages, year_id = 2021, sex_id = c(1,2), release_id = 9)

    asd_draws <- merge(asd_draws, population_region[,.(age_group_id, sex_id, pop_region = population)], by = c("age_group_id", "sex_id"))
    asd_draws[, pop_agg := sum(population), by = c("age_group_id", "sex_id", "draw")]
    asd_draws[, pop_correction := pop_region / pop_agg]
    asd_draws[, `:=` (deaths = sum(deaths), ylls = sum(ylls), asd_deaths = sum(asd_deaths), asd_ylls = sum(asd_ylls)), by = c("age_group_id", "sex_id", "draw")]
    asd_draws[, `:=` (deaths = deaths * pop_correction, ylls = ylls * pop_correction, asd_deaths = asd_deaths * pop_correction, asd_ylls = asd_ylls * pop_correction)]
    asd_draws <- unique(asd_draws[,.(age_group_id, sex_id, location_id = region, draw, deaths, ylls, asd_deaths, asd_ylls, pop_region)])
    region_draws <- rbind(region_draws, asd_draws)
    print(paste0("Finished region ", region))
  }
  population_super_region <- get_population(run_id = 359, location_id = super_region, age_group_id = epi_ages, year_id = 2021, sex_id = c(1,2), release_id = 9)

  region_children <- locations[super_region_id == super_region & !is.na(region_id), unique(region_id)]

  super_region_data <- region_draws[location_id %in% region_children, ]
  super_region_data <-  merge(super_region_data, population_super_region[,.(age_group_id, sex_id, pop_super_region = population)], by = c("age_group_id", "sex_id"))

  super_region_data[, pop_agg := sum(pop_region), by = c("age_group_id", "sex_id", "draw")]
  super_region_data[, pop_correction := pop_super_region / pop_agg]
  super_region_data[, `:=` (deaths = sum(deaths), ylls = sum(ylls), asd_deaths = sum(asd_deaths), asd_ylls = sum(asd_ylls)), by = c("age_group_id", "sex_id", "draw")]
  super_region_data[, `:=` (deaths = deaths * pop_correction, ylls = ylls * pop_correction, asd_deaths = asd_deaths * pop_correction, asd_ylls = asd_ylls * pop_correction)]
  super_region_data <- unique(super_region_data[,.(age_group_id, sex_id, location_id = super_region, draw, deaths, ylls, asd_deaths, asd_ylls, pop_super_region)])
  super_region_draws <- rbind(super_region_draws, super_region_data)
  print(paste0("Finished super region ", super_region))
}

population_global <- get_population(run_id = 359, location_id = 1, age_group_id = epi_ages, year_id = 2021, sex_id = c(1,2), release_id = 9)

global_draws <- copy(super_region_draws)
global_draws <-  merge(global_draws, population_global[,.(age_group_id, sex_id, pop_global = population)], by = c("age_group_id", "sex_id"))
global_draws[, pop_agg := sum(pop_super_region), by = c("age_group_id", "sex_id", "draw")]
global_draws[, pop_correction := pop_global / pop_agg]
global_draws[, `:=` (deaths = sum(deaths), ylls = sum(ylls), asd_deaths = sum(asd_deaths), asd_ylls = sum(asd_ylls)), by = c("age_group_id", "sex_id", "draw")]
global_draws[, `:=` (deaths = deaths * pop_correction, ylls = ylls * pop_correction, asd_deaths = asd_deaths * pop_correction, asd_ylls = asd_ylls * pop_correction)]
global_draws <- unique(global_draws[,.(age_group_id, sex_id, location_id = 1, draw, deaths, ylls, asd_deaths, asd_ylls, pop_global)])

global_draws[, level := 0]
super_region_draws[, level := 1]
region_draws[, level := 2]
setnames(global_draws, "pop_global", "population")
setnames(super_region_draws, "pop_super_region", "population")
setnames(region_draws, "pop_region", "population")
global_draws <- rbind(global_draws, super_region_draws, region_draws)
rm(super_region_draws, region_draws)

########################

global_draws_bothsex <- copy(global_draws)
global_draws_bothsex[, `:=` (deaths = sum(deaths), ylls = sum(ylls), asd_deaths = sum(asd_deaths), asd_ylls = sum(asd_ylls), population = sum(population)), by = c("age_group_id", "location_id", "draw")]
global_draws_bothsex[, sex_id := 3]
global_draws_bothsex <- unique(global_draws_bothsex)

global_draws <- rbind(global_draws, global_draws_bothsex)
rm(global_draws_bothsex)

global_draws_allage <- copy(global_draws)
global_draws_allage[, `:=` (deaths = sum(deaths), ylls = sum(ylls), asd_deaths = sum(asd_deaths), asd_ylls = sum(asd_ylls), population = sum(population)), by = c("sex_id", "location_id", "draw")]
global_draws_allage[, age_group_id := 22]
global_draws_allage <- unique(global_draws_allage)

global_draws <- rbind(global_draws, global_draws_allage)
rm(global_draws_allage)

global_draws[, `:=` (deaths_rate = deaths / population, ylls_rate = ylls / population, asd_deaths_rate = asd_deaths / population, asd_ylls_rate = asd_ylls / population), ]

global_draws[, `:=` (paf_deaths = asd_deaths / deaths, paf_ylls = asd_ylls / ylls)]
global_draws[is.na(paf_deaths), paf_deaths := 0]
global_draws[is.na(paf_ylls), paf_ylls := 0]

global_draws[, `:=` (deaths_count_mean = mean(deaths), ylls_count_mean = mean(ylls), asd_deaths_count_mean = mean(asd_deaths), asd_ylls_count_mean = mean(asd_ylls),
                     deaths_count_lower = quantile(deaths, 0.025), ylls_count_lower = quantile(ylls, 0.025), asd_deaths_count_lower = quantile(asd_deaths, 0.025), asd_ylls_count_lower = quantile(asd_ylls, 0.025),
                     deaths_count_upper = quantile(deaths, 0.975), ylls_count_upper = quantile(ylls, 0.975), asd_deaths_count_upper = quantile(asd_deaths, 0.975), asd_ylls_count_upper = quantile(asd_ylls, 0.975)), by = c("age_group_id", "sex_id", "location_id")]
global_draws[, `:=` (deaths_rate_mean = mean(deaths_rate), ylls_rate_mean = mean(ylls_rate), asd_deaths_rate_mean = mean(asd_deaths_rate), asd_ylls_count_mean = mean(asd_ylls), asd_ylls_rate_mean = mean(asd_ylls_rate),
                     deaths_rate_lower = quantile(deaths_rate, 0.025), ylls_rate_lower = quantile(ylls_rate, 0.025), asd_deaths_rate_lower = quantile(asd_deaths_rate, 0.025), asd_ylls_rate_lower = quantile(asd_ylls_rate, 0.025),
                     deaths_rate_upper = quantile(deaths_rate, 0.975), ylls_rate_upper = quantile(ylls_rate, 0.975), asd_deaths_rate_upper = quantile(asd_deaths_rate, 0.975), asd_ylls_rate_upper = quantile(asd_ylls_rate, 0.975)), by = c("age_group_id", "sex_id", "location_id")]
global_draws[, `:=` (paf_deaths_mean = mean(paf_deaths), paf_deaths_lower = quantile(paf_deaths, 0.025), paf_deaths_upper = quantile(paf_deaths, 0.975),
                     paf_ylls_mean = mean(paf_ylls), paf_ylls_lower = quantile(paf_ylls, 0.025), paf_ylls_upper = quantile(paf_ylls, 0.975)), by = c("age_group_id", "sex_id", "location_id")]


# Reporting results -------------------------------------------------------

### Death estimates
## Global deaths due to suicide (regardless of ASD)
global_allsuicide_draws <- get_draws(gbd_id_type = "cause_id", gbd_id = 718, source = "codcorrect",  measure_id = c(1, 4), age_group_id = 22, metric_id = 1, location_id = 1, year_id = 2021, sex_id = 3, version_id = cod_correction_version, release_id = 9)
allsuicide_death_draws <- melt.data.table(global_allsuicide_draws[measure_id == 1, ], id.vars = names(global_allsuicide_draws)[!(names(global_allsuicide_draws) %like% "draw")], value.name = "deaths", variable.name = "draw")
allsuicide_death_draws[, .(deaths_mean = mean(deaths), deaths_lower = quantile(deaths, 0.025), deaths_upper = quantile(deaths, 0.975))]

## Total deaths due to ASD suicide (global)
unique(global_draws[location_id == 1 & sex_id == 3 & age_group_id == 22, .(asd_deaths_count_mean, asd_deaths_count_lower, asd_deaths_count_upper)])
unique(global_draws[location_id == 1 & sex_id == 1 & age_group_id == 22, .(asd_deaths_count_mean, asd_deaths_count_lower, asd_deaths_count_upper)])
unique(global_draws[location_id == 1 & sex_id == 2 & age_group_id == 22, .(asd_deaths_count_mean, asd_deaths_count_lower, asd_deaths_count_upper)])

## Death rate due to ASD suicide (global)
unique(global_draws[location_id == 1 & sex_id == 3 & age_group_id == 22, .(asd_deaths_rate_mean*1000000, asd_deaths_rate_lower*1000000, asd_deaths_rate_upper*1000000)])

## Death rate per ASD case due to suicide (global)
global_asd_prev <- get_draws(gbd_id_type = "cause_id", gbd_id = 575, source = "como",  measure_id = c(5), age_group_id = 22, metric_id = 3, location_id = unique(global_draws$location_id), year_id = 2021, sex_id = c(1:3), version_id = como_version, release_id = 9)
global_asd_prev <- melt.data.table(global_asd_prev, id.vars = names(global_asd_prev)[!(names(global_asd_prev) %like% "draw")], value.name = "prev", variable.name = "draw")
global_draws <- merge(global_asd_prev, global_draws, all =T, by = c("location_id", "sex_id", "draw", "age_group_id"))
global_draws[, death_per_asd := asd_deaths_rate / prev]
global_draws[location_id == 1 & age_group_id == 22, .(death_per_asd = mean(death_per_asd)*100000, lower = quantile(death_per_asd, 0.025)*100000, upper = quantile(death_per_asd, 0.975)*100000), by = c("sex_id")]


## Proportion of deaths due to ASD suicide (global)
unique(global_draws[location_id == 1 & sex_id == 3 & age_group_id == 22, .(paf_deaths_mean, paf_deaths_lower, paf_deaths_upper)])

## MAP
map_data <- country_overall_draws[,.(death_rate = 1000000*mean(asd_death_rate), yll_rate = 100000*mean(asd_yll_rate)), by = c("location_id")]
map_data[, mapvar := yll_rate]

bin_levels <- 10
bins_start <- quantile(map_data$mapvar, seq(0, 1, 1/bin_levels))
bins <- round(bins_start, 1)
bins[1] <- bins[1] - 1
bins[length(bins)] <- bins[length(bins)] +1
binlabels <- paste("<", round_c(bins[2],1))
for(i in 3:length(bins)){
  if(i != length(bins)){
    binlabels <- c(binlabels, paste(round_c(bins[i-1],1), "to <", round_c(bins[i],1)))
  } else {
    binlabels <- c(binlabels, paste(">=",round_c(bins[i-1],1)))
  }
}

map_title <-  paste0("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nFigure 1: Excess years of life lost due to suicide among persons on the autism spectrum, 2021")
legend_title <- "Excess years of life lost per 100 000 persons"

gbd_map(data = map_data, limits = bins, sub_nat="none", legend=TRUE, inset=F,
        labels=binlabels,
        pattern=NULL,
        col="Spectral", # Spectral or RdYlBu
        col.reverse=TRUE, na.color = "lightgrey",
        title=map_title,
        fname=paste0("/FILEPATH/ylls_map.pdf"),
        legend.title=legend_title, legend.columns = NULL, legend.cex=1.3)

### YLL estimates
## Total YLLs due to ASD suicide (global)
unique(global_draws[location_id == 1 & sex_id == 3 & age_group_id == 22, .(asd_ylls_count_mean, asd_ylls_count_lower, asd_ylls_count_upper)])
unique(global_draws[location_id == 1 & sex_id == 1 & age_group_id == 22, .(asd_ylls_count_mean, asd_ylls_count_lower, asd_ylls_count_upper)])
unique(global_draws[location_id == 1 & sex_id == 2 & age_group_id == 22, .(asd_ylls_count_mean, asd_ylls_count_lower, asd_ylls_count_upper)])

## YLL rate due to ASD suicide (global)
unique(global_draws[location_id == 1 & sex_id == 3 & age_group_id == 22, .(asd_ylls_rate_mean*100000, asd_ylls_rate_lower*100000, asd_ylls_rate_upper*100000)])
unique(global_draws[location_id == 1 & sex_id == 1 & age_group_id == 22, .(asd_ylls_rate_mean*100000, asd_ylls_rate_lower*100000, asd_ylls_rate_upper*100000)])
unique(global_draws[location_id == 1 & sex_id == 2 & age_group_id == 22, .(asd_ylls_rate_mean*100000, asd_ylls_rate_lower*100000, asd_ylls_rate_upper*100000)])

## YLL rate per ASD case due to suicide (global)
global_draws[, yll_per_asd := asd_ylls_rate / prev]
global_draws[location_id == 1 & age_group_id == 22, .(yll_per_asd = mean(yll_per_asd)*100000, lower = quantile(yll_per_asd, 0.025)*100000, upper = quantile(yll_per_asd, 0.975)*100000), by = c("sex_id")]

## New ASD DALYs and Proportion of ASD DALYs due to YLLs
asd_daly_draws <- get_draws(gbd_id_type = "cause_id", gbd_id = c(575), source = "como", measure_id = 3, metric_id = 3, location_id = unique(global_draws$location_id), year_id = 2021, sex_id = c(1, 2, 3), age_group_id =22, version_id = como_version, release_id = 9)
asd_daly_draws <- melt.data.table(asd_daly_draws, id.vars = names(asd_daly_draws)[!(names(asd_daly_draws) %like% "draw")], value.name = "asd_ylds_rate", variable.name = "draw")
global_draws <- merge(asd_daly_draws, global_draws, all = T, by = c("sex_id", "draw", "age_group_id", "location_id"))
global_draws[, asd_ylds := asd_ylds_rate * population]
global_draws[, asd_dalys := asd_ylls + asd_ylds]
global_draws[, prop_ylls := asd_ylls / asd_dalys]
global_draws[location_id == 1 & age_group_id == 22, .(asd_ylds = mean(asd_ylds), lower = quantile(asd_ylds, 0.025), upper = quantile(asd_ylds, 0.975)), by = c("sex_id")]
global_draws[location_id == 1 & age_group_id == 22, .(asd_dalys = mean(asd_dalys), lower = quantile(asd_dalys, 0.025), upper = quantile(asd_dalys, 0.975)), by = c("sex_id")]
global_draws[location_id == 1 & age_group_id == 22, .(prop_ylls = mean(prop_ylls), lower = quantile(prop_ylls, 0.025), upper = quantile(prop_ylls, 0.975)), by = c("sex_id")]

## Change in DALY ranking
print("Be sure to crosscheck the correct compare version")
all_causes <- get_outputs(topic="cause", compare_version_id = compare_version, location_id=1, year_id=c(2021), age_group_id=c(22), sex_id = c(3), measure_id=c(1, 2, 4), metric_id=c(1), cause_id='lvl4', release_id = 9)
all_causes <- all_causes[!(grepl("Total", cause_name)),] # Remove custom aggregate groups that are not distinct causes
all_causes[, old_rank:=rank(-val), by=c("sex_id", "measure_id")]
all_causes[cause_id == 575 & measure_id == 1, val :=  global_draws[location_id == 1 & age_group_id == 22 & sex_id == 3, mean(asd_deaths)]]
all_causes[cause_id == 575 & measure_id == 2, val :=  global_draws[location_id == 1 & age_group_id == 22 & sex_id == 3, mean(asd_dalys)]]
all_causes[cause_id == 575 & measure_id == 4, val :=  global_draws[location_id == 1 & age_group_id == 22 & sex_id == 3, mean(asd_ylls)]]
all_causes[, new_rank:=rank(-val), by=c("sex_id", "measure_id")]
all_causes[cause_id == 575, .(measure_id, sex_id, old_rank, new_rank)]
length(all_causes[measure_id == 1 & !is.na(val), val])
length(all_causes[measure_id == 2 & !is.na(val), val])
length(all_causes[measure_id == 4 & !is.na(val), val])

## Show the disorders these YLLs outrank
new_YLL_rankings <-all_causes[measure_id == 4 & new_rank > all_causes[cause_id == 575 & measure_id == 4, new_rank] & !is.na(val),.(cause_name, val, new_rank)]
new_YLL_rankings[order(new_rank),]

new_death_rankings <-all_causes[measure_id == 1 & new_rank > all_causes[cause_id == 575 & measure_id == 1, new_rank] & !is.na(val),.(cause_name, val, new_rank)]
new_death_rankings[order(new_rank),]

### Make Table

table_data <- global_draws[level < 2 & age_group_id == 22,]
table_data <- table_data[,.(asd_deaths_rate = mean_and_UIs(asd_deaths_rate*100000), paf_deaths = mean_and_UIs(paf_deaths*100), death_per_asd = mean_and_UIs(death_per_asd*100000),
                            asd_ylls_rate = mean_and_UIs(asd_ylls_rate*100000), paf_ylls = mean_and_UIs(paf_ylls*100), yll_per_asd = mean_and_UIs(yll_per_asd*100000), prop_ylls = mean_and_UIs(prop_ylls*100)), by = c("location_id", "sex_id")]

table_data <- merge(table_data, locations[,.(location_id, lancet_label, sort_order)], by = 'location_id')
table_data <- table_data[order(sort_order, -sex_id), .(Location = lancet_label, Sex = ifelse(sex_id == 1, "Male", ifelse(sex_id == 2, "Female", "Both")),
                                                       `Excess suicide deaths (per 100 000 persons)` = asd_deaths_rate, `Proportion of suicide deaths (%)` = paf_deaths,
                                                       `Excess suicide death rate (per 100 000 autistic persons)` = death_per_asd, `Excess YLL rate (per 100 000 persons)` = asd_ylls_rate,
                                                       `Proportion of suicide YLLs (%)` = paf_ylls, `Excess suicide YLLs (per 100 000 autistic persons)` = yll_per_asd,
                                                       `Proportion of autism spectrum DALYs (%)` = prop_ylls)]

write.csv(table_data, "/FILEPATH/table_results_by_superregion.csv", row.names = F)

## Make rate by age plot
asd_daly_draws <- get_draws(gbd_id_type = "cause_id", gbd_id = c(575), source = "como", measure_id = c(5,3), metric_id = 3, location_id = 1, year_id = 2021, sex_id = c(1, 2), version_id = como_version, release_id = 9)
asd_daly_draws <- melt.data.table(asd_daly_draws, id.vars = names(asd_daly_draws)[!(names(asd_daly_draws) %like% "draw")], value.name = "asd_ylds_rate", variable.name = "draw")
asd_daly_draws <- merge(asd_daly_draws[measure_id == 3,.(age_group_id, location_id, sex_id, draw, asd_ylds_rate)], asd_daly_draws[measure_id == 5,.(age_group_id, sex_id, draw, asd_prev = asd_ylds_rate)], by = c("age_group_id", "sex_id", "draw"))
asd_daly_draws <- merge(global_draws[,.(age_group_id, sex_id, location_id, draw, asd_deaths_rate, asd_ylls, population, asd_ylls_rate)], asd_daly_draws, by = c('age_group_id', 'location_id', 'sex_id', 'draw'))
asd_daly_draws <- merge(asd_daly_draws, age_ids, by = 'age_group_id')

asd_daly_draws[, asd_dalys_rate := asd_ylls_rate + asd_ylds_rate]
asd_daly_draws[, asd_total_rate := asd_deaths_rate + asd_prev]
asd_daly_draws[, `:=` (asd_dalys_count = asd_dalys_rate * population, asd_ylds_count = asd_ylds_rate * population, asd_total_count = asd_total_rate * population, asd_prev_count = asd_prev*population)]

under_5s <- asd_daly_draws[age_group_id %in% c(2, 3, 34, 238, 388, 389),]
under_5s[, `:=` (population = sum(population), asd_dalys_count = sum(asd_dalys_count), asd_ylds_count = sum(asd_ylds_count), asd_total_count = sum(asd_total_count), asd_prev_count = sum(asd_prev_count)), by = c('sex_id', 'draw')]
under_5s[, `:=` (asd_dalys_rate = asd_dalys_count / population, asd_ylds_rate = asd_ylds_count / population, asd_total_rate = asd_total_count / population, asd_prev = asd_prev_count / population)]

under_5s <- unique(under_5s[,.(age_group_id = 1, sex_id, draw, age_group_name = "Under 5", age_start = 0, asd_dalys_rate, asd_ylds_rate, asd_total_rate, asd_prev)])

asd_daly_draws <- rbind(asd_daly_draws[!(age_group_id %in% c(2, 3, 34, 238, 388, 389)), .(age_group_id, sex_id, draw, age_group_name, age_start, asd_dalys_rate, asd_ylds_rate, asd_total_rate, asd_prev)], under_5s)

asd_daly_draws[, `:=` (asd_dalys_rate = mean(asd_dalys_rate), asd_dalys_lower = quantile(asd_dalys_rate, 0.025), asd_dalys_upper = quantile(asd_dalys_rate, 0.975),
                       asd_ylds_rate = mean(asd_ylds_rate), asd_ylds_lower = quantile(asd_ylds_rate, 0.025), asd_ylds_upper = quantile(asd_ylds_rate, 0.975),
                       asd_total_rate = mean(asd_total_rate), asd_total_lower = quantile(asd_total_rate, 0.025), asd_total_upper = quantile(asd_total_rate, 0.975),
                       asd_prev = mean(asd_prev), asd_prev_lower = quantile(asd_prev, 0.025), asd_prev_upper = quantile(asd_prev, 0.975)), by = c("age_group_id", "sex_id")]
asd_daly_draws[, draw := NULL]

asd_daly_draws <- unique(asd_daly_draws)

epi_plot <- rbind(asd_daly_draws[, .(Prevalence = "Prevalence of autistic persons", sex = ifelse(sex_id == 1, "Male", "Female"), age_start, prev = asd_prev, lower = asd_prev_lower, upper = asd_prev_upper)],
                                  asd_daly_draws[, .(Prevalence = "Suicide rate of autistic persons", sex = ifelse(sex_id == 1, "Male", "Female"), age_start, prev = asd_total_rate, lower = asd_total_lower, upper = asd_total_upper)])

burden_plot <- rbind(asd_daly_draws[, .(`Health burden` = "YLDs", sex = ifelse(sex_id == 1, "Male", "Female"), age_start, prev = asd_ylds_rate*100000, lower = asd_ylds_lower*100000, upper = asd_ylds_upper*100000)],
                    asd_daly_draws[, .(`Health burden` = "DALYs", sex = ifelse(sex_id == 1, "Male", "Female"), age_start, prev = asd_dalys_rate*100000, lower = asd_dalys_lower*100000, upper = asd_dalys_upper*100000)])

epi_plot$sex <- factor(epi_plot$sex, levels=c("Female", "Male"))
burden_plot$sex <- factor(epi_plot$sex, levels=c("Female", "Male"))

plot_age <- ggplot(data=burden_plot, aes(x=as.numeric(age_start)+2, y=prev, color = `Health burden`, fill = `Health burden`)) +
  geom_ribbon(data= burden_plot, aes(x=as.numeric(age_start)+2, ymin=lower, ymax=upper, color = `Health burden`, fill = `Health burden`), alpha=.3) +
  geom_line(size=1) +
  scale_color_manual(values = c("YLDs" = "#000bab", "DALYs" = "#b00000"))+
  scale_fill_manual(values = c("YLDs" = "#000bab", "DALYs" = "#b00000"))+
  facet_wrap(~sex)+
  ylab("YLDs/DALYs per 100 000 people") +
  xlab("Age (years)") +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0), breaks = c(0, 50, 100, 150, 200, 250, 300, 350), limits = c(0, 350)) +
  labs(caption = "Figure 2: Global health burden of persons on the autism spectrum by age and sex, 2021") +
  theme(axis.line=element_line(colour="black"), plot.caption = element_text(hjust=0.25, size=rel(0.9)))
plot_age

ggsave(plot_age, filename="/FILEPATH/figure_2.pdf", width = 9, height = 6)


## Make pyramid plot
asd_daly_draws <- get_draws(gbd_id_type = "cause_id", gbd_id = c(575), source = "como", measure_id = c(3), metric_id = 3, location_id = 1, year_id = 2021, sex_id = c(1, 2), version_id = como_version, release_id = 9)
asd_daly_draws <- melt.data.table(asd_daly_draws, id.vars = names(asd_daly_draws)[!(names(asd_daly_draws) %like% "draw")], value.name = "asd_ylds_rate", variable.name = "draw")
asd_daly_draws <- merge(global_draws[,.(age_group_id, sex_id, location_id, draw, asd_ylls, population, asd_ylls_rate)], asd_daly_draws[,.(age_group_id, location_id, sex_id, year_id, draw, asd_ylds_rate)], by = c('age_group_id', 'location_id', 'sex_id', 'draw'))
asd_daly_draws <- merge(asd_daly_draws, age_ids, by = 'age_group_id')

asd_daly_draws[, asd_ylds := asd_ylds_rate * population]
asd_daly_draws[, asd_dalys := asd_ylls + asd_ylds]
asd_daly_draws[, `:=` (asd_ylls_mean = mean(asd_ylls), asd_ylds_mean = mean(asd_ylds)), by = c('sex_id', 'age_group_id')]

under_5s <- asd_daly_draws[age_group_id %in% c(2, 3, 34, 238, 388, 389),]
under_5s[, `:=` (asd_ylds = sum(asd_ylds), asd_ylls = sum(asd_ylls)), by = c('sex_id', 'draw')]
under_5s <- unique(under_5s[,.(age_group_id = 1, sex_id, draw, age_group_name = "Under 5", age_start = 0, asd_ylds, asd_ylls)])
under_5s[, `:=` (asd_ylls_mean = mean(asd_ylls), asd_ylds_mean = mean(asd_ylds)), by = c('sex_id', 'age_group_id')]
under_5s[, (c('asd_ylls', 'asd_ylds')) := NULL]

asd_daly_draws <- rbind(asd_daly_draws[!(age_group_id %in% c(2, 3, 34, 238, 388, 389)), .(age_group_id, sex_id, draw, age_group_name, age_start, asd_ylds_mean, asd_ylls_mean)], under_5s)

pyramid_data <- unique(rbind(asd_daly_draws[,.(Cause = "YLDs", sex_id, age_group_id,  Age = age_group_name, age_start, `Disability-adjusted life-years (thousands)`= asd_ylds_mean/1000)],
                      asd_daly_draws[,.(Cause = "YLLs", sex_id, age_group_id,  Age = age_group_name, age_start, `Disability-adjusted life-years (thousands)` = asd_ylls_mean/1000)]))
pyramid_data[sex_id == 2, `:=` (`Disability-adjusted life-years (thousands)` = -`Disability-adjusted life-years (thousands)`)]
pyramid_data$Age <- factor(pyramid_data$Age, levels=c("Under 5", age_ids$age_group_name))
pyramid_plot <- ggplot(data=pyramid_data, aes(x=Age, y=`Disability-adjusted life-years (thousands)`, fill=Cause)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("YLLs" = "#e62e00", "YLDs" = "#0035a8"))+
  scale_y_continuous(breaks = seq(-400, 800, 50), labels = abs(seq(-400, 800, 50)), limits = c(-400, 800))+ 
  theme_minimal() +
  geom_hline(yintercept = 0) +
  geom_text(x=18, y=-225, label="Females", size = 5) +
  geom_text(x=18, y=425, label="Males", size = 5) +
  labs(caption = "Figure 3: Global DALYs attributable to the autism spectrum in 2021 by YLDs, YLLs, age, and sex.") +
  theme(plot.caption = element_text(hjust=0.6, size=rel(0.8)))+
  coord_flip()
pyramid_plot
ggsave(pyramid_plot, filename="/FILEPATH/figure_3.pdf", width = 8, height = 6)


## Appendix table
country_overall_draws[, `:=` (paf_deaths = asd_death_count / death_count, paf_ylls = asd_yll_count / yll_count)]
country_asd_prev <- get_draws(gbd_id_type = "cause_id", gbd_id = 575, source = "como",  measure_id = c(5), age_group_id = 22, metric_id = 3, location_id = unique(country_overall_draws$location_id), year_id = 2021, sex_id = c(1:3), version_id = como_version, release_id = 9)
country_asd_prev <- melt.data.table(country_asd_prev, id.vars = names(country_asd_prev)[!(names(country_asd_prev) %like% "draw")], value.name = "prev", variable.name = "draw")
country_overall_draws <- merge(country_asd_prev, country_overall_draws, all =T, by = c("location_id", "sex_id", "draw"))
country_overall_draws[, death_per_asd := asd_death_rate / prev]

country_asd_daly_draws <- get_draws(gbd_id_type = "cause_id", gbd_id = c(575), source = "como", measure_id = 3, metric_id = 3, location_id = unique(country_overall_draws$location_id), year_id = 2021, sex_id = c(1, 2, 3), age_group_id =22, version_id = como_version, release_id = 9)
country_asd_daly_draws <- melt.data.table(country_asd_daly_draws, id.vars = names(country_asd_daly_draws)[!(names(country_asd_daly_draws) %like% "draw")], value.name = "asd_ylds_rate", variable.name = "draw")
country_overall_draws <- merge(country_asd_daly_draws, country_overall_draws, all = T, by = c("sex_id", "draw", "location_id"))
country_overall_draws[, asd_ylds := asd_ylds_rate * pop]
country_overall_draws[, asd_dalys := asd_yll_count + asd_ylds]
country_overall_draws[, prop_ylls := asd_yll_count / asd_dalys]
country_overall_draws[, yll_per_asd := asd_yll_rate / prev]

country_table <- country_overall_draws[,.(asd_deaths_rate = mean_and_UIs(asd_death_rate*100000), paf_deaths = mean_and_UIs(paf_deaths*100), death_per_asd = mean_and_UIs(death_per_asd*100000),
                                          asd_ylls_rate = mean_and_UIs(asd_yll_rate*100000), paf_ylls = mean_and_UIs(paf_ylls*100), yll_per_asd = mean_and_UIs(yll_per_asd*100000), prop_ylls = mean_and_UIs(prop_ylls*100)), by = c("location_id", "sex_id")]
table_data <- global_draws[level < 3 & age_group_id == 22,]
table_data <- table_data[,.(asd_deaths_rate = mean_and_UIs(asd_deaths_rate*100000), paf_deaths = mean_and_UIs(paf_deaths*100), death_per_asd = mean_and_UIs(death_per_asd*100000),
                            asd_ylls_rate = mean_and_UIs(asd_ylls_rate*100000), paf_ylls = mean_and_UIs(paf_ylls*100), yll_per_asd = mean_and_UIs(yll_per_asd*100000), prop_ylls = mean_and_UIs(prop_ylls*100)), by = c("location_id", "sex_id")]

country_table <- rbind(country_table, table_data)

country_table <- merge(country_table, locations[,.(location_id, lancet_label, sort_order, level)], by = 'location_id')
country_table <- country_table[order(sort_order, -sex_id), .(level, Location = lancet_label, Sex = ifelse(sex_id == 1, "Male", ifelse(sex_id == 2, "Female", "Both")),
                                                       `Excess suicide deaths (per 100 000 persons)` = asd_deaths_rate, `Proportion of suicide deaths (%)` = paf_deaths,
                                                       `Excess suicide death rate (per 100 000 autistic persons)` = death_per_asd, `Excess YLL rate (per 100 000 persons)` = asd_ylls_rate,
                                                       `Proportion of suicide YLLs (%)` = paf_ylls, `Excess suicide YLLs (per 100 000 autistic persons)` = yll_per_asd,
                                                       `Proportion of autism spectrum DALYs (%)` = prop_ylls)]
country_table <- country_table[!(Location %in% c("North Africa and Middle East", "South Asia") & level == 2),]
country_table[, row_id := seq_len(.N)]

wb <- createWorkbook()
addWorksheet(wb, "Sheet 1")
writeData(wb, "Sheet 1", country_table,
          colNames = T, rowNames = F, startCol = "A",
          startRow = 1, borders = "all",
)

for(l in unique(country_table$Location)){
  rows <- country_table[Location == l, row_id+1]
  mergeCells(wb, "Sheet 1", cols = 2, rows = rows)
}

saveWorkbook(wb, "/FILEPATH/table_results_by_country.xlsx", overwrite = TRUE)



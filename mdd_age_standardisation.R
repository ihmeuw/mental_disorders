rm(list=ls())

library(data.table)
library(ggplot2)
source("/FILEPATH/get_draws.R")
source("/FILEPATH/interpolate.R")
source("/FILEPATH/get_age_metadata.R")
source("/FILEPATH/get_population.R")
source("/FILEPATH/get_location_metadata.R")
source("/FILEPATH/get_crosswalk_version.R")
source("/FILEPATH/get_outputs.R")
source("/FILEPATH/get_pct_change.R")

source("/FILEPATH/gbd2021_map.R")

years <- c(2021)
release <- 9 
como_version <- 1471
pop_version <- 359 
compare_version <- 8016

# MDD ranking for introduction --------------------------------------------
data <- get_outputs(topic="cause", compare_version_id = compare_version, location_id=1, year_id=2021, age_group_id=22, sex_id = c(3), measure_id=c(2, 3), metric_id=c(1), cause_id="lvl4", release_id = 9)
data <- data[!grepl("Total ", cause_name), ]
data[, rank := rank(-val), by = c('measure_id')]
data[cause_name == "Major depressive disorder",]

daly_change_1990_2019 <- get_pct_change(gbd_id_type = 'cause_id', gbd_id = 568, year_start_id = 1990, year_end_id = 2019, source = 'dalynator', change_type = 'pct_change_rate', release_id = 9, measure_id = 2, location_id = 1, sex_id = 3, age_group_id = 27, version_id = 78)
daly_change_2019_2020 <- get_pct_change(gbd_id_type = 'cause_id', gbd_id = 568, year_start_id = 2019, year_end_id = 2020, source = 'dalynator', change_type = 'pct_change_rate', release_id = 9, measure_id = 2, location_id = 1, sex_id = 3, age_group_id = 27, version_id = 78)
daly_change_2019_2021 <- get_pct_change(gbd_id_type = 'cause_id', gbd_id = 568, year_start_id = 2019, year_end_id = 2021, source = 'dalynator', change_type = 'pct_change_rate', release_id = 9, measure_id = 2, location_id = 1, sex_id = 3, age_group_id = 27, version_id = 78)

# Generate age-standard for MDD cases globally ----------------------------
ages <- get_age_metadata(release_id = release)
locations <- get_location_metadata(location_set_id = 9, release_id = release)
## MAT DisMod-mr 2.1 model had to be rerun within the GBD 2023 environment due to timing of peer-review feedback, so adjustments made to make location hierarchy comparable
locations_gbd2023 <- get_location_metadata(location_set_id = 9, release_id = 16)
locations_estimated <- locations[is_estimate == 1, location_id]
locations_estimated_gbd2023 <- locations_gbd2023[is_estimate == 1, location_id]
missing_locations <- locations_estimated[!(locations_estimated %in% locations_estimated_gbd2023)]
locations[location_id %in% missing_locations,location_name]
locations_estimated <- c(locations_estimated[!(locations_estimated %in% missing_locations)], 4749, 93, 179)

mdd_prev <- get_draws(version_id = como_version, year_id = years, gbd_id_type = "cause_id", gbd_id = 568, source = "como", release_id = release, metric_id = 3, sex_id = c(1, 2, 3), measure_id = 5, age_group_id = c(22, ages$age_group_id), location_id = 1)

mdd_prev <- rbind(mdd_prev_2021, mdd_prev_2000)
rm(mdd_prev_2021, mdd_prev_2000)

population <- get_population(run_id = pop_version, release_id = release, sex_id = c(1, 2), location_id = unique(mdd_prev$location_id), year_id = c(2000, 2021), age_group_id = c(ages$age_group_id))

mdd_prev <- melt.data.table(mdd_prev, id.vars = names(mdd_prev)[!(names(mdd_prev) %like% "draw")], value.name="prev", variable.name="draw")
mdd_prev <- merge(mdd_prev, population, all.x = T, by = c("location_id", "year_id", "sex_id", "age_group_id"))
mdd_prev[, cases := prev * population]


# Pull estimates of MDD MAT coverage --------------------------------------
mat_coverage <- get_draws(year_id = years, gbd_id_type = "modelable_entity_id", gbd_id = 27266, source = "epi", release_id = 9, metric_id = 3, sex_id = c(1, 2), measure_id = 18, age_group_id = ages$age_group_id)

mat_coverage <- melt.data.table(mat_coverage, id.vars = names(mat_coverage)[!(names(mat_coverage) %like% "draw")], value.name="mat_prop", variable.name="draw")

locations_estimated[!(locations_estimated %in% mat_coverage$location_id)] # All locations present
unique(mat_coverage$location_id[!(mat_coverage$location_id %in% locations_estimated)])

mat_coverage[, mat_cases := mat_prop * cases]

locations[, country_id:=unlist(strsplit(path_to_top_parent , ","))[4], by = 'path_to_top_parent']
locations[, country_id:=unlist(strsplit(path_to_top_parent , ","))[4], by = 'path_to_top_parent']
locations[country_id == 95 & level == 4, country_id := location_id]
locations[country_id == 95 & level > 4, country_id := 4749]

mat_coverage <- merge(mat_coverage, locations[, .(location_id, location_name, super_region_id, super_region_name, region_id, region_name, country_id, is_estimate)], by = 'location_id', all.x = T)

mat_coverage <- merge(mat_coverage[sex_id == 1,.(location_id, location_name, super_region_id, super_region_name, region_id, region_name, country_id, year_id, age_group_id, draw, population_m = population, cases_m = cases, mat_cases_m = mat_cases)],
                      mat_coverage[sex_id == 2,.(location_id, year_id, age_group_id, draw, population_f = population, cases_f = cases, mat_cases_f = mat_cases)], by = c('location_id', 'year_id', 'draw', 'age_group_id'))

mat_coverage[, `:=` (cases = cases_m + cases_f, population = population_m + population_f, mat_cases = mat_cases_m + mat_cases_f)]

mat_countries <- mat_coverage[country_id == location_id, ]

for(i in unique(mat_coverage[country_id != location_id, country_id])){
  temp_data <- mat_coverage[country_id == i, ]
  temp_data[, `:=` (population_m = sum(population_m), cases_m = sum(cases_m), mat_cases_m = sum(mat_cases_m),
                    population_f = sum(population_f), cases_f = sum(cases_f), mat_cases_f = sum(mat_cases_f),
                    population = sum(population), cases = sum(cases), mat_cases = sum(mat_cases)), by = c('year_id', 'age_group_id', 'draw')]
  c_name <- locations[location_id == i, location_name]
  temp_data[, `:=` (location_id = i, location_name = c_name)]
  temp_data <- unique(temp_data)
  mat_countries <- rbind(mat_countries, temp_data)
}

mat_regions <- data.table()

for(i in unique(mat_coverage[, region_id])){
  temp_data <- mat_countries[region_id == i, ]
  temp_data[, `:=` (population_m = sum(population_m), cases_m = sum(cases_m), mat_cases_m = sum(mat_cases_m),
                    population_f = sum(population_f), cases_f = sum(cases_f), mat_cases_f = sum(mat_cases_f),
                    population = sum(population), cases = sum(cases), mat_cases = sum(mat_cases)), by = c('year_id', 'age_group_id', 'draw')]
  temp_data[, `:=` (location_id = i, location_name = region_name, country_id = NA)]
  temp_data <- unique(temp_data)
  mat_regions <- rbind(mat_regions, temp_data)
}

mat_super_regions <- data.table()

for(i in unique(mat_coverage[, super_region_id])){
  temp_data <- mat_regions[super_region_id == i, ]
  temp_data[, `:=` (population_m = sum(population_m), cases_m = sum(cases_m), mat_cases_m = sum(mat_cases_m),
                    population_f = sum(population_f), cases_f = sum(cases_f), mat_cases_f = sum(mat_cases_f),
                    population = sum(population), cases = sum(cases), mat_cases = sum(mat_cases)), by = c('year_id', 'age_group_id', 'draw')]
  temp_data[, `:=` (location_id = i, location_name = super_region_name, region_name = NA, region_id = NA, country_id = NA)]
  temp_data <- unique(temp_data)
  mat_super_regions <- rbind(mat_super_regions, temp_data)
}


mat_global <- copy(mat_super_regions)
mat_global[, `:=` (population_m = sum(population_m), cases_m = sum(cases_m), mat_cases_m = sum(mat_cases_m),
                   population_f = sum(population_f), cases_f = sum(cases_f), mat_cases_f = sum(mat_cases_f),
                   population = sum(population), cases = sum(cases), mat_cases = sum(mat_cases)), by = c('year_id', 'age_group_id', 'draw')]
mat_global[, `:=` (location_id = 1, location_name = "Global", super_region_name = NA, super_region_id = NA, region_name = NA, region_id = NA, country_id = NA)]
mat_global <- unique(mat_global)

# Global all-age ----------------------------------------------------------

mat_global_all_age <- copy(mat_global)
mat_global_all_age[, `:=` (population_m = sum(population_m), cases_m = sum(cases_m), mat_cases_m = sum(mat_cases_m),
                           population_f = sum(population_f), cases_f = sum(cases_f), mat_cases_f = sum(mat_cases_f),
                           population = sum(population), cases = sum(cases), mat_cases = sum(mat_cases)), by = c('year_id', 'draw')]
mat_global_all_age[, `:=` (age_group_id = 22)]
mat_global_all_age <- unique(mat_global_all_age)

mat_global_all_age[, `:=` (mat_prop = mat_cases / cases, mat_prop_m = mat_cases_m / cases_m, mat_prop_f = mat_cases_f / cases_f)]

mat_global_all_age[, .(mean = mean(mat_prop), lower = quantile(mat_prop, 0.025), upper = quantile(mat_prop, 0.975)), by = 'year_id']

# Female
mat_global_all_age[, .(mean = mean(mat_prop_f), lower = quantile(mat_prop_f, 0.025), upper = quantile(mat_prop_f, 0.975)), by = 'year_id']

# Male
mat_global_all_age[, .(mean = mean(mat_prop_m), lower = quantile(mat_prop_m, 0.025), upper = quantile(mat_prop_m, 0.975)), by = 'year_id']

# Super region all-age ----------------------------------------------------
mat_super_regions_all_age <- copy(mat_super_regions)
mat_super_regions_all_age[, `:=` (population_m = sum(population_m), cases_m = sum(cases_m), mat_cases_m = sum(mat_cases_m),
                                  population_f = sum(population_f), cases_f = sum(cases_f), mat_cases_f = sum(mat_cases_f),
                                  population = sum(population), cases = sum(cases), mat_cases = sum(mat_cases)), by = c('year_id', 'draw', 'location_id')]
mat_super_regions_all_age[, `:=` (age_group_id = 22)]
mat_super_regions_all_age <- unique(mat_super_regions_all_age)

mat_super_regions_all_age[, `:=` (mat_prop = mat_cases / cases, mat_prop_m = mat_cases_m / cases_m, mat_prop_f = mat_cases_f / cases_f)]

mat_super_regions_all_age[year_id == 2021, .(mean = mean(mat_prop), lower = quantile(mat_prop, 0.025), upper = quantile(mat_prop, 0.975)), by = c('super_region_name', 'year_id')]

# Region all-age -----------------------------------------------------------
mat_regions_all_age <- copy(mat_regions)
mat_regions_all_age[, `:=` (population_m = sum(population_m), cases_m = sum(cases_m), mat_cases_m = sum(mat_cases_m),
                            population_f = sum(population_f), cases_f = sum(cases_f), mat_cases_f = sum(mat_cases_f),
                            population = sum(population), cases = sum(cases), mat_cases = sum(mat_cases)), by = c('year_id', 'draw', 'location_id')]
mat_regions_all_age[, `:=` (age_group_id = 22)]
mat_regions_all_age <- unique(mat_regions_all_age)

mat_regions_all_age[, `:=` (mat_prop = mat_cases / cases, mat_prop_m = mat_cases_m / cases_m, mat_prop_f = mat_cases_f / cases_f)]

mat_regions_all_age[year_id == 2021, .(mean = mean(mat_prop), lower = quantile(mat_prop, 0.025), upper = quantile(mat_prop, 0.975)), by = c('region_name', 'year_id')]

# Countries all-age -------------------------------------------------------
mat_countries_all_age <- copy(mat_countries)
mat_countries_all_age[, `:=` (population_m = sum(population_m), cases_m = sum(cases_m), mat_cases_m = sum(mat_cases_m),
                              population_f = sum(population_f), cases_f = sum(cases_f), mat_cases_f = sum(mat_cases_f),
                              population = sum(population), cases = sum(cases), mat_cases = sum(mat_cases)), by = c('year_id', 'draw', 'location_id')]
mat_countries_all_age[, `:=` (age_group_id = 22)]
mat_countries_all_age <- unique(mat_countries_all_age)

mat_countries_all_age[, `:=` (mat_prop = mat_cases / cases, mat_prop_m = mat_cases_m / cases_m, mat_prop_f = mat_cases_f / cases_f)]

mat_countries_all_age[year_id == 2021, .(mean = mean(mat_prop), lower = quantile(mat_prop, 0.025), upper = quantile(mat_prop, 0.975)), by = c('location_name', 'year_id')]

# Generate table data and csv -----------------------------------------------------

round_c <- function(x,y){gsub("\\.", "·", as.character(sprintf(paste0("%.",y, "f"), x)))}

###### 2021

## global

table_data_global <- unique(mat_global_all_age[year_id == 2021, .(year_id, location_id, mat_prop = mean(mat_prop), lower = quantile(mat_prop, 0.025), upper = quantile(mat_prop, 0.975))])
table_data_global_males <- unique(mat_global_all_age[year_id == 2021, .(year_id, location_id, mat_prop = mean(mat_prop_m), lower = quantile(mat_prop_m, 0.025), upper = quantile(mat_prop_m, 0.975))])
table_data_global_females <- unique(mat_global_all_age[year_id == 2021, .(year_id, location_id, mat_prop = mean(mat_prop_f), lower = quantile(mat_prop_f, 0.025), upper = quantile(mat_prop_f, 0.975))])

table_data_global <- table_data_global[, .(location_id, Both = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]
table_data_global_males <- table_data_global_males[, .(location_id, Males = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]
table_data_global_females <- table_data_global_females[, .(location_id, Females = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]

table_data_global <- merge(table_data_global, table_data_global_females, by = 'location_id')
table_data_global <- merge(table_data_global, table_data_global_males, by = 'location_id')

## super region

table_data_sr <- mat_super_regions_all_age[year_id == 2021, .(mat_prop = mean(mat_prop), lower = quantile(mat_prop, 0.025), upper = quantile(mat_prop, 0.975)), by = c('location_id')]
table_data_sr_males <- mat_super_regions_all_age[year_id == 2021, .(location_id, mat_prop = mean(mat_prop_m), lower = quantile(mat_prop_m, 0.025), upper = quantile(mat_prop_m, 0.975)), by = c('location_id')]
table_data_sr_females <- mat_super_regions_all_age[year_id == 2021, .(location_id, mat_prop = mean(mat_prop_f), lower = quantile(mat_prop_f, 0.025), upper = quantile(mat_prop_f, 0.975)), by = c('location_id')]

table_data_sr <- table_data_sr[, .(location_id, Both = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]
table_data_sr_males <- table_data_sr_males[, .(location_id, Males = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]
table_data_sr_females <- table_data_sr_females[, .(location_id, Females = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]

table_data_sr <- merge(table_data_sr, table_data_sr_females, by = 'location_id')
table_data_sr <- merge(table_data_sr, table_data_sr_males, by = 'location_id')

## region

table_data_region <- mat_regions_all_age[year_id == 2021, .(mat_prop = mean(mat_prop), lower = quantile(mat_prop, 0.025), upper = quantile(mat_prop, 0.975)), by = c('location_id')]
table_data_region_males <- mat_regions_all_age[year_id == 2021, .(location_id, mat_prop = mean(mat_prop_m), lower = quantile(mat_prop_m, 0.025), upper = quantile(mat_prop_m, 0.975)), by = c('location_id')]
table_data_region_females <- mat_regions_all_age[year_id == 2021, .(location_id, mat_prop = mean(mat_prop_f), lower = quantile(mat_prop_f, 0.025), upper = quantile(mat_prop_f, 0.975)), by = c('location_id')]

table_data_region <- table_data_region[, .(location_id, Both = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]
table_data_region_males <- table_data_region_males[, .(location_id, Males = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]
table_data_region_females <- table_data_region_females[, .(location_id, Females = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]

table_data_region <- merge(table_data_region, table_data_region_females, by = 'location_id')
table_data_region <- merge(table_data_region, table_data_region_males, by = 'location_id')

## country

table_data_country <- mat_countries_all_age[year_id == 2021, .(mat_prop = mean(mat_prop), lower = quantile(mat_prop, 0.025), upper = quantile(mat_prop, 0.975)), by = c('location_id')]
table_data_country_males <- mat_countries_all_age[year_id == 2021, .(location_id, mat_prop = mean(mat_prop_m), lower = quantile(mat_prop_m, 0.025), upper = quantile(mat_prop_m, 0.975)), by = c('location_id')]
table_data_country_females <- mat_countries_all_age[year_id == 2021, .(location_id, mat_prop = mean(mat_prop_f), lower = quantile(mat_prop_f, 0.025), upper = quantile(mat_prop_f, 0.975)), by = c('location_id')]

table_data_country <- table_data_country[, .(location_id, Both = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]
table_data_country_males <- table_data_country_males[, .(location_id, Males = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]
table_data_country_females <- table_data_country_females[, .(location_id, Females = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]

table_data_country <- merge(table_data_country, table_data_country_females, by = 'location_id')
table_data_country <- merge(table_data_country, table_data_country_males, by = 'location_id')

###### 2000

## global

table_data_global_2000 <- unique(mat_global_all_age[year_id == 2000, .(year_id, location_id, mat_prop = mean(mat_prop), lower = quantile(mat_prop, 0.025), upper = quantile(mat_prop, 0.975))])
table_data_global_males_2000 <- unique(mat_global_all_age[year_id == 2000, .(year_id, location_id, mat_prop = mean(mat_prop_m), lower = quantile(mat_prop_m, 0.025), upper = quantile(mat_prop_m, 0.975))])
table_data_global_females_2000 <- unique(mat_global_all_age[year_id == 2000, .(year_id, location_id, mat_prop = mean(mat_prop_f), lower = quantile(mat_prop_f, 0.025), upper = quantile(mat_prop_f, 0.975))])

table_data_global_2000 <- table_data_global_2000[, .(location_id, Both = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]
table_data_global_males_2000 <- table_data_global_males_2000[, .(location_id, Males = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]
table_data_global_females_2000 <- table_data_global_females_2000[, .(location_id, Females = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]

table_data_global_2000 <- merge(table_data_global_2000, table_data_global_females_2000, by = 'location_id')
table_data_global_2000 <- merge(table_data_global_2000, table_data_global_males_2000, by = 'location_id')

## super region

table_data_sr_2000 <- mat_super_regions_all_age[year_id == 2000, .(mat_prop = mean(mat_prop), lower = quantile(mat_prop, 0.025), upper = quantile(mat_prop, 0.975)), by = c('location_id')]
table_data_sr_males_2000 <- mat_super_regions_all_age[year_id == 2000, .(location_id, mat_prop = mean(mat_prop_m), lower = quantile(mat_prop_m, 0.025), upper = quantile(mat_prop_m, 0.975)), by = c('location_id')]
table_data_sr_females_2000 <- mat_super_regions_all_age[year_id == 2000, .(location_id, mat_prop = mean(mat_prop_f), lower = quantile(mat_prop_f, 0.025), upper = quantile(mat_prop_f, 0.975)), by = c('location_id')]

table_data_sr_2000 <- table_data_sr_2000[, .(location_id, Both = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]
table_data_sr_males_2000 <- table_data_sr_males_2000[, .(location_id, Males = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]
table_data_sr_females_2000 <- table_data_sr_females_2000[, .(location_id, Females = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]

table_data_sr_2000 <- merge(table_data_sr_2000, table_data_sr_females_2000, by = 'location_id')
table_data_sr_2000 <- merge(table_data_sr_2000, table_data_sr_males_2000, by = 'location_id')

## region

table_data_region_2000 <- mat_regions_all_age[year_id == 2000, .(mat_prop = mean(mat_prop), lower = quantile(mat_prop, 0.025), upper = quantile(mat_prop, 0.975)), by = c('location_id')]
table_data_region_males_2000 <- mat_regions_all_age[year_id == 2000, .(location_id, mat_prop = mean(mat_prop_m), lower = quantile(mat_prop_m, 0.025), upper = quantile(mat_prop_m, 0.975)), by = c('location_id')]
table_data_region_females_2000 <- mat_regions_all_age[year_id == 2000, .(location_id, mat_prop = mean(mat_prop_f), lower = quantile(mat_prop_f, 0.025), upper = quantile(mat_prop_f, 0.975)), by = c('location_id')]

table_data_region_2000 <- table_data_region_2000[, .(location_id, Both = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]
table_data_region_males_2000 <- table_data_region_males_2000[, .(location_id, Males = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]
table_data_region_females_2000 <- table_data_region_females_2000[, .(location_id, Females = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]

table_data_region_2000 <- merge(table_data_region_2000, table_data_region_females_2000, by = 'location_id')
table_data_region_2000 <- merge(table_data_region_2000, table_data_region_males_2000, by = 'location_id')

## country

table_data_country_2000 <- mat_countries_all_age[year_id == 2000, .(mat_prop = mean(mat_prop), lower = quantile(mat_prop, 0.025), upper = quantile(mat_prop, 0.975)), by = c('location_id')]
table_data_country_males_2000 <- mat_countries_all_age[year_id == 2000, .(location_id, mat_prop = mean(mat_prop_m), lower = quantile(mat_prop_m, 0.025), upper = quantile(mat_prop_m, 0.975)), by = c('location_id')]
table_data_country_females_2000 <- mat_countries_all_age[year_id == 2000, .(location_id, mat_prop = mean(mat_prop_f), lower = quantile(mat_prop_f, 0.025), upper = quantile(mat_prop_f, 0.975)), by = c('location_id')]

table_data_country_2000 <- table_data_country_2000[, .(location_id, Both = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]
table_data_country_males_2000 <- table_data_country_males_2000[, .(location_id, Males = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]
table_data_country_females_2000 <- table_data_country_females_2000[, .(location_id, Females = paste0(round_c(mat_prop*100, 1), " (", round_c(lower*100, 1), "–", round_c(upper*100, 1), ")"))]

table_data_country_2000 <- merge(table_data_country_2000, table_data_country_females_2000, by = 'location_id')
table_data_country_2000 <- merge(table_data_country_2000, table_data_country_males_2000, by = 'location_id')

## combine estimates

table_data <- rbind(table_data_global, table_data_sr, table_data_region, table_data_country)

table_data_appendix <- merge(table_data, locations[level < 4,.(location_id, location_name, sort_order)], by = c("location_id"))
table_data_appendix <- table_data_appendix[order(sort_order),.(location_name, location_id, Both, Females, Males)]

table_data_2000 <- rbind(table_data_global_2000, table_data_sr_2000, table_data_region_2000, table_data_country_2000)
table_data_appendix_2000 <- merge(table_data_2000, locations[level < 4,.(location_id, location_name, sort_order)], by = c("location_id"))
table_data_appendix_2000 <- table_data_appendix_2000[order(sort_order),.(location_name, location_id, Both, Females, Males)]


table_data_paper <- merge(table_data, locations[level < 3,.(location_id, location_name, sort_order)], by = c("location_id"))
table_data_paper <- table_data_paper[order(sort_order),.(location_name, location_id, Both, Females, Males)]

write.xlsx(table_data_appendix, "/FILEPATH/mat_by_location_appendix_2021.xlsx")
write.xlsx(table_data_appendix_2000, "/FILEPATH/mat_by_location_appendix_2000.xlsx")
write.xlsx(table_data_paper, "/FILEPATH/mat_by_location_paper.xlsx")

# Create table with estimates over time -----------------------------------

table_bytime <- rbind(mat_global_all_age, mat_super_regions_all_age)

table_bytime <- merge(table_bytime[year_id == 2000,.(location_id, draw, mdd_cases_2000 = cases/1000000, mdd_mat_2000 = mat_cases/1000000, mat_prop_2000 = mat_prop)],
                      table_bytime[year_id == 2021,.(location_id, draw, mdd_cases_2021 = cases/1000000, mdd_mat_2021 = mat_cases/1000000, mat_prop_2021 = mat_prop)], by = c('location_id', "draw"))

table_bytime[, mat_prop_change := (mat_prop_2021 - mat_prop_2000)/mat_prop_2000]

table_bytime[, `:=` (mdd_cases_2000_c = paste0(round_c(mean(mdd_cases_2000), 1), " (", round_c(quantile(mdd_cases_2000, 0.025), 1), "–", round_c(quantile(mdd_cases_2000, 0.975), 1), ")")), by = 'location_id']
table_bytime[, `:=` (mdd_mat_2000_c = paste0(round_c(mean(mdd_mat_2000), 1), " (", round_c(quantile(mdd_mat_2000, 0.025), 1), "–", round_c(quantile(mdd_mat_2000, 0.975), 1), ")")), by = 'location_id']
table_bytime[, `:=` (mat_prop_2000_c = paste0(round_c(mean(mat_prop_2000)*100, 1), " (", round_c(quantile(mat_prop_2000, 0.025)*100, 1), "–", round_c(quantile(mat_prop_2000, 0.975)*100, 1), ")")), by = 'location_id']
table_bytime[, `:=` (mdd_cases_2021_c = paste0(round_c(mean(mdd_cases_2021), 1), " (", round_c(quantile(mdd_cases_2021, 0.025), 1), "–", round_c(quantile(mdd_cases_2021, 0.975), 1), ")")), by = 'location_id']
table_bytime[, `:=` (mdd_mat_2021_c = paste0(round_c(mean(mdd_mat_2021), 1), " (", round_c(quantile(mdd_mat_2021, 0.025), 1), "–", round_c(quantile(mdd_mat_2021, 0.975), 1), ")")), by = 'location_id']
table_bytime[, `:=` (mat_prop_2021_c = paste0(round_c(mean(mat_prop_2021)*100, 1), " (", round_c(quantile(mat_prop_2021, 0.025)*100, 1), "–", round_c(quantile(mat_prop_2021, 0.975)*100, 1), ")")), by = 'location_id']
table_bytime[, `:=` (mat_prop_change_c = paste0(round_c(mean(mat_prop_change)*100, 1), " (", round_c(quantile(mat_prop_change, 0.025)*100, 1), "–", round_c(quantile(mat_prop_change, 0.975)*100, 1), ")")), by = 'location_id']

table_bytime <- unique(table_bytime[,.(location_id, mdd_cases_2000_c, mdd_mat_2000_c, mat_prop_2000_c, mdd_cases_2021_c, mdd_mat_2021_c, mat_prop_2021_c, mat_prop_change_c)])

table_bytime <- merge(table_bytime, locations[,.(location_id, location_name, sort_order)], by = c("location_id"))
table_bytime <- table_bytime[order(sort_order), .(location_name, `MDD cases in 2000` = mdd_cases_2000_c, `Receiving MAT in 2000 (count)` = mdd_mat_2000_c,
                                                  `Receiving MAT in 2000 (%)` = mat_prop_2000_c, `MDD cases in 2021` = mdd_cases_2021_c, `Receiving MAT in 2021 (count)` = mdd_mat_2021_c,
                                                  `Receiving MAT in 2021 (%)` = mat_prop_2021_c, `Change in MAT utilisation (%)` = mat_prop_change_c)]

write.xlsx(table_bytime, "/FILEPATH/mat_over_time.xlsx")

# Over time in appendix ---------------------------------------------------

table_bytime_appendix <- rbind(mat_global_all_age, mat_super_regions_all_age, mat_regions_all_age, mat_countries_all_age)

table_bytime_appendix <- merge(table_bytime_appendix[year_id == 2000,.(location_id, draw, mdd_cases_2000 = cases/1000, mdd_mat_2000 = mat_cases/1000, mat_prop_2000 = mat_prop)],
                               table_bytime_appendix[year_id == 2021,.(location_id, draw, mdd_cases_2021 = cases/1000, mdd_mat_2021 = mat_cases/1000, mat_prop_2021 = mat_prop)], by = c('location_id', "draw"))

table_bytime_appendix[, mat_prop_change := (mat_prop_2021 - mat_prop_2000)/mat_prop_2000]

table_bytime_appendix[, `:=` (mdd_cases_2000_c = paste0(round_c(mean(mdd_cases_2000), 1), " (", round_c(quantile(mdd_cases_2000, 0.025), 1), "–", round_c(quantile(mdd_cases_2000, 0.975), 1), ")")), by = 'location_id']
table_bytime_appendix[, `:=` (mdd_mat_2000_c = paste0(round_c(mean(mdd_mat_2000), 1), " (", round_c(quantile(mdd_mat_2000, 0.025), 1), "–", round_c(quantile(mdd_mat_2000, 0.975), 1), ")")), by = 'location_id']
table_bytime_appendix[, `:=` (mat_prop_2000_c = paste0(round_c(mean(mat_prop_2000)*100, 1), " (", round_c(quantile(mat_prop_2000, 0.025)*100, 1), "–", round_c(quantile(mat_prop_2000, 0.975)*100, 1), ")")), by = 'location_id']
table_bytime_appendix[, `:=` (mdd_cases_2021_c = paste0(round_c(mean(mdd_cases_2021), 1), " (", round_c(quantile(mdd_cases_2021, 0.025), 1), "–", round_c(quantile(mdd_cases_2021, 0.975), 1), ")")), by = 'location_id']
table_bytime_appendix[, `:=` (mdd_mat_2021_c = paste0(round_c(mean(mdd_mat_2021), 1), " (", round_c(quantile(mdd_mat_2021, 0.025), 1), "–", round_c(quantile(mdd_mat_2021, 0.975), 1), ")")), by = 'location_id']
table_bytime_appendix[, `:=` (mat_prop_2021_c = paste0(round_c(mean(mat_prop_2021)*100, 1), " (", round_c(quantile(mat_prop_2021, 0.025)*100, 1), "–", round_c(quantile(mat_prop_2021, 0.975)*100, 1), ")")), by = 'location_id']
table_bytime_appendix[, `:=` (mat_prop_change_c = paste0(round_c(mean(mat_prop_change)*100, 1), " (", round_c(quantile(mat_prop_change, 0.025)*100, 1), "–", round_c(quantile(mat_prop_change, 0.975)*100, 1), ")")), by = 'location_id']

table_bytime_appendix <- unique(table_bytime_appendix[,.(location_id, mdd_cases_2000_c, mdd_mat_2000_c, mat_prop_2000_c, mdd_cases_2021_c, mdd_mat_2021_c, mat_prop_2021_c, mat_prop_change_c)])

table_bytime_appendix <- merge(table_bytime_appendix, locations[level < 4,.(location_id, location_name, sort_order)], by = c("location_id"))
table_bytime_appendix <- table_bytime_appendix[order(sort_order), .(location_name, `MDD cases in 2000` = mdd_cases_2000_c, `Receiving MAT in 2000 (count)` = mdd_mat_2000_c,
                                                                    `Receiving MAT in 2000 (%)` = mat_prop_2000_c, `MDD cases in 2021` = mdd_cases_2021_c, `Receiving MAT in 2021 (count)` = mdd_mat_2021_c,
                                                                    `Receiving MAT in 2021 (%)` = mat_prop_2021_c, `Change in MAT utilisation (%)` = mat_prop_change_c)]

write.xlsx(table_bytime_appendix, "/FILEPATH/mat_over_time_appendix.xlsx")


# In-text data by regions  -------------------------------------------
table_data_country <- mat_countries_all_age[year_id == 2021, .(mat_prop = mean(mat_prop), lower = quantile(mat_prop, 0.025), upper = quantile(mat_prop, 0.975)), by = c('location_id')]
table_data_country_males <- mat_countries_all_age[year_id == 2021, .(location_id, mat_prop = mean(mat_prop_m), lower = quantile(mat_prop_m, 0.025), upper = quantile(mat_prop_m, 0.975)), by = c('location_id')]
table_data_country_females <- mat_countries_all_age[year_id == 2021, .(location_id, mat_prop = mean(mat_prop_f), lower = quantile(mat_prop_f, 0.025), upper = quantile(mat_prop_f, 0.975)), by = c('location_id')]


## Highest country
table_data_country[mat_prop == max(mat_prop), .(location_id,mat_prop, lower, upper)]
locations[location_id == 101]

## Lowest country
table_data_country[mat_prop == min(mat_prop), .(location_id,mat_prop, lower, upper)]
locations[location_id == 213]

in_text_data <- mat_countries_all_age[year_id == 2021, .(mat_prop = mean(mat_prop), lower = quantile(mat_prop, 0.025), upper = quantile(mat_prop, 0.975)), by = c('location_name', 'super_region_name', 'region_name')]

in_text_data[super_region_name == "High-income" & mat_prop > 0.3, location_name]
length(in_text_data[super_region_name != "High-income" & mat_prop < 0.05, location_name])
length(in_text_data[super_region_name != "High-income" & mat_prop >= 0.05 & mat_prop <= 0.09, location_name])

map_data <- mat_countries_all_age[year_id == 2021, .(mat_prop = mean(mat_prop), lower = quantile(mat_prop, 0.025), upper = quantile(mat_prop, 0.975)), by = c('location_id')]
map_data[, mapvar := mat_prop*100]
map_data[location_id == 4749, location_id := 95] 
map_data <- map_data[!(location_id %in% c(433, 434, 4636)),]


bins_start <- seq(0, max(map_data$mapvar), max(map_data$mapvar)/10)
bins <- round(bins_start, 1)
bins[1] <- bins[1]-1
bins[11] <- bins[11]+1
binlabels <- c(paste(round(min(map_data$mapvar),1), "to <", bins[2]), paste(bins[2], "to <", bins[3]), paste(bins[3], "to <", bins[4]), paste(bins[4], "to <", bins[5]), paste(bins[5], "to <", bins[6]), paste(bins[6], "to <", bins[7]), paste(bins[7], "to <", bins[8]), paste(bins[8], "to <", bins[9]), paste(bins[9], "to <", bins[10]), paste(bins[10], "to", round(max(map_data$mapvar),1)))

gbd_map(data = map_data, limits = bins, sub_nat="none", legend=TRUE, inset=T,
        labels=binlabels,
        pattern=NULL,
        col="RdYlBu", # Spectral or RdYlBu
        col.reverse=F, na.color = "lightgrey",
        title="\nMDD cases with access to minimally adequate treatment, 2021", fname="/FILEPATH/map.pdf",
        legend.title="Percentage (%)", legend.columns = NULL, legend.cex=2, legend.shift=c(0,-2))

write.csv(map_data, "/FILEPATH/map_data.csv", row.names = F)

# Plot by age -------------------------------------------------------------
mat_global_by_age <- mat_global[year_id == 2021,]

mat_global_by_age <- merge(mat_global_by_age, ages[,.(age_group_id, age_group_name, age_start = age_group_years_start, age_end = age_group_years_end)], by = 'age_group_id', all.x = T)

# make bar graph stackable appropriately by removing 
mat_global_by_age[, `:=` (cases_m = cases_m - mat_cases_m, cases_f = cases_f - mat_cases_f)] 

mat_under5 <- mat_global_by_age[age_end < 6, ]
#mat_under5[, `:=` (cases = prev * population, mat_cases = mat_prop_adj * population)]
mat_under5[, `:=` (cases_m = sum(cases_m), mat_cases_m = sum(mat_cases_m), population_m = sum(population_m),
                   cases_f = sum(cases_f), mat_cases_f = sum(mat_cases_f), population_f = sum(population_f),
                   cases = sum(cases), mat_cases = sum(mat_cases), population = sum(population)), by = c("location_id", "draw")]
#mat_under5[, `:=` (prev = cases / population, mat_prop_adj = mat_cases / population)]

mat_under5[, `:=` (age_group_id = 2, age_group_name = "Under 5", age_start = 0, age_end = 5)]

mat_under5 <- unique(mat_under5)

mat_global_by_age <- rbind(mat_under5, mat_global_by_age[age_end > 5, ])

mat_global_by_age <- rbind(mat_global_by_age[,.(Sex = "Male", age_group_id, Estimate = "All MDD cases", draw, Age = age_group_name, age_start, rate = cases_m / population_m, population = population_m)],
                           mat_global_by_age[,.(Sex = "Female", age_group_id, Estimate = "All MDD cases", draw, Age = age_group_name, age_start, rate = cases_f / population_f, population = population_f)],
                           mat_global_by_age[,.(Sex = "Male", age_group_id, Estimate = "MDD cases receiving MAT", draw, Age = age_group_name, age_start, rate = mat_cases_m / population_m, population = population_m)],
                           mat_global_by_age[,.(Sex = "Female", age_group_id,  Estimate = "MDD cases receiving MAT", draw, Age = age_group_name, age_start, rate = mat_cases_f / population_f, population = population_f)])


mat_global_by_age[, `:=` (rate = mean(rate*100000), lower = quantile(rate*100000, 0.025), upper = quantile(rate*100000, 0.975)), by = c("Sex", "Estimate", "age_group_id")]


mat_global_by_age <- unique(mat_global_by_age[,.(Sex, Estimate, Age, age_start, rate, lower, upper)])

plot <- ggplot(data=mat_global_by_age, aes(x=as.numeric(age_start), y=rate, color = Estimate, fill = Estimate)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(x=as.numeric(age_start), ymin=lower, ymax=upper), colour="black")+
  
  scale_color_manual(values = c("All MDD cases" = "#e62e00", "MDD cases receiving MAT" = "#0035a8"))+
  scale_fill_manual(values = c("All MDD cases" = "#e62e00", "MDD cases receiving MAT" = "#0035a8"))+
  facet_grid(~ Sex ,scales="free_y")+
  ylab("Rate per 100 000 persons") +
  xlab("Age (years)") +
  scale_x_continuous(expand=c(0,0), breaks = unique(mat_global_by_age$age_start), labels = unique(mat_global_by_age$age_start))+
  scale_y_continuous(expand=c(0,0), breaks = seq(0, 6500, 500), limits = c(0, 6500)) +
   theme(axis.line=element_line(colour="black"))
plot


ggsave(plot, filename="/FILEPATH/mat_byage_global_bar.pdf", width = 10, height = 6)

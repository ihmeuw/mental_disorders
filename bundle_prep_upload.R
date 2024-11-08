###########################################################################################
#### Make changes to bundle, save new bundle version, and update bundle meta-data file ####
###########################################################################################

rm(list=ls())

library(data.table)
library(openxlsx)
library(msm)
library(reticulate)
reticulate::use_python("/FILEPATH/python")
mr <- import("mrtool")

acause <- 'mental_unipolar_mdd'

cause_meta_data <- data.table(acause_label = 'mental_unipolar_mdd', bundle_id = 10051)

bundle_id <- cause_meta_data[acause_label == acause, bundle_id]

## Load required functions
source("/FILEPATH/get_bundle_data.R")
source("/FILEPATH/upload_bundle_data.R")
source("/FILEPATH/save_bundle_version.R")
source("/FILEPATH/get_bundle_version.R")
source("/FILEPATH/get_crosswalk_version.R")
source("/FILEPATH/save_crosswalk_version.R")
source("/FILEPATH/get_covariate_estimates.R")
source("/FILEPATH/get_population.R")
rlogit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}
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

# Prep data upload --------------------------------------------------------
data <- data.table(read.xlsx(('/FILEPATH/input_data.xlsx'), sheet = "cleaned"))
data <- data[cv_psychiatrist_only == 0 & cv_psychologist_only == 0 & cv_antidep_psychiatrist == 0, ]

data[, cv_count := cv_antidep + cv_antidep_psychiatrist + cv_psychiatrist_only + cv_psychologist_only + cv_any_mental + cv_mat]
data[cv_count > 1, ]

data[is.na(standard_error) & !is.na(lower), standard_error := (upper-lower)/3.92]

data[, `:=` (location_name = country, location_id = loc_id)]
beijing_data <- data[site %in% c("Beijing and Shanghai", "Beijing, Shanghai"),]
shanghai_data <- data[site %in% c("Beijing and Shanghai", "Beijing, Shanghai"),]
beijing_data[, `:=` (location_name = "Beijing", location_id = 492, standard_error = NA, denominator = denominator * (2633/5201), numerator = NA)]
shanghai_data[, `:=` (location_name = "Shanghai", location_id = 514, standard_error = NA, denominator = denominator * (2568/5201), numerator = NA)]
beijing_data <- rbind(beijing_data, shanghai_data)
beijing_data[, `:=` (numerator = denominator * mean)]
data <- rbind(data[!(site %in% c("Beijing and Shanghai", "Beijing, Shanghai")),], beijing_data)

data[, unique_label := paste0(site, "_", country, "_", sex, "_", age_start, "_", age_end, "_", year_start, "_", year_end)]
data[, unique_study := paste0(site, "_", country, "_", year_start, "_", year_end)]

## For validation
data[study_id > 100, nid := study_id] 

data[urbanicity == "Mixed", urbanicity := "Mixed/both"]

data[recall_service_months == 12, `Recall period` := "12 months"]
data[, `:=` (cv_mat_lenient = 0)]

data <- data[!(cv_mat == 1 & mat_definition_new == "moderate")]

data[cv_mat == 1 & mat_definition_new != "stringent", `:=` (cv_mat_lenient = 1, cv_mat = 0)]

data <- data[cv_mat_lenient == 0, ]

data[, `:=` (recall_type_value = recall_service_months)]
data[recall_type_value == 0, `:=` (recall_type_value = 1)]
data[, cv_recall_1yr := ifelse(recall_type_value < 12, 0, 1)]


data <- data[,.(underlying_nid = NA, nid, field_citation_value = citation, underlying_field_citation_value = NA, file_path = NA, page_num = NA, table_num = NA,	source_type	= "Survey - cross-sectional",
                location_name, location_id, ihme_loc_id = NA, smaller_site_unit = ifelse(location_name == site, 0, 1), site_memo = site, sex, prop_female = pcent_female, sex_issue = 0, year_start, year_end, year_issue = 0,
                age_start, age_end, age_issue = 0, age_demographer = 1, measure = "proportion", mean, lower, upper,
                standard_error, cases = numerator, sample_size = denominator, unit_type = "Person", unit_value_as_published = 1, measure_issue = 0, measure_adjustment	= 0,
                uncertainty_type = NA, uncertainty_type_value = ifelse(is.na(lower), NA, 95), representative_name = ifelse(coverage == "National", "Nationally representative only", "Representative for subnational location only"),
                urbanicity_type = urbanicity, recall_type = "Period: months", recall_type_value, sampling_type = "Multistage", response_rate, case_name	= Unit.of.disaggregation_broad_new,
                case_definition = Diagnostic.clasification,	case_diagnostics = Diagnostic.instrument, group = study_id,	specificity = "not_specified",
                group_review = 1, note_modeler = NA, note_sr = notes, extractor = Extractor, bundle_id = 10051, bundle_name = "Major depressive disorder treatment coverage (any mental health service)",
                seq = NA, effective_sample_size = denominator, design_effect = NA, is_outlier = 0,	cv_mat, cv_mat_lenient, cv_antidep_psychiatrist, cv_antidep, cv_psychiatrist_only, cv_psychologist_only, cv_any_mental, cv_whs = whs_flag,
                cv_wmhs = wmhs_flag, cv_recall_1yr, modelable_entity_name = NA, modelable_entity_id	= NA, seq_parent = NA, input_type = NA)]

data[, note_sr := gsub("\x91", "", note_sr)]
data[, note_sr := gsub("\x93", "", note_sr)]
data[, note_sr := gsub("\x94", "", note_sr)]
data[, note_sr := gsub("\x92", "", note_sr)]

data[sex == "Persons", sex := "Both"]


non_proper_char <- c("Ã", "ƒ", "Æ", "’", "‚", "Â", "¢", "â", "¬", "Å", "¡", "¾", "†", "€", "™", "„", "š", "ž", "¦", "…", "œ")
for(n in non_proper_char){
  data[, `:=` (site_memo = gsub(n, "", site_memo))]
  for(note in names(data)[names(data) %like% "note_"]){
    data[, paste0(note) := gsub(n, "", get(note))]
    data[nchar(get(note)) > 1999, paste0(note) := substring(get(note), 1, 1999)]
  }
}

bundle_save_folder <- paste0("/FILEPATH/")
bundle_save_file <- paste0(bundle_save_folder, "bundle_", gsub(":", "-", gsub(" ", "-", Sys.time())), ".xlsx")
write.xlsx(data, bundle_save_file, sheetName = "extraction")

## upload edits
upload_bundle_data(bundle_id = 10051, filepath = bundle_save_file)

v_id <- save_bundle_version(bundle_id = bundle_id)

v_id <- v_id$bundle_version_id

v_id  
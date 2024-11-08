source("/FILEPATH/get_crosswalk_version.R")
source("/FILEPATH/save_crosswalk_version.R")
source("/FILEPATH/get_location_metadata.R")
source("/FILEPATH/get_bundle_version.R")

v_id <- 42865

raw_data <- get_bundle_version(v_id)
raw_data[, id := paste0(nid, year_start, year_end, age_start, age_end, location_id)]

data <- get_crosswalk_version(45853)
locatins <- get_location_metadata(9, release_id = 9)

set.seed(212015452) ## mashed numpad with fist 3 times 16/9/2024 - Damian

data[, country := location_name]

data[location_name == "Beijing", country := "China"]
data[location_name == "Shanghai", country := "China"]
data[location_name == "New York", country := "United States of America"]
data[location_name == "Maryland", country := "United States of America"]
data[location_name == "São Paulo", country := "Brazil"]
data[location_name == "Minas Gerais", country := "Brazil"]

countries_to_keep <- sample(unique(data$country), 31*0.8, replace=F)

data_retained <- data[country %in% countries_to_keep,]
data_excluded <- data[!(country %in% countries_to_keep),]

data_retained[, note_sr := ""]

data_retained[, id := paste0(nid, year_start, year_end, age_start, age_end, location_id)]
data_retained <- merge(data_retained, raw_dat[, .(seq, id)], by = c('id'), all.x = T)
data_retained[seq.x != seq.y, crosswalk_parent_seq := seq.y]
data_retained[, seq.y := NULL]
setnames(data_retained, "seq.x", "seq")

crosswalk_save_folder <- paste0("/FILEPATH/")
dir.create(file.path(crosswalk_save_folder), showWarnings = FALSE)
crosswalk_save_file <- paste0(crosswalk_save_folder, "crosswalk_mat_out_of_country_validation.xlsx")
write.xlsx(data_retained, crosswalk_save_file, sheetName = "extraction")

## Upload crosswalked dataset to database
save_results <- save_crosswalk_version(v_id, crosswalk_save_file, description = paste0("Out of country validation", Sys.time()))

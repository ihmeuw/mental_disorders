library(data.table)
library(msm)

source("/FILEPATH/mr_brt_functions.R")

dataset <- fread("/FILEPATH/anxiety_mrbrt.csv")
dataset[, se := sqrt(vi)]

length(unique(dataset[d_behavioural_therapy == 1 & d_cognitive_therapy == 1, Study.Name])) ## Number of CBT studies
length(unique(dataset[d_behavioural_therapy == 0 & d_cognitive_therapy == 1, Study.Name])) ## Number of cognitive studies
length(unique(dataset[d_behavioural_therapy == 1 & d_cognitive_therapy == 0, Study.Name])) ## Number of behavioural studies

## Run network meta-analysis in MR-BRT ##

fit1 <- run_mr_brt(
  output_dir = "/FILEPATH/treatment_effect_size/",
  model_label = "treatment",
  data = dataset,
  mean_var = "yi",
  se_var = "se",
  covs = list(
    cov_info("d_cognitive_therapy", "X"),
    cov_info("d_behavioural_therapy", "X"),
    cov_info("d_antidepressants", "X"),
    cov_info("d_antipsychotics", "X"),
    cov_info("d_azapirones", "X"),
    cov_info("d_anticonvulsants", "X"),
    cov_info("d_benzodiazepines", "X"),
    cov_info("d_neurokinin", "X"),
    cov_info("d_betablocker", "X"),
    cov_info("d_lithium", "X"),
    cov_info("d_antihistamine", "X"),
    cov_info("d_cycloserine", "X"),
    cov_info("d_psychodynamic", "X"),
    cov_info("d_supportive", "X")),
  remove_x_intercept = TRUE,
  method = "trim_maxL",
  trim_pct = 0.1,
  study_id = "Study.Name",
  overwrite_previous = TRUE,
  lasso = FALSE)

check_for_outputs(fit1)
results1 <- load_mr_brt_outputs(fit1)
df_pred1 <- expand.grid(d_cognitive_therapy = c(0,1), d_behavioural_therapy = c(0, 1), d_antidepressants = c(0,1),
                        d_antipsychotics = c(0, 1), d_azapirones = c(0, 1), d_anticonvulsants = c(0, 1),
                        d_benzodiazepines = c(0, 1), d_neurokinin = c(0, 1), d_betablocker = c(0, 1), d_lithium = c(0, 1),
                        d_antihistamine = c(0, 1), d_cycloserine = c(0, 1), d_psychodynamic = c(0, 1), d_supportive = c(0, 1)
)
pred1 <- as.data.table(predict_mr_brt(fit1, newdata = df_pred1)["model_summaries"])

df_optimal <- expand.grid(d_cognitive_therapy = 1, d_behavioural_therapy = 1, d_antidepressants = 1,
                d_antipsychotics = 0, d_azapirones = 0, d_anticonvulsants = 0,
                d_benzodiazepines = 0, d_neurokinin = 0, d_betablocker = 0, d_lithium = 0,
                d_antihistamine = 0, d_cycloserine = 0, d_psychodynamic = 0, d_supportive = 0)

pred_optimal <- as.data.table(predict_mr_brt(fit1, newdata = df_optimal)["model_summaries"])
pred_optimal[,.(mean = model_summaries.Y_mean, lower = model_summaries.Y_mean_lo, upper = model_summaries.Y_mean_hi)]

names(pred1) <- gsub("model_summaries.", "", names(pred1))
names(pred1) <- gsub("X_", "", names(pred1))

write.csv(pred1, "/FILEPATH/treatment_effects_predicted.csv", row.names=F)

pred1[, count := apply(.SD, 1, sum), .SDcols = names(pred1)[grepl("d_", names(pred1))]]
pred1 <- pred1[(count == 1 & d_cognitive_therapy == 0 & d_behavioural_therapy == 0 )| (count == 2 & d_cognitive_therapy == 1 & d_behavioural_therapy == 1), ]
for(n in names(pred1)[grepl("d_", names(pred1))]){
  pred1[get(n) == 1, treatment := gsub("d_", "", n)]
  pred1[, paste0(n) := NULL]
}

write.csv(pred1, "/FILEPATH/treatment_effect_summaries.csv", row.names=F)



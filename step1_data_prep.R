
library(data.table)
library(openxlsx)
library(metafor)

dataset <- data.table(read.xlsx("/FILEPATH/anxiety_extraction_sheet.xlsx"))
dataset[Improvement.Direction == 1, `:=` (Intervention.M = -Intervention.M, Comparison.M = -Comparison.M)]

dataset <- data.table(escalc(measure = "SMD", m1i = Intervention.M, sd1i = Intervention.SD, n1i = Intervention.N,
                                m2i = Comparison.M, sd2i = Comparison.SD, n2i = Comparison.N, data = dataset))

## If only SMD and SE was reported
dataset[is.na(yi), `:=` (yi =Intervention.M, vi = Intervention.SD^2)]

## Exclude duplicate effect sizes
dataset[, `:=` (yi_round = round(yi, 1), vi_round = round(vi, 1))]
dataset <- unique(dataset, by = c("Study.Name", "yi_round", "vi_round"))

## Create intervention dummy variables ##
dataset[, `:=` (d_cognitive_therapy = 0, d_behavioural_therapy = 0, d_antidepressants = 0, 
                d_antipsychotics = 0, d_azapirones = 0, d_anticonvulsants = 0, 
                d_benzodiazepines = 0, d_neurokinin = 0, d_betablocker = 0, d_lithium = 0, 
                d_antihistamine = 0, d_cycloserine = 0, d_psychodynamic = 0, d_supportive = 0)]

# Psychotherapy
dataset[grepl("CBT", Intervention), `:=` (d_cognitive_therapy = 1, d_behavioural_therapy = 1)]
dataset[grepl("Cognitive therapy", Intervention), `:=` (d_cognitive_therapy = 1)]
dataset[grepl("Behavioural therapy", Intervention), `:=` (d_behavioural_therapy = 1)]
dataset[grepl("Psychodynamic therapy", Intervention), `:=` (d_psychodynamic = 1)]
dataset[grepl("Exposure therapy", Intervention), `:=` (d_behavioural_therapy = 1)]

# Medications
antidepressants <- c("Noradrenalin Reuptake Inhibitors", "Tetracyclic antidepressants", "SARI", "SNRI", "SSRI", "Monoamine oxidase inhibitors", "Tricyclic antidepressants")
for(i in antidepressants){
  dataset[grepl(i, Intervention), `:=` (d_antidepressants = 1)]  
}
dataset[grepl("Antipsychotics", Intervention), `:=` (d_antipsychotics = 1)]
dataset[grepl("Azapirones", Intervention), `:=` (d_azapirones = 1)]
dataset[grepl("Anticonvulsants", Intervention), `:=` (d_anticonvulsants = 1)]
dataset[grepl("Benzodiazepines", Intervention), `:=` (d_benzodiazepines = 1)]
dataset[grepl("Neurokinin", Intervention), `:=` (d_neurokinin  = 1)]
dataset[grepl("Beta", Intervention), `:=` (d_betablocker  = 1)]
dataset[grepl("Lithium", Intervention), `:=` (d_lithium  = 1)]
dataset[grepl("Antihistamine", Intervention), `:=` (d_antihistamine  = 1)]
dataset[grepl("cycloserine", Intervention), `:=` (d_cycloserine  = 1)]

# Code comparison dummies #
dataset[grepl("Benzodiazepines", Comparison), `:=` (d_benzodiazepines = d_benzodiazepines-1)]
dataset[grepl("SSRI", Comparison), `:=` (d_antidepressants = d_antidepressants-1)]
dataset[grepl("Psychodynamic therapy", Comparison), `:=` (d_psychodynamic = d_psychodynamic-1)]
dataset[grepl("Behavioural therapy", Comparison), `:=` (d_behavioural_therapy = d_behavioural_therapy-1)]
dataset[grepl("Prolonged exposure", Comparison), `:=` (d_behavioural_therapy = d_behavioural_therapy-1)]
dataset[grepl("CBT", Comparison), `:=` (d_cognitive_therapy = d_cognitive_therapy - 1, d_behavioural_therapy = d_behavioural_therapy-1)]
dataset[grepl("Azapirones", Comparison), `:=` (d_azapirones = d_azapirones - 1)]
dataset[grepl("Supportive therapy", Comparison), `:=` (d_supportive = d_supportive - 1)]

# Write out datafile #
write.csv(dataset, "/FILEPATH/anxiety_mrbrt.csv", row.names=F)

length(dataset[, Study.Name]) # Number of estimates
length(dataset[, unique(Study.Name)]) # Number of studies



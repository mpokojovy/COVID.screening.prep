####################################################################################################
# Copyright (C) 2022 Michael Pokojovy                                                              #
#                                                                                                  #
# Data preparation for COVID-19 patient screening using machine learning                           #                                                                                              #
####################################################################################################

raw_tables_dir = "raw.tables"

all.vitals.df       = read.csv(file = paste(raw_tables_dir, "/", "all_vitals.tsv", sep = ""),            header = TRUE, stringsAsFactors = FALSE, sep = "\t")
cxr_covid_films.df  = read.csv(file = paste(raw_tables_dir, "/", "cxr_covid_films.tsv", sep = ""),       header = TRUE, stringsAsFactors = FALSE, sep = "\t")
all.diagnoses.df    = read.csv(file = paste(raw_tables_dir, "/", "all_diagnoses.tsv", sep = ""),         header = TRUE, stringsAsFactors = FALSE, sep = "\t")
patient.outcomes.df = read.csv(file = paste(raw_tables_dir, "/", "study_patient_outcome.tsv", sep = ""), header = TRUE, stringsAsFactors = FALSE, sep = "\t")
demographics.df     = read.csv(file = paste(raw_tables_dir, "/", "all_demographics.tsv", sep = ""),      header = TRUE, stringsAsFactors = FALSE, sep = "\t")

## Coerce data types

all.diagnoses.df$DATE         = as.Date(all.diagnoses.df$DATE)
cxr_covid_films.df$study_date = as.Date(cxr_covid_films.df$study_date)
all.vitals.df$RECORDED_TIME   = as.Date(all.vitals.df$RECORDED_TIME)

cxr_covid_films.df$tdelta_days_testing_to_cxr[is.na(cxr_covid_films.df$tdelta_days_testing_to_cxr)] = 0.0

## Prepare vitals and outcomes

all.diagnoses.df$CURRENT_ICD10_LIST = trimws(all.diagnoses.df$CURRENT_ICD10_LIST, which = "both")
ICD10.unique = sort(unique(all.diagnoses.df$CURRENT_ICD10_LIST))
ICD10.unique.names = rep("", length(ICD10.unique))

# ## Remove all COVID-related ICD codes
# COVID.ICD10 = c()
# for (i in 1:nrow(all.diagnoses.df)) {
#   disease.name = toupper(all.diagnoses.df[i, "DIAGNOSIS.NAME"])
#   if (length(grep("COVID", disease.name)) || length(grep("CORONAVIRUS", disease.name))) {
#     COVID.ICD10 = union(COVID.ICD10, all.diagnoses.df[i, "CURRENT_ICD10_LIST"])
#   }
# }
# [1] "U07.1"    "B34.2"    "Z20.828"  " J12.89"  " J40"     " J22"     " J06.9"   "Z86.19"   "Z01.812"  " Z20.828" " J20.8"   " J80"     " J98.8"   " J98.4"  
# [15] "Z71.89"   "B97.29"   " A08.39"  "Z91.89"   "J12.81"  
# 
# U07.1 COVID-19
# B34.2 Coronavirus infection, unspecified
# Z20.828 Contact with and (suspected) exposure to other viral communicable diseases
# J12.89 Other viral pneumonia
# J40 Bronchitis, not specified as acute or chronic
# J22 Unspecified acute lower respiratory infection
# J06.9 Acute upper respiratory infection, unspecified
# Z86.19 Personal history of other infectious and parasitic diseases
# Z01.812 Encounter for preprocedural laboratory examination
# Z20.828 Contact with and (suspected) exposure to other viral communicable diseases
# J20.8 Acute bronchitis due to other specified organisms
# J80 Acute respiratory distress syndrome
# J98.8 Other specified respiratory disorders
# J98.4 Other disorders of lung
# Z71.89 Other specified counseling
# B97.29 Other coronavirus as the cause of diseases classified elsewhere
# A08.39 Other viral enteritis
# Z91.89 Other specified personal risk factors, not elsewhere classified
# J12.81 Pneumonia due to SARS-associated coronavirus
#
# ## Exclude:
# U07.1 COVID-19
# B34.2 Coronavirus infection, unspecified
# J12.89 Other viral pneumonia
# B97.29 Other coronavirus as the cause of diseases classified elsewhere
# J12.81 Pneumonia due to SARS-associated coronavirus

COVID.ICD10 = c("U07.1", "B34.2", "J12.89", "B97.29", "J12.81")
ICD10.unique = setdiff(ICD10.unique, COVID.ICD10)

for (i in 1:length(ICD10.unique)) {
  ICD10.unique.names[i] = all.diagnoses.df$DIAGNOSIS.NAME[which(all.diagnoses.df$CURRENT_ICD10_LIST == ICD10.unique[i])[1]]
}

cxr_covid_films.df$Value = factor(cxr_covid_films.df$Value, levels = c("Not Detected", "Positive", "Not Tested"))

I = !is.na(cxr_covid_films.df$tdelta_days_testing_to_cxr)
cxr_covid_films.df$pcr_date = as.Date(cxr_covid_films.df$study_date)
cxr_covid_films.df$pcr_date[I] = cxr_covid_films.df$pcr_date[I] + cxr_covid_films.df$tdelta_days_testing_to_cxr[I]

history.ICD    = 365.0 # 356 days = 1 year
history.vitals = 14.0  # 14 days = 2 weeks
future  = 1.0          # 1 day

prepare.observations = function(IDs) {
  VITALS_DATE = all.vitals.df$RECORDED_TIME
  ICD_DATE    = all.diagnoses.df$DATE
  
  varnames    = c("BMI", "BP_DIASTOLIC", "BP_SYSTOLIC", "PULSE", "PULSE.OXIMETRY", "RESPIRATIONS", "TEMPERATURE")
  varnames.x5 = as.vector(sapply(varnames, function(varname) sapply(1:5, function(ind) paste(varname, ind, sep = ""))))
  
  res = data.frame(matrix(0.0, nrow = length(IDs), ncol = 1L + 4L + 5L*length(varnames) + length(ICD10.unique)))
  colnames(res) = c("hashed_mrn", "SEX", "ETHNICITY", "AGE", "PCR", varnames.x5, ICD10.unique)

  res$hashed_mrn = IDs
  
  for (i in 1:length(IDs)) {
    ID = IDs[i]
    
    I = which(cxr_covid_films.df$hashed_mrn == ID)
    stopifnot(length(I) > 0L)
    
    PCR.Result = cxr_covid_films.df$Value[I[which.max(cxr_covid_films.df$study_date[I] + cxr_covid_films.df$tdelta_days_testing_to_cxr[I])]]
    PCR.Date   = max(cxr_covid_films.df$study_date[I] + cxr_covid_films.df$tdelta_days_testing_to_cxr[I])
    
    # Response: PCR
    res[i, "PCR"] = PCR.Result
    
    # Demographics
    I = which(demographics.df$hashed_mrn == ID)
    
    stopifnot(length(I) <= 1)
    
    if (length(I) == 0) {
      res[i, "SEX"]       = NA
      res[i, "ETHNICITY"] = "UNKNOWN"
      res[i, "AGE"]       = NA
    } else {
      res[i, "SEX"]       = demographics.df[I, "SEX"]
      res[i, "ETHNICITY"] = demographics.df[I, "ETHNICITY"]
      res[i, "AGE"]       = demographics.df[I, "AGE"]
      
      if (res[i, "ETHNICITY"] == "UNKNOWN") {
        res[i, "ETHNICITY"] = NA
      }
    }
    
    # Vitals
    I = which(all.vitals.df$hashed_mrn == ID)
    
    I.Dates = I[which((VITALS_DATE[I] <= PCR.Date + future) & (VITALS_DATE[I] >= PCR.Date - history.vitals))]
    
    if (length(I.Dates) > 0L) {
      res[i, varnames.x5] = as.vector(sapply(1:length(varnames), function(varind) fivenum(all.vitals.df[I.Dates, varnames[varind]], na.rm = TRUE)))
    } else {
      res[i, varnames.x5] = NA
    }
    
    # ICD
    I = which(all.diagnoses.df$hashed_mrn == ID)
    
    I.Dates = I[which((ICD_DATE[I] <= PCR.Date + future) & (ICD_DATE[I] >= PCR.Date - history.ICD))]
    res[i, intersect(ICD10.unique, unique(all.diagnoses.df$CURRENT_ICD10_LIST[I.Dates]))] = TRUE
  }
  
  for (ICD in ICD10.unique) {
    res[, ICD] = factor(res[, ICD], levels = c(0, 1), labels = c("FALSE", "TRUE"))
  }
  
  res$PCR = factor(res$PCR, levels = c(1L, 2L, 3L), labels = c("Not Detected", "Positive", "Not Tested"))
  
  return(res)
}

pat.IDs = unique(cxr_covid_films.df[, "hashed_mrn"])

COVID.diagnosis.df = prepare.observations(pat.IDs)

write.csv(COVID.diagnosis.df, file = "preproc.tables/COVID.diagnosis.df.csv", row.names = FALSE)
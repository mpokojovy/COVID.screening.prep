####################################################################################################
# Copyright (C) 2022 Michael Pokojovy                                                              #
#                                                                                                  #
# Data preparation for COVID-19 patient screening using machine learning                           #                                                                                              #
####################################################################################################


## !!!
# Set working directory
path = NULL # Replace with actual directory path!
setwd(path)
## !!!

githib.path = "https://raw.githubusercontent.com/pmccaffrey6/COVID-LOS/main/"
github.raw.files = c("April_2021_encounters_film.tsv",
                     "all_demographics.tsv",
                     "all_diagnoses.tsv",
                     "all_vitals.tsv",
                     "cxr_covid_films.tsv",
                     "study_patient_outcome.tsv",
                     "study_vitals_at_exam.tsv")

if (!dir.exists("raw.tables")) {
  dir.create("raw.tables")
}

if (!dir.exists("preproc.tables")) {
  dir.create("preproc.tables")
}

if (!dir.exists("imputed.tables")) {
  dir.create("imputed.tables")
}

for (file in github.raw.files) {
  url = paste(githib.path, file, sep = "")
  destfile = paste("raw.tables/", file, sep = "")
  download.file(url, destfile)
}

source("auxil.R")

cat("Preprocessing...\n")
source("preprocess.R")

cat("\nImputing...\n")
source("impute.R")
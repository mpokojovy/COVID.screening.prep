####################################################################################################
# Copyright (C) 2022 Michael Pokojovy                                                              #
#                                                                                                  #
# Data preparation for COVID-19 patient screening using machine learning                           #                                                                                              #
####################################################################################################

COVID.diagnosis.df = read.csv(file = "preproc.tables/COVID.diagnosis.df.csv", stringsAsFactors = FALSE)

## Remove "Not Tested"
I = which(COVID.diagnosis.df$PCR == "Not Tested")
COVID.diagnosis.df = COVID.diagnosis.df[-I, ]

## Coerce factors

COVID.diagnosis.df$PCR       = factor(COVID.diagnosis.df$PCR,       levels = c("Not Detected", "Positive"), labels = c("Not Detected", "Positive"))
COVID.diagnosis.df$SEX       = factor(COVID.diagnosis.df$SEX,       levels = c("MALE", "FEMALE"), labels = c("MALE", "FEMALE"))
COVID.diagnosis.df$ETHNICITY = factor(COVID.diagnosis.df$ETHNICITY, levels = c("NOT HISPANIC OR LATINO", "HISPANIC OR LATINO"), labels = c("NOT HISPANIC OR LATINO", "HISPANIC OR LATINO"))

varnames    = c("BMI", "BP_DIASTOLIC", "BP_SYSTOLIC", "PULSE", "PULSE.OXIMETRY", "RESPIRATIONS", "TEMPERATURE")
varnames.x5 = as.vector(sapply(varnames, function(varname) sapply(1:5, function(ind) paste(varname, ind, sep = ""))))

basic.var.names = c("SEX", "ETHNICITY", "AGE", "PCR", varnames.x5)
ICD10.var.names = setdiff(colnames(COVID.diagnosis.df), c("hashed_mrn", basic.var.names))

for (ICD in ICD10.var.names) {
  COVID.diagnosis.df[, ICD] = factor(COVID.diagnosis.df[, ICD], levels = c(FALSE, TRUE), labels = c("FALSE", "TRUE"))
}

## Split into training & test datasets
train2test.ratio = 0.8

set.seed(1)

n = nrow(COVID.diagnosis.df)

n.train = ceiling(n*train2test.ratio)
n.test  = n - n.train

I.train = sample(1:n, n.train, replace = FALSE)
I.test  = setdiff(1:n, I.train)

COVID.diagnosis.train.df = COVID.diagnosis.df[I.train, ]
COVID.diagnosis.test.df  = COVID.diagnosis.df[I.test, ]

write.csv(COVID.diagnosis.train.df, file = "preproc.tables/COVID.diagnosis.train.df.csv", row.names = FALSE)
write.csv(COVID.diagnosis.test.df,  file = "preproc.tables/COVID.diagnosis.test.df.csv",  row.names = FALSE)

## Detect top ICD10 codes

prev = colSums(ifelse(COVID.diagnosis.train.df[, ICD10.var.names] == "TRUE", 1.0, 0.0))/nrow(COVID.diagnosis.train.df)
top  = 20L
top.ICD10.var.names = ICD10.var.names[order(abs(prev - 0.5))[1:top]]

cat("Top ICD 10 codes: \n")
print(top.ICD10.var.names)

## Vars to impute

imp.var.names = c(basic.var.names, top.ICD10.var.names)

## Impute train

set.seed(2)

for (method in c("cart", "rf")) {
  imp.train = mice::mice(COVID.diagnosis.train.df[, imp.var.names], m = 5L, method = method, maxit = 5L,
                         remove.collinear = FALSE, remove.constant = FALSE)

  COVID.diagnosis.train.imp.df = COVID.diagnosis.train.df
  COVID.diagnosis.train.imp.df[, imp.var.names] = mice::complete(imp.train)

  write.csv(COVID.diagnosis.train.imp.df, file = paste("imputed.tables/COVID.diagnosis.train.imp.", method, ".df.csv", sep = ""),
            row.names = FALSE)
}

## Impute test

set.seed(3)

for (method in c("cart", "rf")) {
  COVID.diagnosis.train.imp.df = read.csv(file = paste("imputed.tables/COVID.diagnosis.train.imp.", method, ".df.csv", sep = ""), 
                                          header = TRUE, stringsAsFactors = FALSE, sep = ",")
  
  COVID.diagnosis.test.imp.df = rbind(COVID.diagnosis.train.imp.df, COVID.diagnosis.test.df)
  
  COVID.diagnosis.test.imp.df$PCR       = factor(COVID.diagnosis.test.imp.df$PCR,       levels = c("Not Detected", "Positive"), labels = c("Not Detected", "Positive"))
  COVID.diagnosis.test.imp.df$SEX       = factor(COVID.diagnosis.test.imp.df$SEX,       levels = c("MALE", "FEMALE"), labels = c("MALE", "FEMALE"))
  COVID.diagnosis.test.imp.df$ETHNICITY = factor(COVID.diagnosis.test.imp.df$ETHNICITY, levels = c("NOT HISPANIC OR LATINO", "HISPANIC OR LATINO"), labels = c("NOT HISPANIC OR LATINO", "HISPANIC OR LATINO"))
  
  for (ICD in ICD10.var.names) {
    COVID.diagnosis.test.imp.df[, ICD] = factor(COVID.diagnosis.test.imp.df[, ICD], levels = c(FALSE, TRUE), labels = c("FALSE", "TRUE"))
  }
  
  imp.test = mice::mice(COVID.diagnosis.test.imp.df[, imp.var.names], m = 5L, method = method, maxit = 5L,
                        remove.collinear = FALSE, remove.constant = FALSE)
    
  COVID.diagnosis.test.imp.df[, imp.var.names] = mice::complete(imp.test)[, imp.var.names]
  COVID.diagnosis.test.imp.df = COVID.diagnosis.test.imp.df[(n.train + 1):n, ]
  
  write.csv(COVID.diagnosis.test.imp.df, file = paste("imputed.tables/COVID.diagnosis.test.imp.", method, ".df.csv", sep = ""),
            row.names = FALSE)
}
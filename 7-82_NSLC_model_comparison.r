
library(tidyverse)
library(data.table)
library(glue)
library(xlsx)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################################################################
######################## Data import and formatting ############################
################################################################################

# Clinicopathological data - keeping only adenocarcinomas
clin_data <- as.data.table(read.xlsx("Data/NCI_NeverSmoker_n92_20210812_TMB_drivers.xlsx", sheetIndex = 1))
clin_data <- clin_data[histology == 3, ]
Patients_list <- gsub("NSLC-", "",sort(clin_data[, Patient_ID]))

# Fetching and merging all patients variants
VEP_data_all_patients <- rbindlist(lapply(Patients_list, function(x) fread(glue("Data/VEP_82_NSLC_TMB/VEP_NSLC-{x}.csv"))), use.names = TRUE)

# TMB data for whole genome and specific regions
TMB <- as.data.table(read.xlsx("Data/VEP_82_NSLC_TMB.xlsx", sheetIndex = 1))

# Merged dataset for survival models
surv.dt <- merge.data.table(clin_data, TMB, by = "Patient_ID")
surv.dt <- surv.dt %>% mutate(TMB_high_low_old = case_when(complete_WGS_TMB >= 1.70 ~ "High",
                                                           complete_WGS_TMB < 1.70 ~ "Low"),
                              pathological_stage_refactor = case_when(pathological_stage %in% c("1A1","1A2","1A3","1B") ~ "I",
                                                                      pathological_stage %in% c("2A", "2B") ~ "II",
                                                                      pathological_stage %in% c("3A", "3B", "4A") ~ "III & IV"),
                              recurrence.two = case_when(time_RpFS >= 730 ~ ">=2",
                                                         time_RpFS < 730 ~ "<2"),
                              recurrence.five = case_when(time_RpFS >= 1825 ~ ">=5",
                                                          time_RpFS < 1825 ~ "<5"),
                              os.two = case_when(time_os >= 730 ~ ">=2",
                                                 time_os < 730 ~ "<2"),
                              os.five = case_when(time_os >= 1825 ~ ">=5",
                                                  time_os < 1825 ~ "<5"))
# Changing comorbidities NA to None
surv.dt[, comorbidities:=replace_na(comorbidities, "None")]


################################################################################
############################## Model building ##################################
################################################################################






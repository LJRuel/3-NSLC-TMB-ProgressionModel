
library(xlsx)
library(data.table)
library(tidyverse)
library(survival)
library(survminer)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Clinicopathological data - keeping only adenocarcinomas
clin_data <- as.data.table(read.xlsx("Data/NCI_NeverSmoker_n92_20210812_TMB_drivers.xlsx", sheetIndex = 1))
clin_data <- clin_data[histology == 3, ]
Patients_list <- gsub("NSLC-", "",sort(clin_data[, Patient_ID]))

# TMB data for whole genome and specific regions
TMB <- as.data.table(read.xlsx("Data/VEP_82_NSLC_TMB.xlsx", sheetIndex = 1))

# Merged dataset for survival models
surv.dt <- merge.data.table(clin_data, TMB, by = "Patient_ID")
surv.dt <- surv.dt %>% mutate(TMB_high_low_old = case_when(complete_WGS_TMB >= 1.70 ~ "High",
                                                           complete_WGS_TMB < 1.70 ~ "Low"),
                              pathological_stage_refactor = case_when(pathological_stage %in% c("1A1","1A2","1A3","1B") ~ "I",
                                                                      pathological_stage %in% c("2A", "2B") ~ "II",
                                                                      pathological_stage %in% c("3A", "3B", "4A") ~ "III"))

# Exporting surv.dt for use with SAS "newsurv" macro.
#colnames(surv.dt) <- gsub("[.]", "_", colnames(surv.dt))
#write.xlsx(surv.dt, "n82_NSLC_surv_data.xlsx", row.names = FALSE)

# Survival models with Overall survival
surv.OS.wg <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.genome_TMB)
surv.OS.exons <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.exons_TMB)
surv.OS.introns <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.introns_TMB)
surv.OS.reg <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.regulatory_TMB)
surv.OS.intergenic <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.intergenic_TMB)

surv.di.simple <-  coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_high_low_old)
surv.di.pstage <-  coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_high_low_old + pathological_stage_refactor)
surv.di.pstage.strata <-  coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_high_low_old + strata(pathological_stage_refactor))


summary(surv.OS.wg)
summary(surv.di.simple)
summary(surv.di.pstage)

plot(predict(surv.di.pstage), 
     residuals(surv.di.pstage, type="martingale"))
abline(h=0)
lines(smooth.spline(predict(surv.di.pstage), 
                    residuals(surv.di.pstage, type="martingale")), col = "red")

plot(cox.zph(surv.di.pstage))
abline(h=0, col="red")

summary(surv.OS.exons)
summary(surv.OS.introns)
summary(surv.OS.reg)
summary(surv.OS.intergenic)


# Survival models with Recurrence free 


################################################################################
############################# Distribution plot ################################
################################################################################

ggplot(surv.dt, aes(x=TMB_per_region.genome_TMB)) +
  geom_density()


################################################################################
############################### Model testing ##################################
################################################################################

model.cox.simple <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.genome_TMB)

model.cox.simple.spline <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ pspline(TMB_per_region.genome_TMB, df=4))
ptemp <- termplot(model.cox.simple.spline, se=TRUE, plot=FALSE)

model.cox.pstage <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.genome_TMB + pathological_stage_refactor)

model.cox.pstage.strata <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.genome_TMB + strata(pathological_stage_refactor))

model.cox.age <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.genome_TMB + age)

model.cox.pstage_age <- coxph(data = surv.dt, Surv(time=time_os, VitalStatus) ~ TMB_per_region.genome_TMB + age + pathological_stage_refactor)

anova(model.cox.simple, model.cox.pstage, test="LRT") # The two models do not seem different, therefore we can exclude pathological stage

## Check for linearity for the simple model
# Using Martingale residuals
plot(predict(model.cox.simple), 
     residuals(model.cox.simple, type="martingale"))
abline(h=0)
lines(smooth.spline(predict(model.cox.simple), 
                    residuals(model.cox.simple, type="martingale")), col = "red")

# Using deviance residuals
plot(predict(model.cox.simple), 
     residuals(model.cox.simple, type="deviance"))
abline(h=0)
lines(smooth.spline(predict(model.cox.simple), 
                    residuals(model.cox.simple, type="deviance")), col = "red")



## Check for proportional hazards assumption for the simple model
plot(cox.zph(model.cox.simple))
abline(h=0, col="red")

## Check for proportional hazards assumption for the model with pathological stage
plot(cox.zph(model.cox.pstage)[1])
abline(h=0, col="red")
plot(cox.zph(model.cox.pstage)[2]) # The hazards are not proportional for the pathological stage
abline(h=0, col="red")

## Check for proportional hazards assumption for the model with stratified pathological stage
plot(cox.zph(model.cox.pstage.strata))
abline(h=0, col="red")



## TODO
## - Correct the assumptions of the cox model: the proportional hazards assumption is met for the simple model, but not the linearity. ** Correct the linearity
## - Make Cox models for drivers
## - Look at survival and reccurence (<2 years, >= 2 years)
## - Check the count of unique mutations throughout the cohort. Then, separate the cohort into high and low OS to check the ratio of mutations (Sebastien's idea).



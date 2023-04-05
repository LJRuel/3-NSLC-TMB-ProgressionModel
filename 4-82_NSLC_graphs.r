
library(GenomicFeatures)
library(tidyverse)
library(data.table)
library(glue)
library(xlsx)
library(karyoploteR)
library(regioneR)
library(corrplot)
library(Hmisc)
library(ggsci)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
############################# Distribution plot ################################
################################################################################

ggplot(surv.dt, aes(x=TMB_per_region.genome_TMB)) +
  geom_density() +
  geom_vline(xintercept = 1.7)

################################################################################
############################# Correlation matrix ###############################
################################################################################

# Creating new data table where non-numerical variables are made binary
surv.dt.split <- surv.dt %>% mutate(comorbidities.Hypertension = case_when(comorbidities=="Hypertension" ~ 1,
                                                                           comorbidities!="Hypertension" ~ 0),
                                    comorbidities.None = case_when(comorbidities=="None" ~ 1,
                                                                   comorbidities!="None" ~ 0),
                                    comorbidities.Diabetes = case_when(comorbidities=="Diabetes" ~ 1,
                                                                       comorbidities!="Diabetes" ~ 0),
                                    comorbidities.Asthma = case_when(comorbidities=="Asthma" ~ 1,
                                                                     comorbidities!="Asthma" ~ 0),
                                    comorbidities.Emphysema = case_when(comorbidities=="Emphysema" ~ 1,
                                                                        comorbidities!="Emphysema" ~ 0),
                                    pathological_stage_refactor.I = case_when(pathological_stage_refactor=="I" ~ 1,
                                                                              pathological_stage_refactor!="I" ~ 0),
                                    pathological_stage_refactor.II = case_when(pathological_stage_refactor=="II" ~ 1,
                                                                               pathological_stage_refactor!="II" ~ 0),
                                    pathological_stage_refactor.III.IV = case_when(pathological_stage_refactor=="III & IV" ~ 1,
                                                                                   pathological_stage_refactor!="III & IV" ~ 0),
                                    Surgery_type.Lobectomy = case_when(Surgery_type=="Lobectomy" ~ 1,
                                                                       Surgery_type!="Lobectomy" ~ 0),
                                    Surgery_type.Pneumonectomy = case_when(Surgery_type=="Pneumonectomy" ~ 1,
                                                                           Surgery_type!="Pneumonectomy" ~ 0),
                                    Surgery_type.Bilobectomy = case_when(Surgery_type=="Bilobectomy" ~ 1,
                                                                         Surgery_type!="Bilobectomy" ~ 0),
                                    Surgery_type.Wedge_resection = case_when(Surgery_type=="Wedge resection" ~ 1,
                                                                             Surgery_type!="Wedge resection" ~ 0),
                                    Surgery_type.Segmentectomy = case_when(Surgery_type=="Segmentectomy" ~ 1,
                                                                           Surgery_type!="Segmentectomy" ~ 0),
                                    TMB.High = case_when(TMB_high_low_old=="High" ~ 1,
                                                         TMB_high_low_old!="High" ~ 0),
                                    TMB.Low = case_when(TMB_high_low_old=="Low" ~ 1,
                                                        TMB_high_low_old!="Low" ~ 0),
                                    recurrence.two.above = case_when(recurrence.two == ">=2" ~ 1,
                                                                     recurrence.two != ">=2" ~ 0),
                                    recurrence.two.below = case_when(recurrence.two == "<2" ~ 1,
                                                                     recurrence.two != "<2" ~ 0),
                                    recurrence.five.above = case_when(recurrence.five == ">=5" ~ 1,
                                                                      recurrence.five != ">=5" ~ 0),
                                    recurrence.five.below = case_when(recurrence.five == "<5" ~ 1,
                                                                      recurrence.five != "<5" ~ 0))

surv.dt.split <- surv.dt.split[, .(time_os, VitalStatus, TMB_per_region.genome_TMB,
                                   TMB.High, TMB.Low, sex, age, BMI,
                                   recurrence.two.above, recurrence.two.below,
                                   recurrence.five.above, recurrence.five.below,
                                   passive.smoking, pathological_stage_refactor.I,
                                   pathological_stage_refactor.II, pathological_stage_refactor.III.IV,
                                   Surgery_type.Lobectomy, Surgery_type.Pneumonectomy, 
                                   Surgery_type.Bilobectomy, Surgery_type.Wedge_resection,
                                   Surgery_type.Segmentectomy, tumor_size, comorbidities.None, 
                                   comorbidities.Hypertension, comorbidities.Diabetes, 
                                   comorbidities.Asthma, comorbidities.Emphysema, 
                                   EGFR, ERBB2, KRAS, BRAF, TP53, PIK3CA, MET, driver)]

# Correlation excluding RpFS data, because 2 patients don't have RpFS follow-up
cor.m <- rcorr(as.matrix(surv.dt.split),
               type = "pearson")
cor.m$P <- cor.m$P %>% replace(is.na(.), 0)

png("Results/Misc plots/correlation_matrix.png", width = 10, height = 10, units = "in", res = 1200)
corrplot(cor.m$r, tl.col = "black", p.mat = cor.m$P, insig = "blank")
dev.off()

################################################################################
######################### Pathological stage plots #############################
################################################################################

## Pathological stage (merge) x TMB
ggplot(surv.dt, aes(x=pathological_stage_refactor, y=TMB_per_region.genome_TMB, color = pathological_stage_refactor)) +
  geom_boxplot(width = 0.25, color = "black") +
  geom_point() +
  geom_text(data=surv.dt[pathological_stage == "4A"], label = "IV", hjust=-1.5, color = "black") +
  theme_classic() +
  ylab("Whole genome TMB") +
  xlab("Pathological stage") +
  scale_color_discrete(name = "Pathological stage")

ggsave("Results/Misc plots/TMBxpath_stage.png", width=7, height=7, dpi=300)

## Pathological stage (all) x TMB X Age
ggplot(surv.dt, aes(x=TMB_per_region.genome_TMB, y=age, color=pathological_stage)) +
  geom_point(size=2) +
  theme_classic() +
  scale_y_continuous(breaks = seq(25, 85,by = 5), limits = c(25, 85)) + 
  xlab("Whole genome TMB") +
  ylab("Age") +
  scale_color_discrete(name = "Pathological stage")


ggsave("Results/Misc plots/TMB x Path stage x Age/path_stage_all.png", width=7, height=7, dpi=300)

## Pathological stage (merge) x TMB X Age
ggplot(surv.dt, aes(x=TMB_per_region.genome_TMB, y=age, color=pathological_stage_refactor)) +
  geom_point(size=2) +
  theme_classic() +
  scale_y_continuous(breaks = seq(25, 85,by = 5), limits = c(25, 85)) + 
  xlab("Whole genome TMB") +
  ylab("Age") +
  scale_color_discrete(name = "Pathological stage")

ggsave("Results/Misc plots/TMB x Path stage x Age/path_stage_merge.png", width=7, height=7, dpi=300)

## Pathological stage (individual) x TMB X Age
# Path stage I
ggplot(surv.dt[pathological_stage_refactor == "I"], aes(x=TMB_per_region.genome_TMB, y=age)) +
  geom_point(size=2, color = "#F8766D") +
  theme_classic() +
  scale_y_continuous(breaks = seq(25, 85,by = 5), limits = c(25, 85)) + 
  xlab("Whole genome TMB") +
  ylab("Age") +
  scale_color_discrete(name = "Pathological stage")

ggsave("Results/Misc plots/TMB x Path stage x Age/path_stage_I.png", width=7, height=7, dpi=300)

# Path stage II
ggplot(surv.dt[pathological_stage_refactor == "II"], aes(x=TMB_per_region.genome_TMB, y=age, color = pathological_stage_refactor)) +
  geom_point(size=2, color = "#00BA38") +
  theme_classic() +
  scale_y_continuous(breaks = seq(25, 85,by = 5), limits = c(25, 85)) + 
  xlab("Whole genome TMB") +
  ylab("Age") +
  scale_color_discrete(name = "Pathological stage")

ggsave("Results/Misc plots/TMB x Path stage x Age/path_stage_II.png", width=7, height=7, dpi=300)

# Path stage III & IV
ggplot(surv.dt[pathological_stage_refactor == "III & IV"], aes(x=TMB_per_region.genome_TMB, y=age, color = pathological_stage_refactor)) +
  geom_point(size=2, color = "#619CFF") +
  theme_classic() +
  scale_y_continuous(breaks = seq(25, 85,by = 5), limits = c(25, 85)) + 
  xlab("Whole genome TMB") +
  ylab("Age") +
  scale_color_discrete(name = "Pathological stage")

ggsave("Results/Misc plots/TMB x Path stage x Age/path_stage_III_IV.png", width=7, height=7, dpi=300)

# All path stage and passive smoking status
ggplot(data=surv.dt, aes(x=factor(passive.smoking), y=complete_WGS_TMB, color = pathological_stage_refactor)) + 
  geom_point() +
  scale_x_discrete(breaks = c(0,1)) +
  theme_classic() +
  theme(axis.text.x = element_text()) +
  ylab("Whole genome TMB") +
  xlab("Passive smoking") +
  scale_color_discrete(name = "Pathological stage")

ggsave("Results/Misc plots/passive_smoking_path_stage.png", width=5, height=7, dpi=300)


################################################################################
########################### SBS signature plots ################################
################################################################################

SBS <- fread("Data/COSMIC_SBS96_Activities.txt", sep = "\t", header = TRUE)
SBS[, Samples:=SBS[, gsub("_","-",as.character(Samples))]]
SBS.df <- SBS[Samples %in% clin_data[, Patient_ID]]
SBS.df <- merge(SBS.df, clin_data[, .(Patient_ID, age)], by.x="Samples", by.y="Patient_ID", )
SBS.df <- as.data.table(reshape2::melt(SBS.df, id.vars=c("Samples", "age"), variable.name = "SBS_type", value.name = "Count"))

# Order by age
ggplot(SBS.df, aes(x=reorder(gsub("NSLC-","",as.character(SBS.df[, Samples])), -age), y=Count, fill=SBS_type)) +
  geom_bar(position="stack", stat="identity", width = 0.9) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10), axis.ticks.x=element_blank(),
        panel.grid.major.y = element_line(colour="black", linewidth=0.5), legend.text=element_text(size=10), axis.title = element_text(size=12), 
        panel.spacing = unit(0.5, "lines"), strip.placement = "outside", strip.text.x = element_text(angle = 90, size = 10), 
        plot.margin=unit(c(5.5,5.5,35,5.5),"pt"), panel.grid.minor = element_blank()) +
  scale_y_continuous(name="Mutation count",  limits = c(0, 30000), expand = c(0, 0)) +
  scale_x_discrete(name="Patients") +
  scale_fill_ucscgb(name="Signature")

ggsave("Results/Misc plots/SBS_ordered_age.png", width=8, height=8, dpi=300)

################################################################################
############### Variants plotted according to their CADD score #################
################################################################################

# CADD score info: https://cadd.gs.washington.edu/info

# Top 10% PHRED scaled CADD score cutoff for hotspot regions
hotspot_cutoff.cadd <- 10

# Arranging chromosome names for ggplot
chrom_name <- c(1:22, 'X', 'Y')
VEP_data_all_patients$CHROM <- factor(VEP_data_all_patients$CHROM, levels=chrom_name)

# Plot
ggplot(VEP_data_all_patients, 
       aes(x=POS, y=CADD_PHRED, color = CHROM)) +
  geom_line(linewidth = 0.2) +
  geom_hline(yintercept = hotspot_cutoff.cadd, color = "red") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.placement='outside', # for facet_grid label position
        panel.spacing = unit(0.3, "lines"),
        legend.position = "none") + # for facet_grid space between panels
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,5)) +
  xlab("Chromosomes") + 
  ylab("PHRED CADD score") +
  facet_grid(~CHROM,
             space = "free_x",
             scales = "free_x",
             switch = "x")

ggsave(file=glue("Results/Misc plots/CADD_score_frequency.png"), width=15, height=7, dpi=300)

################################################################################
############### Variants plotted according to their frequency ##################
################################################################################

############ With the frequency pre calculated by sliding windows ##############

# Window size used to compute sliding windows
window_size <- 5000
Mb_length <-  window_size/1e6

## Merging all chromosome sliding windows densities together
chrom_name <- c(1:22, 'X', 'Y')
# Reading frequency data of chromosomes with a specific window size
chrom_frequency <- rbindlist(lapply(chrom_name, 
                                    function(x) fread(glue("Data/Sliding windows and frequency/chromosome_{x}_sw_frequency_{Mb_length}Mb.csv"))[, chrom:=glue("{x}")]),
                             use.names = TRUE)
# Arranging chromosome names for ggplot
chrom_frequency$chrom <- factor(chrom_frequency$chrom, levels=chrom_name)

# 3rd quartile cutoff for hotspot regions
hotspot_cutoff.third_quartile <- quantile(chrom_frequency$variant_frequency)[4]
hotspot_cutoff.cadd <- quantile(chrom_frequency$CADD_PHRED)[4]

# Plotting the frequency
ggplot(data = chrom_frequency, aes(x = window_start, y = variant_frequency, color = chrom)) + 
  geom_line(linewidth = 0.2) +
  geom_area() +
  geom_hline(yintercept = hotspot_cutoff, color = "red") + 
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.placement='outside', # for facet_grid label position
        panel.spacing = unit(0.3, "lines"),
        legend.position = "none") + # for facet_grid space between panels
  scale_color_hue() +
  xlab("Chromosomes") + 
  ylab("Variant frequency") +
  facet_grid(~chrom,
             space = "free_x",
             scales = "free_x",
             switch = "x")

ggsave(file=glue("Results/Frequency plots/sw_{Mb_length}Mb_frequency.png"), width=15, height=7, dpi=300)


## Merging all chromosome sliding windows densities together and combining all sliding windows sizes
chrom_name <- c(1:22, 'X', 'Y')
window_size.all <- list(0.005, 0.05, 0.5, 1, 10)

chrom_frequency.all <- rbindlist(lapply(chrom_name, 
                                        function(x) rbindlist(lapply(window_size.all, 
                                                                     function(i) fread(glue("Data/Sliding windows and frequency/chromosome_{x}_sw_frequency_{i}Mb.csv"))[, ':='(chrom=glue("{x}"), window=glue("{i}Mb"))]))))
chrom_frequency.all$chrom <- factor(chrom_frequency.all$chrom, levels=chrom_name)

for(size in window_size.all){
  chrom_frequency.all[window == glue("{size}Mb"), variant_frequency_norm:=((variant_frequency-min(variant_frequency))/(max(variant_frequency)-min(variant_frequency)))]
}

# One combined plot for all normalized window size
frequency.all.gp <- ggplot(data = chrom_frequency.all, aes(x = window_start, y = variant_frequency_norm, fill = window, color = window)) + 
  geom_line(linewidth=0.1) +
  #geom_hline(yintercept = hotspot_cutoff, color = "red") + 
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.placement='outside', # for facet_grid label position
        panel.spacing = unit(0.3, "lines")) + # for facet_grid space between panels
  scale_fill_hue(name = "Sliding window size") +
  scale_color_hue(name = "Sliding window size") +
  xlab("Chromosomes") + 
  ylab("Variant frequency") +
  facet_grid(~chrom,
             space = "free_x",
             scales = "free_x",
             switch = "x")

ggsave(frequency.all.gp, file=glue("Results/frequency plots/all_sw_frequency.png"), width=15, height=7, dpi=300)


################################################################################
################################ RpFS graphs ###################################
################################################################################

# Correlation between WGS TMB and RFS
ggplot(na.omit(surv.dt, cols = "time_RpFS"), 
       aes(y = complete_WGS_TMB, x = time_RpFS, color = factor(RpFS_indicator))) +
  geom_point(size=2) +
  theme_classic() +
  ylab("Whole-genome TMB") +
  xlab("Relapse-free survival") +
  scale_color_discrete(name = "Relapse") +
  geom_vline(xintercept = 730, lty = "dashed") +
  geom_vline(xintercept = 1825, lty = "dashed") +
  annotate("text", x=730, y = 10, label = "2 years", hjust = -0.1) + 
  annotate("text", x=1825, y = 10, label = "5 years", hjust = -0.1)

ggsave(file=glue("Results/RFS plots/RFS_WGS_correlation.png"), width=7, height=7, dpi=300)

# Correlation between WES TMB and RFS
ggplot(na.omit(surv.dt, cols = "time_RpFS"), 
       aes(y = complete_WES_TMB, x = time_RpFS, color = factor(RpFS_indicator))) +
  geom_point(size=2) +
  theme_classic() +
  ylab("Whole-exome TMB") +
  xlab("Relapse-free survival") +
  scale_color_discrete(name = "Relapse") +
  geom_vline(xintercept = 730, lty = "dashed") +
  geom_vline(xintercept = 1825, lty = "dashed") +
  annotate("text", x=730, y = 10, label = "2 years", hjust = -0.1) + 
  annotate("text", x=1825, y = 10, label = "5 years", hjust = -0.1)

ggsave(file=glue("Results/RFS plots/RFS_WES_correlation.png"), width=7, height=7, dpi=300)

## 2 years RFS with path stage annotation
ggplot(na.omit(surv.dt, cols="recurrence.two"), aes(x=complete_WGS_TMB)) +
  geom_density() +
  lims(x = c(0,10), y=c(0,0.6)) +
  facet_grid(~recurrence.two)
shapiro.test(surv.dt[recurrence.two == "<2", complete_WGS_TMB])
shapiro.test(surv.dt[recurrence.two == ">=2", complete_WGS_TMB])

ggplot(na.omit(surv.dt, cols="recurrence.two"), aes(x=complete_WES_TMB)) +
  geom_density() +
  lims(x = c(0,8), y =c(0, 0.7)) +
  facet_grid(~recurrence.two)
shapiro.test(surv.dt[recurrence.two == "<2", complete_WES_TMB])
shapiro.test(surv.dt[recurrence.two == ">=2", complete_WES_TMB])

# 2 years RFS x continuous WGS TMB and path stage
WGS.two.WRS <- wilcox.test(complete_WGS_TMB ~ recurrence.two, data = surv.dt, # Comparing the distribution between the two RFS survival subgroups
                           na.rm = TRUE, paired = FALSE, exact = FALSE, conf.int = TRUE)
WGS.two.WRS.p <- round(WGS.two.WRS$p.value, digits = 3)

ggplot(na.omit(surv.dt, cols="recurrence.two"), aes(x=recurrence.two, y=complete_WGS_TMB)) +
  geom_boxplot(width=0.5, outlier.colour = "white") +
  geom_jitter(aes(color = pathological_stage_refactor), 
              size = 2, width = 0.15) +
  geom_hline(yintercept = 1.7, color = "red") +
  geom_vline(xintercept = 1.5, lty = "dashed") +
  theme_classic() +
  xlab("Relapse free survival") +
  ylab("Whole-genome TMB") + 
  scale_x_discrete(labels = c("<2" = "<2 years\nn=20",
                              ">=2" = "\u22652 years\nn=60")) +
  scale_y_continuous(breaks = c(seq(0,10,1), 1.7),
                     limits = c(0,10), expand = c(0,0)) +
  scale_color_discrete(name = "Pathological stage") +
  theme(axis.text.y = element_text(color = c(rep("black", 11), "red")),
        axis.ticks.y = element_line(color = c(rep("black", 11), "red"))) +
  annotate("text", x=0.5, y = 1.55, label = "n=8", hjust = "left") +
  annotate("text", x=0.5, y = 1.90, label = "n=12", hjust = "left") +
  annotate("text", x=2.4, y = 1.55, label = "n=42", hjust = "left") +
  annotate("text", x=2.4, y = 1.90, label = "n=18", hjust = "left") +
  annotate("text", x=0, y = 10, label = glue("Wilcoxon Rank Sum p={WGS.two.WRS.p}"), hjust = -0.05) +
  coord_cartesian(clip = "off")

ggsave(file=glue("Results/RFS plots/RFS_WGS_2years_path.png"), width=7, height=7, dpi=300)

# 2 years RFS x continuous WES TMB and path stage
WES.two.WRS <- wilcox.test(complete_WES_TMB ~ recurrence.two, data = surv.dt, # Comparing the distribution between the two RFS survival subgroups
                           na.rm = TRUE, paired = FALSE, exact = FALSE, conf.int = TRUE)
WES.two.WRS.p <- round(WES.two.WRS$p.value, digits = 3)

ggplot(na.omit(surv.dt, cols="recurrence.two"), aes(x=recurrence.two, y=complete_WES_TMB)) +
  geom_boxplot(width=0.5, outlier.colour = "white") +
  geom_jitter(aes(color = pathological_stage_refactor), 
              size = 2, width = 0.15) +
  geom_hline(yintercept = 1.2, color = "red") +
  geom_vline(xintercept = 1.5, lty = "dashed") +
  theme_classic() +
  xlab("Relapse free survival") +
  ylab("Whole-exome TMB") +
  scale_x_discrete(labels = c("<2" = "<2 years\nn=20",
                              ">=2" = "\u22652 years\nn=60")) +
  scale_y_continuous(breaks = c(seq(0,8,1), 1.2),
                     limits = c(0,8), expand = c(0,0)) +
  scale_color_discrete(name = "Pathological stage") +
  theme(axis.text.y = element_text(color = c(rep("black", 9), "red")),
        axis.ticks.y = element_line(color = c(rep("black", 9), "red"))) +
  annotate("text", x=0.5, y = 1.05, label = "n=9", hjust = "left") +
  annotate("text", x=0.5, y = 1.40, label = "n=11", hjust = "left") +
  annotate("text", x=2.4, y = 1.05, label = "n=41", hjust = "left") +
  annotate("text", x=2.4, y = 1.40, label = "n=19", hjust = "left") +
  annotate("text", x=0, y = 8, label = glue("Wilcoxon Rank Sum p={WES.two.WRS.p}"), hjust = -0.05) +
  coord_cartesian(clip = "off")

ggsave(file=glue("Results/RFS plots/RFS_WES_2years_path.png"), width=7, height=7, dpi=300)

## 5 years RFS with path stage annotation
ggplot(na.omit(surv.dt, cols="recurrence.five"), aes(x=complete_WGS_TMB)) +
  geom_density() +
  lims(x = c(0,10), y =c(0, 0.6)) +
  facet_grid(~recurrence.five)
shapiro.test(surv.dt[recurrence.five == "<5", complete_WGS_TMB])
shapiro.test(surv.dt[recurrence.five == ">=5", complete_WGS_TMB])

ggplot(na.omit(surv.dt, cols="recurrence.five"), aes(x=complete_WES_TMB)) +
  geom_density() +
  lims(x = c(0,8), y =c(0, 0.7)) +
  facet_grid(~recurrence.five)
shapiro.test(surv.dt[recurrence.five == "<5", complete_WES_TMB])
shapiro.test(surv.dt[recurrence.five == ">=5", complete_WES_TMB])

# 5 years RFS x continuous WGS TMB and path stage
WGS.five.WRS <- wilcox.test(complete_WGS_TMB ~ recurrence.five, data = surv.dt, # Comparing the distribution between the two RFS survival subgroups
                            na.rm = TRUE, paired = FALSE, exact = FALSE, conf.int = TRUE)
WGS.five.WRS.p <- round(WGS.five.WRS$p.value, digits = 3)

ggplot(na.omit(surv.dt, cols="recurrence.five"), aes(x=recurrence.five, y=complete_WGS_TMB)) +
  geom_boxplot(width=0.5, outlier.colour = "white") +
  geom_jitter(aes(color = pathological_stage_refactor), 
              size = 2, width = 0.15) +
  geom_hline(yintercept = 1.7, color = "red") +
  geom_vline(xintercept = 1.5, lty = "dashed") +
  theme_classic() +
  xlab("Relapse free survival") +
  ylab("Whole-genome TMB") + 
  scale_x_discrete(labels = c("<5" = "<5 years\nn=48",
                              ">=5" = "\u22655 years\nn=32")) +
  scale_y_continuous(breaks = c(seq(0,10,1), 1.7),
                     limits = c(0,10), expand = c(0,0)) +
  scale_color_discrete(name = "Pathological stage") +
  theme(axis.text.y = element_text(color = c(rep("black", 11), "red")),
        axis.ticks.y = element_line(color = c(rep("black", 11), "red"))) +
  annotate("text", x=0.5, y = 1.55, label = "n=28", hjust = "left") +
  annotate("text", x=0.5, y = 1.90, label = "n=20", hjust = "left") +
  annotate("text", x=2.4, y = 1.55, label = "n=22", hjust = "left") +
  annotate("text", x=2.4, y = 1.90, label = "n=10", hjust = "left") +
  annotate("text", x=0, y = 10, label = glue("Wilcoxon Rank Sum p={WGS.five.WRS.p}"), hjust = -0.05) +
  coord_cartesian(clip = "off")

ggsave(file=glue("Results/RFS plots/RFS_WGS_5years_path.png"), width=7, height=7, dpi=300)

# 5 years RFS x continuous WES TMB and path stage
WES.five.WRS <- wilcox.test(complete_WES_TMB ~ recurrence.five, data = surv.dt, # Comparing the distribution between the two RFS survival subgroups
                            na.rm = TRUE, paired = FALSE, exact = FALSE, conf.int = TRUE)
WES.five.WRS.p <- round(WES.five.WRS$p.value, digits = 3)

ggplot(na.omit(surv.dt, cols="recurrence.five"), aes(x=recurrence.five, y=complete_WES_TMB)) +
  geom_boxplot(width=0.5, outlier.colour = "white") +
  geom_jitter(aes(color = pathological_stage_refactor), 
              size = 2, width = 0.15) +
  geom_hline(yintercept = 1.2, color = "red") +
  geom_vline(xintercept = 1.5, lty = "dashed") +
  theme_classic() +
  xlab("Relapse free survival") +
  ylab("Whole-exome TMB") +
  scale_x_discrete(labels = c("<5" = "<5 years\nn=48",
                              ">=5" = "\u22655 years\nn=32")) +
  scale_y_continuous(breaks = c(seq(0,8,1), 1.2),
                     limits = c(0,8), expand = c(0,0)) +
  scale_color_discrete(name = "Pathological stage") +
  theme(axis.text.y = element_text(color = c(rep("black", 9), "red")),
        axis.ticks.y = element_line(color = c(rep("black", 9), "red"))) +
  annotate("text", x=0.5, y = 1.05, label = "n=28", hjust = "left") +
  annotate("text", x=0.5, y = 1.40, label = "n=20", hjust = "left") +
  annotate("text", x=2.4, y = 1.05, label = "n=22", hjust = "left") +
  annotate("text", x=2.4, y = 1.40, label = "n=10", hjust = "left") +
  annotate("text", x=0, y = 8, label = glue("Wilcoxon Rank Sum p={WES.five.WRS.p}"), hjust = -0.05) +
  coord_cartesian(clip = "off")

ggsave(file=glue("Results/RFS plots/RFS_WES_5years_path.png"), width=7, height=7, dpi=300)

## 2 years RFS with relapse annotation
# 2 years RFS x continuous WGS TMB and event
WGS.two.WRS.event <- wilcox.test(complete_WGS_TMB ~ recurrence.two, data = surv.dt[RpFS_indicator == 1],
                                 na.rm = TRUE, paired = FALSE, exact = TRUE, conf.int = TRUE) # Comparing the event (replase or not) distribution between the two RFS survival subgroups
WGS.two.WRS.event.p <- round(WGS.two.WRS.event$p.value, digits = 3)

ggplot(na.omit(surv.dt, cols="recurrence.two"), aes(x=recurrence.two, y=complete_WGS_TMB)) +
  geom_boxplot(width=0.5, outlier.colour = "white") +
  geom_jitter(aes(color = factor(RpFS_indicator)), 
              size = 2, width = 0.15) +
  geom_hline(yintercept = 1.7, color = "red") +
  geom_vline(xintercept = 1.5, lty = "dashed") +
  theme_classic() +
  xlab("Relapse free survival") +
  ylab("Whole-genome TMB") + 
  scale_x_discrete(labels = c("<2" = "<2 years\nn=20",
                              ">=2" = "\u22652 years\nn=60")) +
  scale_y_continuous(breaks = c(seq(0,10,1), 1.7),
                     limits = c(0,10), expand = c(0,0)) +
  scale_color_discrete(name = "Relapse") +
  theme(axis.text.y = element_text(color = c(rep("black", 11), "red")),
        axis.ticks.y = element_line(color = c(rep("black", 11), "red"))) +
  annotate("text", x=0.5, y = 1.55, label = "n=8", hjust = "left") +
  annotate("text", x=0.5, y = 1.90, label = "n=12", hjust = "left") +
  annotate("text", x=2.4, y = 1.55, label = "n=42", hjust = "left") +
  annotate("text", x=2.4, y = 1.90, label = "n=18", hjust = "left") +
  annotate("text", x=0, y = 10, label = glue("Wilcoxon Rank Sum p={WGS.two.WRS.p}"), hjust = -0.05) +
  annotate("text", x=0, y = 9.5, label = glue("Event Wilcoxon Rank Sum p={WGS.two.WRS.event.p}"), colour = "#00b0c4", hjust = -0.05) +
  coord_cartesian(clip = "off")

ggsave(file=glue("Results/RFS plots/RFS_WGS_2years_event.png"), width=7, height=7, dpi=300)

# 2 years RFS x continuous WES TMB and event
WES.two.WRS.event <- wilcox.test(complete_WES_TMB ~ recurrence.two, data = surv.dt[RpFS_indicator == 1],
                                 na.rm = TRUE, paired = FALSE, exact = TRUE, conf.int = TRUE) # Comparing the event (relapse or not) distribution between the two RFS survival subgroups
WES.two.WRS.event.p <- round(WES.two.WRS.event$p.value, digits = 3)

ggplot(na.omit(surv.dt, cols="recurrence.two"), aes(x=recurrence.two, y=complete_WES_TMB)) +
  geom_boxplot(width=0.5, outlier.colour = "white") +
  geom_jitter(aes(color = factor(RpFS_indicator)), 
              size = 2, width = 0.15) +
  geom_hline(yintercept = 1.2, color = "red") +
  geom_vline(xintercept = 1.5, lty = "dashed") +
  theme_classic() +
  xlab("Relapse free survival") +
  ylab("Whole-exome TMB") +
  scale_x_discrete(labels = c("<2" = "<2 years\nn=20",
                              ">=2" = "\u22652 years\nn=60")) +
  scale_y_continuous(breaks = c(seq(0,8,1), 1.2),
                     limits = c(0,8), expand = c(0,0)) +
  scale_color_discrete(name = "Relapse") +
  theme(axis.text.y = element_text(color = c(rep("black", 9), "red")),
        axis.ticks.y = element_line(color = c(rep("black", 9), "red"))) +
  annotate("text", x=0.5, y = 1.05, label = "n=9", hjust = "left") +
  annotate("text", x=0.5, y = 1.40, label = "n=11", hjust = "left") +
  annotate("text", x=2.4, y = 1.05, label = "n=41", hjust = "left") +
  annotate("text", x=2.4, y = 1.40, label = "n=19", hjust = "left") +
  annotate("text", x=0, y = 8, label = glue("Wilcoxon Rank Sum p={WES.two.WRS.p}"), hjust = -0.05) +
  annotate("text", x=0, y = 7.6, label = glue("Event Wilcoxon Rank Sum p={WES.two.WRS.event.p}"), colour = "#00b0c4", hjust = -0.05) +
  coord_cartesian(clip = "off")

ggsave(file=glue("Results/RFS plots/RFS_WES_2years_event.png"), width=7, height=7, dpi=300)

## 5 years RFS with path stage annotation
# 5 years RFS x continuous WGS TMB and event
WGS.five.WRS.event <- wilcox.test(complete_WGS_TMB ~ recurrence.five, data = surv.dt[RpFS_indicator == 1],
                                  na.rm = TRUE, paired = FALSE, exact = TRUE, conf.int = TRUE) # Comparing the event (relapse or not) distribution between the two RFS survival subgroups
WGS.five.WRS.event.p <- round(WGS.five.WRS.event$p.value, digits = 3)

ggplot(na.omit(surv.dt, cols="recurrence.five"), aes(x=recurrence.five, y=complete_WGS_TMB)) +
  geom_boxplot(width=0.5, outlier.colour = "white") +
  geom_jitter(aes(color = factor(RpFS_indicator)), 
              size = 2, width = 0.15) +
  geom_hline(yintercept = 1.7, color = "red") +
  geom_vline(xintercept = 1.5, lty = "dashed") +
  theme_classic() +
  xlab("Relapse free survival") +
  ylab("Whole-genome TMB") + 
  scale_x_discrete(labels = c("<5" = "<5 years\nn=48",
                              ">=5" = "\u22655 years\nn=32")) +
  scale_y_continuous(breaks = c(seq(0,10,1), 1.7),
                     limits = c(0,10), expand = c(0,0)) +
  scale_color_discrete(name = "Relapse") +
  theme(axis.text.y = element_text(color = c(rep("black", 11), "red")),
        axis.ticks.y = element_line(color = c(rep("black", 11), "red"))) +
  annotate("text", x=0.5, y = 1.55, label = "n=28", hjust = "left") +
  annotate("text", x=0.5, y = 1.90, label = "n=20", hjust = "left") +
  annotate("text", x=2.4, y = 1.55, label = "n=22", hjust = "left") +
  annotate("text", x=2.4, y = 1.90, label = "n=10", hjust = "left") +
  annotate("text", x=0, y = 10, label = glue("Wilcoxon Rank Sum p={WGS.five.WRS.p}"), hjust = -0.05) +
  #annotate("text", x=0, y = 9.5, label = glue("Event Wilcoxon Rank Sum p={WGS.five.WRS.event.p}"), colour = "#00b0c4", hjust = -0.05) +
  coord_cartesian(clip = "off")

ggsave(file=glue("Results/RFS plots/RFS_WGS_5years_event.png"), width=7, height=7, dpi=300)

# 5 years RFS x continuous WES TMB and event
WES.five.WRS.event <- wilcox.test(complete_WES_TMB ~ recurrence.five, data = surv.dt[RpFS_indicator == 1], # Comparing the distribution between the two RFS survival subgroups
                                  na.rm = TRUE, paired = FALSE, exact = TRUE, conf.int = TRUE)
WES.five.WRS.event.p <- round(WES.five.WRS.event$p.value, digits = 3)

ggplot(na.omit(surv.dt, cols="recurrence.five"), aes(x=recurrence.five, y=complete_WES_TMB)) +
  geom_boxplot(width=0.5, outlier.colour = "white") +
  geom_jitter(aes(color = factor(RpFS_indicator)), 
              size = 2, width = 0.15) +
  geom_hline(yintercept = 1.2, color = "red") +
  geom_vline(xintercept = 1.5, lty = "dashed") +
  theme_classic() +
  xlab("Relapse free survival") +
  ylab("Whole-exome TMB") +
  scale_x_discrete(labels = c("<5" = "<5 years\nn=48",
                              ">=5" = "\u22655 years\nn=32")) +
  scale_y_continuous(breaks = c(seq(0,8,1), 1.2),
                     limits = c(0,8), expand = c(0,0)) +
  scale_color_discrete(name = "Relapse") +
  theme(axis.text.y = element_text(color = c(rep("black", 9), "red")),
        axis.ticks.y = element_line(color = c(rep("black", 9), "red"))) +
  annotate("text", x=0.5, y = 1.05, label = "n=28", hjust = "left") +
  annotate("text", x=0.5, y = 1.40, label = "n=20", hjust = "left") +
  annotate("text", x=2.4, y = 1.05, label = "n=22", hjust = "left") +
  annotate("text", x=2.4, y = 1.40, label = "n=10", hjust = "left") +
  annotate("text", x=0, y = 8, label = glue("Wilcoxon Rank Sum p={WES.five.WRS.p}"), hjust = -0.05) +
  #annotate("text", x=0, y = 7.6, label = glue("Event Wilcoxon Rank Sum p={WES.five.WRS.event.p}"), colour = "#00b0c4", hjust = -0.05) +
  coord_cartesian(clip = "off")

ggsave(file=glue("Results/RFS plots/RFS_WES_5years_event.png"), width=7, height=7, dpi=300)


################################################################################
################################# OS graphs ####################################
################################################################################

## 5 years OS with path stage annotation
ggplot(na.omit(surv.dt, cols="os.five"), aes(x=complete_WGS_TMB)) +
  geom_density() +
  lims(x = c(0,10), y =c(0, 0.6)) +
  facet_grid(~os.five)
shapiro.test(surv.dt[os.five == "<5", complete_WGS_TMB])
shapiro.test(surv.dt[os.five == ">=5", complete_WGS_TMB])

ggplot(na.omit(surv.dt, cols="os.five"), aes(x=complete_WES_TMB)) +
  geom_density() +
  lims(x = c(0,8), y =c(0, 0.7)) +
  facet_grid(~os.five)
shapiro.test(surv.dt[os.five == "<5", complete_WES_TMB])
shapiro.test(surv.dt[os.five == ">=5", complete_WES_TMB])

# 5 years OS x continuous WGS TMB and path stage
WGS.os.five.WRS <- wilcox.test(complete_WGS_TMB ~ os.five, data = surv.dt, # Comparing the distribution between the two OS survival subgroups
                            na.rm = TRUE, paired = FALSE, exact = FALSE, conf.int = TRUE)
WGS.os.five.WRS.p <- round(WGS.os.five.WRS$p.value, digits = 3)

ggplot(na.omit(surv.dt, cols="os.five"), aes(x=os.five, y=complete_WGS_TMB)) +
  geom_boxplot(width=0.5, outlier.colour = "white") +
  geom_jitter(aes(color = pathological_stage_refactor), 
              size = 2, width = 0.15) +
  geom_hline(yintercept = 1.7, color = "red") +
  geom_vline(xintercept = 1.5, lty = "dashed") +
  theme_classic() +
  xlab("Overall survival") +
  ylab("Whole-genome TMB") + 
  scale_x_discrete(labels = c("<5" = "<5 years\nn=33",
                              ">=5" = "\u22655 years\nn=49")) +
  scale_y_continuous(breaks = c(seq(0,10,1), 1.7),
                     limits = c(0,10), expand = c(0,0)) +
  scale_color_discrete(name = "Pathological stage") +
  theme(axis.text.y = element_text(color = c(rep("black", 11), "red")),
        axis.ticks.y = element_line(color = c(rep("black", 11), "red"))) +
  annotate("text", x=0.5, y = 1.55, label = "n=18", hjust = "left") +
  annotate("text", x=0.5, y = 1.90, label = "n=15", hjust = "left") +
  annotate("text", x=2.4, y = 1.55, label = "n=32", hjust = "left") +
  annotate("text", x=2.4, y = 1.90, label = "n=17", hjust = "left") +
  annotate("text", x=0, y = 10, label = glue("Wilcoxon Rank Sum p={WGS.os.five.WRS.p}"), hjust = -0.05) +
  coord_cartesian(clip = "off")

ggsave(file=glue("Results/OS plots/OS_WGS_5years_path.png"), width=7, height=7, dpi=300)

# 5 years OS x continuous WES TMB and path stage
WES.os.five.WRS <- wilcox.test(complete_WES_TMB ~ os.five, data = surv.dt, # Comparing the distribution between the two OS survival subgroups
                            na.rm = TRUE, paired = FALSE, exact = FALSE, conf.int = TRUE)
WES.os.five.WRS.p <- round(WES.os.five.WRS$p.value, digits = 3)

ggplot(na.omit(surv.dt, cols="os.five"), aes(x=os.five, y=complete_WES_TMB)) +
  geom_boxplot(width=0.5, outlier.colour = "white") +
  geom_jitter(aes(color = pathological_stage_refactor), 
              size = 2, width = 0.15) +
  geom_hline(yintercept = 1.2, color = "red") +
  geom_vline(xintercept = 1.5, lty = "dashed") +
  theme_classic() +
  xlab("Overall survival") +
  ylab("Whole-exome TMB") +
  scale_x_discrete(labels = c("<5" = "<5 years\nn=33",
                              ">=5" = "\u22655 years\nn=49")) +
  scale_y_continuous(breaks = c(seq(0,8,1), 1.2),
                     limits = c(0,8), expand = c(0,0)) +
  scale_color_discrete(name = "Pathological stage") +
  theme(axis.text.y = element_text(color = c(rep("black", 9), "red")),
        axis.ticks.y = element_line(color = c(rep("black", 9), "red"))) +
  annotate("text", x=0.5, y = 1.05, label = "n=18", hjust = "left") +
  annotate("text", x=0.5, y = 1.40, label = "n=15", hjust = "left") +
  annotate("text", x=2.4, y = 1.05, label = "n=32", hjust = "left") +
  annotate("text", x=2.4, y = 1.40, label = "n=17", hjust = "left") +
  annotate("text", x=0, y = 8, label = glue("Wilcoxon Rank Sum p={WES.os.five.WRS.p}"), hjust = -0.05) +
  coord_cartesian(clip = "off")

ggsave(file=glue("Results/OS plots/OS_WES_5years_path.png"), width=7, height=7, dpi=300)

## 5 years OS with path stage annotation
# 5 years OS x continuous WGS TMB and event
WGS.os.five.WRS.event <- wilcox.test(complete_WGS_TMB ~ os.five, data = surv.dt[RpFS_indicator == 1],
                                  na.rm = TRUE, paired = FALSE, exact = TRUE, conf.int = TRUE) # Comparing the event (relapse or not) distribution between the two OS survival subgroups
WGS.os.five.WRS.event.p <- round(WGS.os.five.WRS.event$p.value, digits = 3)

ggplot(na.omit(surv.dt, cols="os.five"), aes(x=os.five, y=complete_WGS_TMB)) +
  geom_boxplot(width=0.5, outlier.colour = "white") +
  geom_jitter(aes(color = factor(VitalStatus)), 
              size = 2, width = 0.15) +
  geom_hline(yintercept = 1.7, color = "red") +
  geom_vline(xintercept = 1.5, lty = "dashed") +
  theme_classic() +
  xlab("Overall survival") +
  ylab("Whole-genome TMB") + 
  scale_x_discrete(labels = c("<5" = "<5 years\nn=33",
                              ">=5" = "\u22655 years\nn=49")) +
  scale_y_continuous(breaks = c(seq(0,10,1), 1.7),
                     limits = c(0,10), expand = c(0,0)) +
  scale_color_discrete(name = "Death") +
  theme(axis.text.y = element_text(color = c(rep("black", 11), "red")),
        axis.ticks.y = element_line(color = c(rep("black", 11), "red"))) +
  annotate("text", x=0.5, y = 1.55, label = "n=18", hjust = "left") +
  annotate("text", x=0.5, y = 1.90, label = "n=15", hjust = "left") +
  annotate("text", x=2.4, y = 1.55, label = "n=32", hjust = "left") +
  annotate("text", x=2.4, y = 1.90, label = "n=17", hjust = "left") +
  annotate("text", x=0, y = 10, label = glue("Wilcoxon Rank Sum p={WGS.five.WRS.p}"), hjust = -0.05) +
  annotate("text", x=0, y = 9.5, label = glue("Event Wilcoxon Rank Sum p={WGS.os.five.WRS.event.p}"), colour = "#00b0c4", hjust = -0.05) +
  coord_cartesian(clip = "off")

ggsave(file=glue("Results/OS plots/OS_WGS_5years_event.png"), width=7, height=7, dpi=300)

# 5 years OS x continuous WES TMB and event
WES.os.five.WRS.event <- wilcox.test(complete_WES_TMB ~ os.five, data = surv.dt[RpFS_indicator == 1], # Comparing the distribution between the two OS survival subgroups
                                  na.rm = TRUE, paired = FALSE, exact = TRUE, conf.int = TRUE)
WES.os.five.WRS.event.p <- round(WES.os.five.WRS.event$p.value, digits = 3)

ggplot(na.omit(surv.dt, cols="os.five"), aes(x=os.five, y=complete_WES_TMB)) +
  geom_boxplot(width=0.5, outlier.colour = "white") +
  geom_jitter(aes(color = factor(VitalStatus)), 
              size = 2, width = 0.15) +
  geom_hline(yintercept = 1.2, color = "red") +
  geom_vline(xintercept = 1.5, lty = "dashed") +
  theme_classic() +
  xlab("Overall survival") +
  ylab("Whole-exome TMB") +
  scale_x_discrete(labels = c("<5" = "<5 years\nn=33",
                              ">=5" = "\u22655 years\nn=49")) +
  scale_y_continuous(breaks = c(seq(0,8,1), 1.2),
                     limits = c(0,8), expand = c(0,0)) +
  scale_color_discrete(name = "Death") +
  theme(axis.text.y = element_text(color = c(rep("black", 9), "red")),
        axis.ticks.y = element_line(color = c(rep("black", 9), "red"))) +
  annotate("text", x=0.5, y = 1.05, label = "n=18", hjust = "left") +
  annotate("text", x=0.5, y = 1.40, label = "n=15", hjust = "left") +
  annotate("text", x=2.4, y = 1.05, label = "n=32", hjust = "left") +
  annotate("text", x=2.4, y = 1.40, label = "n=17", hjust = "left") +
  annotate("text", x=0, y = 8, label = glue("Wilcoxon Rank Sum p={WES.five.WRS.p}"), hjust = -0.05) +
  annotate("text", x=0, y = 7.6, label = glue("Event Wilcoxon Rank Sum p={WES.os.five.WRS.event.p}"), colour = "#00b0c4", hjust = -0.05) +
  coord_cartesian(clip = "off")

ggsave(file=glue("Results/OS plots/OS_WES_5years_event.png"), width=7, height=7, dpi=300)

################################################################################
############################### Rainfall plots #################################
################################################################################

## karyoploteR mutations plots for the entire cohort. Reference: https://bernatgel.github.io/karyoploter_tutorial/

# # Preemptive rainfall plot to get max distance value for the final rainfall plot
# kp <- plotKaryotype(plot.type=4, ideogram.plotter = NULL, labels.plotter = NULL)
# kpRnfall.all <- kpPlotRainfall(kp, data = somatic.mutations.all)
# max_distance_rounded.all <- ceiling(max(unlist(kpRnfall.all$latest.plot$computed.values$distances))*2)/2
# 
# # Making the rainfall plot
# pp <- getDefaultPlotParams(plot.type = 4)
# pp$data1inmargin <- 0
# pp$bottommargin <- 20
# 
# png(file=glue("Results/Rainfall plots/all_NSLC_patients_rainfall.png"),
#     width=465, height=225, units = "mm", res=300)
# 
# kp <- plotKaryotype(plot.type=4, ideogram.plotter = NULL,
#                     labels.plotter = NULL, plot.params = pp)
# kpAddCytobandsAsLine(kp)
# kpAddChromosomeNames(kp, srt=45)
# variant.colors <- getVariantsColors(somatic.mutations.all$REF, somatic.mutations.all$ALT)
# kpPlotRainfall(kp, data = somatic.mutations.all, col = variant.colors, r0=0.7, r1=0, ymin = 0, ymax = max_distance_rounded.all)
# kpAxis(kp, ymax = 0, ymin = max_distance_rounded.all, tick.pos = floor(max_distance_rounded.all):0, r0=0, r1=0.7)
# kpAddLabels(kp, labels = c("Distance between mutations (log10)"), srt=90, pos=1, label.margin = 0.04, r0=0, r1=0.7)
# kpPlotDensity(kp, data = somatic.mutations.all, r0=0.72, r1=1)
# kpAddLabels(kp, labels = c("Mutation frequency"), srt=90, pos=1, label.margin = 0.04, r0=0.71, r1=1)


################################################################################
################################################################################
### Section for all individual patients variants

# karyoploteR mutations graphs for each patient. Reference: https://bernatgel.github.io/karyoploter_tutorial/
# ** Warning: slow loop, redo this only if needed.
# Loop for making Rainfall plots for each patient.

# for(Patient in Patients_list) {    
#   # Preparing data for rainfall plot 
#   somatic.mutations <- toGRanges(get(glue("VEP_data_{Patient}"))[, list(CHROM, POS, END_POS, ID, REF, ALT, region_type, mutation_type)])
#   seqlevelsStyle(somatic.mutations) <- "UCSC"
#   
#   # Preemptive rainfall plot to get max distance value for the final rainfall plot
#   kp <- plotKaryotype(plot.type=4, ideogram.plotter = NULL, labels.plotter = NULL)
#   kpRnfall <- kpPlotRainfall(kp, data = somatic.mutations)
#   max_distance_rounded <- ceiling(max(unlist(kpRnfall$latest.plot$computed.values$distances))*2)/2
#   
#   # Making the rainfall plot
#   pp <- getDefaultPlotParams(plot.type = 4)
#   pp$data1inmargin <- 0
#   pp$bottommargin <- 20
#   
#   png(file=glue("Results/Rainfall plots/NSLC_{Patient}_rainfall.png"),
#       width=465, height=225, units = "mm", res=300)
#   
#   kp <- plotKaryotype(plot.type=4, ideogram.plotter = NULL,
#                       labels.plotter = NULL, plot.params = pp)
#   kpAddCytobandsAsLine(kp)
#   kpAddChromosomeNames(kp, srt=45)
#   variant.colors <- getVariantsColors(somatic.mutations$REF, somatic.mutations$ALT)
#   kpPlotRainfall(kp, data = somatic.mutations, col = variant.colors, r0=0.7, r1=0, ymin = 0, ymax = max_distance_rounded)
#   kpAxis(kp, ymax = 0, ymin = max_distance_rounded, tick.pos = floor(max_distance_rounded):0, r0=0, r1=0.7)
#   kpAddLabels(kp, labels = c("Distance between mutations (log10)"), srt=90, pos=1, label.margin = 0.04, r0=0, r1=0.7)
#   kpPlotDensity(kp, data = somatic.mutations, r0=0.72, r1=1)
#   kpAddLabels(kp, labels = c("Mutation frequency"), srt=90, pos=1, label.margin = 0.04, r0=0.71, r1=1)
#   
#   dev.off()
# }

################################################################################
################################### ARCHIVES ###################################
################################################################################

# ## Using ChromoMap to show variant frequency on chromosomes.
#
# library(VariantAnnotation)
# library(GenomicFeatures)
# library(AnnotationHub)
# library(chromoMap) # reference: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04556-z
# 
# ah = AnnotationHub()
# #AnnotationHub::query(ah, c('gtf', 'Homo_sapiens', 'GRCh37'))
# GRCh37.gtf<- ah[['AH10684']]
# 
# 
# ## subset the gtf files for only protein_coding genes and lincRNAs
# GRCh37.gtf<- GRCh37.gtf[GRCh37.gtf$gene_biotype %in% c('protein_coding', 'lincRNA')]
# 
# ## make a txdb and keep conventional chromosomes
# GRCh37.txdb <- makeTxDbFromGRanges(GRCh37.gtf)
# GRCh37.txdb <- keepSeqlevels(GRCh37.txdb, c(1:22, 'X', 'Y'), pruning.mode = 'coarse')
# chrom_grngs <- as(seqinfo(GRCh37.txdb), 'GRanges')
# chrom_data <- data.table()
# chrom_data[, CHROM:=seqnames(chrom_grngs)%>%as.character()]
# chrom_data[, POS:=start(chrom_grngs)%>%as.numeric()]
# chrom_data[, END_POS:=end(chrom_grngs)%>%as.numeric()]
# 
# anno_data <- VEP_data_all_patients[, list(ID, CHROM, POS, END_POS)]
# 
# chromo <- chromoMap(list(chrom_data), list(anno_data))
# 
# VEP_data_all_patients[, loci:= chromo$x$chData[[1]]$loci]
# window_count <- as.data.table(table(VEP_data_all_patients[, loci]))
# VEP_data_all_patients <- merge(VEP_data_all_patients, window_count, by.x = "loci", by.y = "V1")
# setnames(VEP_data_all_patients, "N", "loci_var_count")
# VEP_data <- VEP_data %>% mutate(loci_var_cat = case_when(loci_var_count < 5 ~ "[1, 5]",
#                                                          loci_var_count >= 5 & loci_var_count <10 ~ "[5, 10[",
#                                                          loci_var_count >=10 ~ "[10, ["))
# 
# anno_data_2 <- VEP_data_all_patients[, list(ID, CHROM, POS, END_POS, loci_var_count)]
# 
# chromoMap(list(chrom_data), list(anno_data_2),
#           data_based_color_map = T,
#           data_type = "numeric",
#           plots = "scatter",
#           plot_filter = list(c("col","byCategory")),
#           ch2D.colors = c("red3","orange3","purple"),
#           remove.last.window = FALSE)
# 
# chromoMap(list(chrom_data), list(anno_data_2),
#           data_based_color_map = T,
#           data_type = "numeric",
#           plots = "bar",
#           remove.last.window = FALSE)
# 
# 
# window_length <- chromo$x$chData[[1]]$loci_end[1] - chromo$x$chData[[1]]$loci_start[1]
# 
# for(chromosome in unique(VEP_data[, CHROM])) {
#   chr_indices <- c(first(grep(glue("chromap-{chromosome}-"),  chromo$x$chData[[1]]$loci)), last(grep(glue("chromap-{chromosome}-"),  chromo$x$chData[[1]]$loci)))
#   assign(glue("chr_{chromosome}"), chr_indices)
#   assign(glue("chr_{chromosome}_wd_length"), chromo$x$chData[[1]]$loci_end[chr_indices[1]] - chromo$x$chData[[1]]$loci_start[chr_indices[1]])
# }
# 
# for(chromosome in chrom_name) {
# assign(glue("chrom_{chromosome}_frequency"), fread(glue("Data/Sliding windows and frequency/chromosome_{chromosome}_sw_frequency.csv")))
# assign(glue("chrom_{chromosome}_frequency.gr"), makeGRangesFromDataFrame(get(glue("chrom_{chromosome}_frequency"))[, chrom:=glue("chr{chromosome}")], 
#                                                                        keep.extra.columns = TRUE, 
#                                                                        ignore.strand = TRUE, 
#                                                                        seqnames.field = "chrom", 
#                                                                        start.field = "window_start",
#                                                                        end.field = "window_end"))
# }
# 
# ############# With the frequency directly calculated by karyoteR ###############
# 
# # Preparing the GRanges for the karyoploteR plots
# somatic.mutations.all <- toGRanges(VEP_data_all_patients[, list(CHROM, POS, END_POS, ID, REF, ALT, region_type, mutation_type)])
# seqlevelsStyle(somatic.mutations.all) <- "UCSC"
# 
# # frequency plot for all combined chromosomes. The frequency is directly calculated by karyoteR
# #png(file="Results/Frequency plots/all_patients_frequency.png", width=465, height=225, units = "mm", res=300)
# 
# pp <- getDefaultPlotParams(plot.type = 4)
# pp$data1inmargin <- 0
# pp$bottommargin <- 20
# 
# kp <- plotKaryotype(plot.type=4, ideogram.plotter = NULL,
#                     labels.plotter = NULL, plot.params = pp)
# kpAddCytobandsAsLine(kp)
# kpAddChromosomeNames(kp, srt=45)
# #kpPlotDensity(kp, data = somatic.mutations.all, window.size = 10e6)
# kpfrequency <- kpPlotDensity(kp, data = somatic.mutations.all, col = "#F6D55C", window.size = 10e6)
# kpfrequency$latest.plot$computed.values$frequency
# kpfrequency$latest.plot$computed.values$windows
# kpAddLabels(kp, labels = c("Mutation frequency"), srt=90, pos=1, label.margin = 0.04)
# legend("topleft", 
#        title = "Window size",
#        legend = "1 Mb",
#        #legend = c("1 Mb", "10 Mb"),
#        xjust = 1,
#        lty = 1,
#        col = "#F6D55C",
#        #col = c("#F6D55C", "black"),
#        lwd = 2,
#        seg.len = 1,
#        xpd = TRUE,
#        bty = "n",
#        x.intersp = 0.5)
# 
# #dev.off()
# 
# # frequency plot for each chromosome
# for(chromosome in seqlevels(somatic.mutations.all)) {
#   # png(file=glue("Results/frequency plots/{chromosome}_frequency.png"),
#   #     width=465, height=225, units = "mm", res=300)
#   
#   kp <- plotKaryotype(plot.type=4, ideogram.plotter = NULL,
#                       labels.plotter = NULL, plot.params = pp, chromosomes = chromosome)
#   kpAddCytobandsAsLine(kp)
#   kpAddChromosomeNames(kp, srt=45)
#   kpPlotDensity(kp, data = somatic.mutations.all)
#   kpAddLabels(kp, labels = c("Mutation frequency"), srt=90, pos=1, label.margin = 0.04)
#   
#   # dev.off()
# }




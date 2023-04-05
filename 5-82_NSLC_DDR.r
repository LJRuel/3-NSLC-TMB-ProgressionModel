
## TODO
## - Compare TMBs with SBS signatures
## - Make a list of cancer protectors, cancer drivers and DNA damage repair genes to find mutations

library(GenomicFeatures)
library(tidyverse)
library(data.table)
library(glue)
library(xlsx)
library(biomaRt)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Clinicopathological data - keeping only adenocarcinomas
clin_data <- as.data.table(read.xlsx("Data/NCI_NeverSmoker_n92_20210812_TMB_drivers.xlsx", sheetIndex = 1))
clin_data <- clin_data[histology == 3, ]
Patients_list <- gsub("NSLC-", "",sort(clin_data[, Patient_ID]))

# Fetching and merging all patients variants
VEP_data_all_patients <- rbindlist(lapply(Patients_list, function(x) fread(glue("Data/VEP_82_NSLC_TMB/VEP_NSLC-{x}.csv"))), use.names = TRUE)
# Adding patient ID (for easier access than the excisting "ID" field)
VEP_data_all_patients[, Patient_ID := gsub(":.*", "", ID)]

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
                                                          time_RpFS < 1825 ~ "<5"))
# Changing comorbidities NA to None
surv.dt[, comorbidities:=replace_na(comorbidities, "None")]


# DNA damage response genes
DDR.genes <- setDT(readxl::read_excel("Data/NIHMS962902_DDR_genes.xlsx", range= "A4:AC280"))
DDR.gene.mutations <- setNames(lapply(DDR.genes[,`Gene Symbol`], function(gene) VEP_data_all_patients[SYMBOL == gene]), 
                               nm= DDR.genes[,`Gene Symbol`])


lapply(names(DDR.gene.mutations), function(x) mutate(surv.dt, ))

################################################################################
#################################### Archive ###################################
################################################################################

# Searching for mutated DDR genes by genome coordinates
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
DDR.genes.info  <- setDT(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                               filters=c('entrezgene_id', 'hgnc_symbol'),
                               values=list(DDR.genes[, `Entrez Gene ID`], DDR.genes[, `Gene Symbol`]),
                               mart=ensembl))
# Order chromosome names
DDR.genes.info <- DDR.genes.info[str_order(chromosome_name, numeric = TRUE),]
# Remove non standard chromosome names
DDR.genes.info <- DDR.genes.info[chromosome_name %in% c(1:22, 'X', 'Y')]

DDR.genes.mutated <-  setNames(lapply(1:nrow(DDR.genes.info), function(x) VEP_data_all_patients[POS >= DDR.genes.info[x, start_position] &
                                                                                                  END_POS <= DDR.genes.info[x, end_position] &
                                                                                                  CHROM == DDR.genes.info[x, chromosome_name]]),
                               nm = DDR.genes.info[, hgnc_symbol])



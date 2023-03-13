
library(GenomicFeatures)
library(tidyverse)
library(data.table)
library(glue)
library(xlsx)
library(karyoploteR)
library(regioneR)
library(ggbio)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Clinicopathological data - keeping only adenocarcinomas
clin_data <- as.data.table(read.xlsx("Data/NCI_NeverSmoker_n92_20210812_TMB_drivers.xlsx", sheetIndex = 1))
clin_data <- clin_data[histology == 3, ]
Patients_list <- gsub("NSLC-", "",sort(clin_data[, Patient_ID]))

# Fetching and merging all patients variants
VEP_data_all_patients <- rbindlist(lapply(Patients_list, function(x) fread(glue("Data/VEP_82_NSLC_TMB/VEP_NSLC-{x}.csv"))), use.names = TRUE)

# Window size used to compute sliding windows
window_size <- 5000
Mb_length <-  window_size/1e6

## Merging all chromosome sliding windows densities together
chrom_name <- c(1:22, 'X', 'Y')
# Reading density data of chromosomes with a specific window size
chrom_density <- rbindlist(lapply(chrom_name, 
                                  function(x) fread(glue("Data/Sliding windows and density/chromosome_{x}_sw_density_{Mb_length}Mb.csv"))[, chrom:=glue("{x}")]),
                           use.names = TRUE)

# 
attr(chrom_density$chrom, 'glue') <- NULL

# Get the windows with the most variants
setorder(chrom_density, -variant_density)
hotspot.cutoff <- ceiling(nrow(chrom_density)*0.001) # Keep top 0.1% of windows with the most variants
hotspot.windows <- chrom_density[1:hotspot.cutoff]
names(hotspot.windows) <- c("POS", "END_POS", "FREQUENCY", "CHROM")
hotspot.windows[, CHROM:=as.character(CHROM)] # Editing 

# Merge the neighboring windows if regions are wanted (combined neighboring windows)
#hotspot.windows.reduced <- GRanges(hotspot.windows) %>% GenomicRanges::reduce() %>% as.data.table() %>% .[, strand:=NULL]
#names(hotspot.windows.reduced) <- c("CHROM", "POS", "END_POS", "WIDTH")

# Extract the variants inside each window
hotspots.list <- lapply(1:nrow(hotspot.windows), 
                        function(x){rbindlist(list(hotspot.windows[x], 
                                                   setorder(VEP_data_all_patients[POS >= hotspot.windows[x, POS] & 
                                                                                    END_POS <= hotspot.windows[x, END_POS] & 
                                                                                    CHROM == hotspot.windows[x, CHROM]], POS)),
                                              fill=TRUE)
                        }
)

# Order the windows based on the number of variants they include
hotspots.list.nrow <- sapply(1:length(hotspots.list), function(x) nrow(hotspots.list[[x]]))
hotspots.list.ordered <- hotspots.list[order(hotspots.list.nrow, decreasing = TRUE)]

################################################################################
################################ HOTSPOT PLOTS #################################
################################################################################

# Reference: https://stackoverflow.com/questions/58696329/is-it-possible-to-over-ride-the-x-axis-range-in-r-package-ggbio-when-using-autop
library(ggbio)
library(EnsDb.Hsapiens.v75)
ensdb <- EnsDb.Hsapiens.v75

hotspots.list.ordered[1]
mut.7.142465001 <- hotspots.list.ordered[1] %>% makeGRangesFromDataFrame(., ignore.strand = TRUE,
                                                               start.field="POS",
                                                               end.field="END_POS",
                                                               strand.field="CHROM")
gene <- autoplot(ensdb, which = mut.7.142465001)
gene +
  #geom_point(data = as.data.table(mut.7.142465001)[2:length(mut.7.142465001)], aes(x=start, y=1)) +
  geom_histogram(data = as.data.table(mut.7.142465001)[2:length(mut.7.142465001)], aes(x=start), binwidth = 1)

ggplot() +
  geom_genemodel()


library(GenomicFeatures)
library(AnnotationHub)
library(tidyverse)
library(data.table)
library(glue)
library(xlsx)
library(progress)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Chromosome ranges for sliding window
chrom_grngs <- as(Seqinfo(genome="GRCh37"), "GRanges")
chrom_grngs <- keepSeqlevels(chrom_grngs, c(1:22, 'X', 'Y'), pruning.mode = 'coarse')

# Clinicopathological data - keeping only adenocarinomas
clin_data <- as.data.table(read.xlsx("Data/NCI_NeverSmoker_n92_20210812_TMB_drivers.xlsx", sheetIndex = 1))
clin_data <- clin_data[histology == 3, ]
Patients_list <- gsub("NSLC-", "",sort(clin_data[, Patient_ID]))

# Fetching variant data for all patients
data_list <- lapply(Patients_list, function(x) fread(glue("Data/VEP_82_NSLC_TMB/VEP_NSLC-{x}.csv")))
NSLC_variants <- rbindlist(data_list, use.names = TRUE)
NSLC_variants <- NSLC_variants[order(NSLC_variants[,CHROM])]

# Define window size and step
window_size <- as.integer(5000)
step <- as.integer(1000)

# Initialize vectors to store results
for(chrom in as.character(seqnames(chrom_grngs))) {
  start_positions <- seq(1, width(chrom_grngs[chrom]) + window_size - 1, by = step)
  end_positions <- start_positions + window_size - 1
  assign(glue("positions.{chrom}"), data.table(window_start=start_positions, window_end=end_positions))
}

# Compute variant density for each sliding window
chrom_list <- objects()[grep("positions.", objects())] %>% str_sort(., numeric = TRUE)
for (chrom in chrom_list) {
  chrom_id <- gsub("positions.", "", chrom)
  variant_density <- numeric(nrow(get(chrom)))
  message(glue("Chromosome {chrom_id} progression"))
  pb <- progress_bar$new(format = "[:bar] :percent",
                         total = nrow(get(chrom)),
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  for(i in 1:nrow(get(chrom))) {
    pb$tick()
    variant_density[i] <- length(NSLC_variants[which(POS>=get(chrom)[i, window_start] & END_POS<=get(chrom)[i, window_end] & CHROM==chrom_id), ID])
  }
  get(chrom)[,variant_density:=variant_density]
}

# Exporting variant density for each chromosome. All patients are included for each chromosome
Mb_length <-  window_size/1e6
lapply(as.character(seqnames(chrom_grngs)), 
       function(x) fwrite(get(glue("positions.{x}")), 
                          file=glue("Data/Sliding windows and density/chromosome_{x}_sw_density_{Mb_length}Mb.csv"),
                          scipen = 10))



################################################################################
################################## ARCHIVE #####################################
################################################################################

#https://www.google.com/search?q=r+genomic+sliding+window&oq=r+genomic+sliding+window&aqs=edge..69i57j0i30i546i625.11434j0j1&sourceid=chrome&ie=UTF-8
# GRanges(VEP_data)
# gr <- makeGRangesFromDataFrame(VEP_data, keep.extra.columns = TRUE, ignore.strand = TRUE, 
#                                start.field = "POS", end.field = "END_POS")
# 
# start_positions <- setNames(sapply(seqnames(chrom_grngs), 
#                                    function(x){
#                                      seq(1, width(chrom_grngs[x]) + window_size - 1, by = step)
#                                    }), 
#                             nm=seqnames(chrom_grngs))
# 
# end_positions <- sapply(names(start_positions),
#                         function(x){
#                           start_positions[[x]] + window_size - 1
#                         })
# 
# for (i in seq_along(start_positions[[chrom]])) {
#   window_start <- start_positions[[chrom]][i]
#   window_end <- end_positions[[chrom]][i]
#   window_variants <- VEP_data[which(POS>=window_start & END_POS<=window_end & CHROM==chrom), ID]
#   variant_density[[chrom]][i] <- length(window_variants)
#   progress(i, max(seq_along(start_positions[[chrom]])))
  
  
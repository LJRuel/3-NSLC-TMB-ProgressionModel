
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

# Merge the neighboring windows if regions are wanted (combined neighboring windows)
hotspot.windows.reduced <- GRanges(hotspot.windows) %>% GenomicRanges::reduce() %>% as.data.table() %>% .[, strand:=NULL]
names(hotspot.windows.reduced) <- c("CHROM", "POS", "END_POS", "WIDTH")


### Unmerged neighboring windows
# If windows are not be merged, change the names of the columns to avoid names matching problems later on
names(hotspot.windows) <- c("POS", "END_POS", "FREQUENCY", "CHROM")
hotspot.windows[, CHROM:=as.character(CHROM)] # Editing 

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

### Merged neighboring windows
# Extract the variants inside each merged window
hotspots.list.merge <- lapply(1:nrow(hotspot.windows.reduced), 
                              function(x){rbindlist(list(hotspot.windows.reduced[x], 
                                                         setorder(VEP_data_all_patients[POS >= hotspot.windows.reduced[x, POS] & 
                                                                                          END_POS <= hotspot.windows.reduced[x, END_POS] & 
                                                                                          CHROM == hotspot.windows.reduced[x, CHROM]], POS)),
                                                    fill=TRUE)
                              }
)

# Order the merged windows based on the number of variants they include
hotspots.list.merged.nrow <- sapply(1:length(hotspots.list.merge), function(x) nrow(hotspots.list.merge[[x]]))
hotspots.list.merged.ordered <- hotspots.list.merge[order(hotspots.list.merged.nrow, decreasing = TRUE)]

################################################################################
################################ HOTSPOT PLOTS #################################
################################################################################

# Reference: https://stackoverflow.com/questions/58696329/is-it-possible-to-over-ride-the-x-axis-range-in-r-package-ggbio-when-using-autop
library(ggbio)
library(ggplotify)
library(gridExtra)
library(EnsDb.Hsapiens.v75)
ensdb <- EnsDb.Hsapiens.v75

# Some functions for plot x axis ticks
nearest_lower_limit <- function(value, divider) {
  value - (value %% divider)
}
nearest_upper_limit <- function(value, divider) {
  value + (value %% divider)
}


### Unmerged neighboring windows
mut.7.142465001 <- hotspots.list.ordered[1] %>% makeGRangesFromDataFrame(., ignore.strand = TRUE,
                                                                         start.field="POS",
                                                                         end.field="END_POS",
                                                                         strand.field="CHROM")
mut.7.142465001.dt <- as.data.table(mut.7.142465001)

freq <- ggplot(data = mut.7.142465001.dt[c(2:.N)], aes(x=start)) +
  #geom_point(data = as.data.table(mut.7.142465001)[2:length(mut.7.142465001)], aes(x=start, y=1)) +
  geom_histogram(binwidth = 1) +
  theme_classic()

endb.anno <- autoplot(ensdb, which = mut.7.142465001) + theme_classic()

cowplot::plot_grid(freq, endb.anno@ggplot + xlim(mut.7.142465001.dt[1, c(2,3)]), ncol=1, align = "v")

### Merged neighboring windows
window.1 <- hotspots.list.merged.ordered[1] %>% makeGRangesFromDataFrame(., ignore.strand = TRUE,
                                                                         start.field="POS",
                                                                         end.field="END_POS",
                                                                         strand.field="CHROM",
                                                                         keep.extra.columns = TRUE)
window.1.dt <- as.data.table(window.1)
freq <- ggplot(data = window.1.dt[c(2:.N)], aes(x=start, fill = region_type)) +
  geom_bar(width = 25, position = "identity") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("#ff595e", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93"),
                    name = "Variant type",
                    labels = c("exonic", "intergenic", "intronic", "splice site")) +
  ylab("Variant Frequency") +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(window.1.dt[1, start]-1, window.1.dt[1, end]+1),
                  # here, max(table()) gets the computed frequency value by geom_bar need for upper y_lim
                  ylim = c(0, max(table(window.1.dt[, start]))), 
                  clip = "off")

endb.anno <- autoplot(ensdb, 
                      which = window.1, 
                      columns = c("tx_biotype", "symbol"),
                      names.expr = "symbol:tx_biotype")@ggplot + 
  theme_classic() +
  xlim(window.1.dt[1, start]-1, window.1.dt[1, end]) +
  scale_x_continuous(n.breaks = length(unique(window.1.dt[, start])), 
                     breaks = seq(nearest_lower_limit(window.1.dt[1, start]-1, 5000), 
                                  nearest_upper_limit(window.1.dt[1, end], 5000), 
                                  by = 5000)) +
  xlab(NULL) +
  ylab("Ensembl Transcripts")

cowplot::plot_grid(freq, endb.anno, ncol=1, align = "v", axis = "b")

## TODO
## - Align the axes from top and bottom plots when the legend is on
## - Loop to plot the other hotspots

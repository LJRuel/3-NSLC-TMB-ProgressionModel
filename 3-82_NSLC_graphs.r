
## TODO
# - Make a graph of all variants with a Manhattan-like style plot (one for each chromosome?).
#   One color per region type.
# - Make a similar graph that combines all the patients.
# - Compare TMBs with SBS signatures

library(GenomicFeatures)
library(tidyverse)
library(data.table)
library(glue)
library(xlsx)
library(karyoploteR)
library(regioneR)
library(RColorBrewer)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Clinicopathological data - keeping only adenocarcinomas
clin_data <- as.data.table(read.xlsx("Data/NCI_NeverSmoker_n92_20210812_TMB_drivers.xlsx", sheetIndex = 1))
clin_data <- clin_data[histology == 3, ]
Patients_list <- gsub("NSLC-", "",sort(clin_data[, Patient_ID]))

# Fetching data for each patient.
for(Patient in Patients_list) {
  assign(glue("VEP_data_{Patient}"), fread(glue("Data/VEP_82_NSLC_TMB/VEP_NSLC-{Patient}.csv")))
}

# Merging all patients variants
VEP_data_list <- objects()[grep("VEP_data_0", objects())]
VEP_data_all_patients <- rbindlist(lapply(1:length(VEP_data_list), function(x) get(VEP_data_list[x])), use.names = TRUE)
VEP_data_all_patients <- VEP_data_all_patients[order(VEP_data_all_patients[,CHROM])]

# Preparing the GRanges for the karyoploteR plots
somatic.mutations.all <- toGRanges(VEP_data_all_patients[, list(CHROM, POS, END_POS, ID, REF, ALT, region_type, mutation_type)])
seqlevelsStyle(somatic.mutations.all) <- "UCSC"

################################################################################
############################### Density plots ##################################
################################################################################

############# With the density pre calculated by sliding windows ###############

# Window size used to compute sliding windows
window_size <- 1e6

## Merging all chromosome sliding windows densities together
chrom_name <- c(1:22, 'X', 'Y')
# Reading data for each chrosome
for(chromosome in chrom_name) {
  assign(glue("chrom_{chromosome}_density_{window_size}Mb"), fread(glue("Data/Sliding windows and density/chromosome_{chromosome}_sw_density_{window_size}Mb.csv"))[, chrom:=glue("{chromosome}")])
}
# Extracting the data read in the previous loop
expr <- regexpr("chrom_.*_density_.*Mb", objects()) # get the objects indices related to the read data
chrom_density.list <- str_sort(regmatches(objects(), expr), numeric = TRUE) # extract the objects names
chrom_density.all <- rbindlist(lapply(chrom_density.list, function(x) get(x))) # merging the chromosomes densities data

chrom_density.all$chrom <- factor(chrom_density.all$chrom, levels=chrom_name)

# 3rd quartile cutoff for hotspot regions
hotspot_cutoff <- quantile(chrom_density.all$variant_density)[4]

# Plotting the density
ggplot(data = chrom_density.all, aes(x = window_start, y = variant_density, color = chrom)) + 
  geom_line() + 
  geom_hline(yintercept = hotspot_cutoff, color = "red") + 
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.placement='outside') +
  scale_color_hue() +
  xlab("Chromosomes") + 
  ylab("Variant Density") +
  facet_grid(~chrom,
             space = "free_x", 
             switch = "x")


############## With the density directly calculated by karyoteR ################

# Density plot for all combined chromosomes. The density is directly calculated by karyoteR
#png(file="Results/Density plots/all_patients_density.png", width=465, height=225, units = "mm", res=300)

pp <- getDefaultPlotParams(plot.type = 4)
pp$data1inmargin <- 0
pp$bottommargin <- 20

kp <- plotKaryotype(plot.type=4, ideogram.plotter = NULL,
                    labels.plotter = NULL, plot.params = pp)
kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, srt=45)
#kpPlotDensity(kp, data = somatic.mutations.all, window.size = 10e6)
kpDensity <- kpPlotDensity(kp, data = somatic.mutations.all, col = "#F6D55C", window.size = 10e6)
kpDensity$latest.plot$computed.values$density
kpDensity$latest.plot$computed.values$windows
kpAddLabels(kp, labels = c("Mutation density"), srt=90, pos=1, label.margin = 0.04)
legend("topleft", 
       title = "Window size",
       legend = "1 Mb",
       #legend = c("1 Mb", "10 Mb"),
       xjust = 1,
       lty = 1,
       col = "#F6D55C",
       #col = c("#F6D55C", "black"),
       lwd = 2,
       seg.len = 1,
       xpd = TRUE,
       bty = "n",
       x.intersp = 0.5)

#dev.off()

# Density plot for each chromosome
for(chromosome in seqlevels(somatic.mutations.all)) {
  # png(file=glue("Results/Density plots/{chromosome}_density.png"),
  #     width=465, height=225, units = "mm", res=300)
  
  kp <- plotKaryotype(plot.type=4, ideogram.plotter = NULL,
                      labels.plotter = NULL, plot.params = pp, chromosomes = chromosome)
  kpAddCytobandsAsLine(kp)
  kpAddChromosomeNames(kp, srt=45)
  kpPlotDensity(kp, data = somatic.mutations.all)
  kpAddLabels(kp, labels = c("Mutation density"), srt=90, pos=1, label.margin = 0.04)
  
  # dev.off()
}

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
# kpAddLabels(kp, labels = c("Mutation density"), srt=90, pos=1, label.margin = 0.04, r0=0.71, r1=1)


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
#   kpAddLabels(kp, labels = c("Mutation density"), srt=90, pos=1, label.margin = 0.04, r0=0.71, r1=1)
#   
#   dev.off()
# }

################################################################################
################################### ARCHIVES ###################################
################################################################################

# ## Using ChromoMap to show variant density on chromosomes.
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
# assign(glue("chrom_{chromosome}_density"), fread(glue("Data/Sliding windows and density/chromosome_{chromosome}_sw_density.csv")))
# assign(glue("chrom_{chromosome}_density.gr"), makeGRangesFromDataFrame(get(glue("chrom_{chromosome}_density"))[, chrom:=glue("chr{chromosome}")], 
#                                                                        keep.extra.columns = TRUE, 
#                                                                        ignore.strand = TRUE, 
#                                                                        seqnames.field = "chrom", 
#                                                                        start.field = "window_start",
#                                                                        end.field = "window_end"))
# }
# 





## TODO
# - Make a graph of all variants with a Manhattan-like style plot (one for each chromosome?).
#   One color per region type.
# - Make a similar graph that combines all the patients.
# - Compare TMBs with SBS signatures

BiocManager::install("karyoploteR")

library(GenomicFeatures)
library(tidyverse)
library(data.table)
library(glue)
library(xlsx)
library(karyoploteR)
library(regioneR)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Clinicopathological data - keeping only adenocarinomas
clin_data <- as.data.table(read.xlsx("Data/NCI_NeverSmoker_n92_20210812_TMB_drivers.xlsx", sheetIndex = 1))
clin_data <- clin_data[histology == 3, ]
Patients_list <- gsub("NSLC-", "",sort(clin_data[, Patient_ID]))

for(Patient in Patients_list) {
  assign(glue("VEP_data_{Patient}"), fread(glue("Data/VEP_82_NSLC_TMB/VEP_NSLC-{Patient}.csv")))
}

VEP_data_list <- objects()[grep("VEP_data_0", objects())]

VEP_data_all_patients <- rbindlist(lapply(1:length(VEP_data_list), function(x) get(VEP_data_list[x])), use.names = TRUE)

VEP_data_all_patients <- VEP_data_all_patients[order(VEP_data_all_patients[,CHROM])]

VEP_data_all_patients <- VEP_data_all_patients %>%  mutate(END_POS = case_when(mutation_type == "SNV" ~ POS,
                                                                               mutation_type == "insertion" ~ POS,
                                                                               mutation_type == "deletion" ~ POS+nchar(REF))) %>% relocate(END_POS, .after = POS)

# karyoploteR mutations graphs. Reference: https://bernatgel.github.io/karyoploter_tutorial/

somatic.mutations <- toGRanges(VEP_data_all_patients[, list(CHROM, POS, END_POS, ID, REF, ALT, region_type, mutation_type)])
seqlevelsStyle(somatic.mutations) <- "UCSC"
bla <- VEP_data_all_patients[grep("0001", VEP_data_all_patients[, ID], fixed = TRUE)]
somatic.mutations.0001 <- toGRanges(VEP_data_all_patients[VEP_data_all_patients[, ID] %in% VEP_data_0001[, ID],
                                                          list(CHROM, POS, END_POS, ID, REF, ALT, region_type, mutation_type)])

kp <- plotKaryotype(plot.type=4)
kpPlotRainfall(kp, data = somatic.mutations.0001)
variant.colors <- getVariantsColors(somatic.mutations$REF, somatic.mutations$ALT)






?str_split

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





## TODO
# - Make a graph of all variants with a Manhattan-like style plot (one for each chromosome?).
#   One color per region type.
# - Make a similar graph that combines all the patients.
# - Compare TMBs with SBS signatures

#BiocManager::install("MutationalPatterns")

library(VariantAnnotation)
library(GenomicFeatures)
library(AnnotationHub)
library(tidyverse)
library(data.table)
library(glue)
library(chromoMap) # reference: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04556-z

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

Patient="0001"

VEP_data <- fread(glue("Data/VEP_82_NSLC_TMB/VEP_NSLC-{Patient}.csv"))

VEP_data <- VEP_data %>%  mutate(END_POS = case_when(mutation_type == "SNV" ~ POS,
                                                     mutation_type == "insertion" ~ POS,
                                                     mutation_type == "deletion" ~ POS+nchar(REF))) %>% relocate(END_POS, .after = POS)

ah = AnnotationHub()
#AnnotationHub::query(ah, c('gtf', 'Homo_sapiens', 'GRCh37'))
GRCh37.gtf<- ah[['AH10684']]


## subset the gtf files for only protein_coding genes and lincRNAs
GRCh37.gtf<- GRCh37.gtf[GRCh37.gtf$gene_biotype %in% c('protein_coding', 'lincRNA')]

## make a txdb and keep conventional chromosomes
GRCh37.txdb <- makeTxDbFromGRanges(GRCh37.gtf)
GRCh37.txdb <- keepSeqlevels(GRCh37.txdb, c(1:22, 'X', 'Y'), pruning.mode = 'coarse')
chrom_grngs <- as(seqinfo(GRCh37.txdb), 'GRanges')
chrom_data <- data.table()
chrom_data[, CHROM:=seqnames(chrom_grngs)%>%as.character()]
chrom_data[, POS:=start(chrom_grngs)%>%as.numeric()]
chrom_data[, END_POS:=end(chrom_grngs)%>%as.numeric()]

anno_data <- VEP_data[, list(ID, CHROM, POS, END_POS)]

chromo <- chromoMap(list(chrom_data), list(anno_data))

VEP_data[, loci:= chromo$x$chData[[1]]$loci]
window_count <- as.data.table(table(VEP_data[, loci]))
VEP_data <- merge(VEP_data, window_count, by.x = "loci", by.y = "V1")
setnames(VEP_data, "N", "loci_var_count")
VEP_data <- VEP_data %>% mutate(loci_var_cat = case_when(loci_var_count < 5 ~ "[1, 5]",
                                                         loci_var_count >= 5 & loci_var_count <10 ~ "[5, 10[",
                                                         loci_var_count >=10 ~ "[10, ["))

anno_data_2 <- VEP_data[, list(ID, CHROM, POS, END_POS, loci_var_count, loci_var_cat)]

chromoMap(list(chrom_data), list(anno_data_2),
          data_based_color_map = T,
          data_type = "numeric",
          plots = "scatter",
          plot_filter = list(c("col","byCategory")),
          ch2D.colors = c("pink3","orange3","purple"))

bla <- window_count[order(window_count[, N], decreasing = TRUE),]
median(bla$N)

################################################################################
################################### ARCHIVES ###################################
################################################################################

window_length <- chromo$x$chData[[1]]$loci_end[1] - chromo$x$chData[[1]]$loci_start[1]

for(chromosome in unique(VEP_data[, CHROM])) {
  chr_indices <- c(first(grep(glue("chromap-{chromosome}-"),  chromo$x$chData[[1]]$loci)), last(grep(glue("chromap-{chromosome}-"),  chromo$x$chData[[1]]$loci)))
  assign(glue("chr_{chromosome}"), chr_indices)
  assign(glue("chr_{chromosome}_wd_length"), chromo$x$chData[[1]]$loci_end[chr_indices[1]] - chromo$x$chData[[1]]$loci_start[chr_indices[1]])
}





## TODO
# - Make a graph of all variants with a Manhattan-like style plot (one for each chromosome?).
#   One color per region type.
# - Make a similar graph that combines all the patients.
# - Compare TMBs with SBS signatures

library(VariantAnnotation)
library(tidyverse)
library(data.table)
library(glue)
library(chromoMap)

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
table(GRCh37.gtf$gene_biotype)

## make a txdb and keep conventional chromosomes
GRCh37.txdb <- makeTxDbFromGRanges(GRCh37.gtf)
GRCh37.txdb <- keepSeqlevels(GRCh37.txdb, c(1:22, 'X', 'Y'), pruning.mode = 'coarse')
chrom_grngs <- as(seqinfo(GRCh37.txdb), 'GRanges')
chrom_data <- data.table()
chrom_data[, CHROM:=seqnames(chrom_grngs)%>%as.character()]
chrom_data[, POS:=start(chrom_grngs)%>%as.numeric()]
chrom_data[, END_POS:=end(chrom_grngs)%>%as.numeric()]

anno_data <- VEP_data[, list(ID, CHROM, POS, END_POS)]

## TODO
# - Déterminer les windows
# - Compter le nb de variants par window
# - Donner ce nombre pour le mapping

chromoMap(list(chrom_data), list(anno_data),
          data_based_color_map = T,
          data_type = "numeric",
          plots = "bar",
          aggregate_func = "count",
          heat_map = F)


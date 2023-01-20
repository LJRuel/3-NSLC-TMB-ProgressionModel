if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("VariantAnnotation")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(BSgenome.Hsapiens.UCSC.hg19)
#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#BSgenome <- BSgenome.Hsapiens.UCSC.hg19

library(readr)
library(GenomicRanges)
library(VariantAnnotation)
library(stringr)
library(tidyr)
library(dplyr)
library(glue)


# Reference : https://bioconductor.org/packages/devel/bioc/vignettes/VariantAnnotation/inst/doc/VariantAnnotation.pdf

Patient <- "0001"

# Linux
if(Sys.info()['sysname'] == "Linux") {
  # Clinicopathological data - keeping only adenocarinomas
  clin_data <- read.xlsx("/mnt/sda2/TMB/Data/VEP_92_NSLC//NCI_NeverSmoker_n92_20210812_TMB_drivers.xlsx", sheetIndex = 1)
  adeno_list <- clin_data$Patient_ID[clin_data$histology == "3"]
  
  
  VEP.first_line <- grep('##', readLines(glue("/mnt/sda2/TMB/Data/VEP_92_NSLC/NSLC_{Patient}.vcf")), invert = TRUE, fixed = TRUE, )[1]
  VCF.info_lines <- grep('##INFO=', readLines(glue("/mnt/sda2/TMB/Data/VEP_92_NSLC/NSLC_{Patient}.vcf")), fixed = TRUE)
  VEP.format_line <- last(VCF.info_lines)
  VCF.info <- readLines(glue("/mnt/sda2/TMB/Data/VEP_92_NSLC/NSLC_{Patient}.vcf"))[c(VCF.info_lines[1:length(VCF.info_lines)-1])] %>% gsub(".*ID=", "", .) %>% gsub(",.*", "", .)
  VEP.format <- readLines(glue("/mnt/sda2/TMB/Data/VEP_92_NSLC/NSLC_{Patient}.vcf"))[VEP.format_line] %>% gsub(".*Format: ", "", .) %>% gsub("\">", "", .) %>% str_split(., "\\|") %>% unlist()
  VEP_data <- read.delim(glue("/mnt/sda2/TMB/Data/VEP_92_NSLC/NSLC_{Patient}.vcf"), sep = "\t", skip = VEP.first_line-1, header = TRUE)
  VEP_data <- VEP_data[VEP_data$X.CHROM != "MT",]
}
  
# Windows
if(Sys.info()['sysname'] == "Windows") {
  # Clinicopathological data - keeping only adenocarinomas
  clin_data <- read.xlsx("Data/NCI_NeverSmoker_n92_20210812_TMB_drivers.xlsx", sheetIndex = 1)
  adeno_list <- clin_data$Patient_ID[clin_data$histology == "3"]
  
  
  VEP.header <- grep('##', readLines(glue("Data/VEP_92_NSLC/NSLC_{Patient}.vcf")), invert = TRUE, fixed = TRUE, )[1]
  VCF.info_lines <- grep('##INFO=', readLines(glue("Data/VEP_92_NSLC/NSLC_{Patient}.vcf")), fixed = TRUE)
  VEP.format_line <- last(VCF.info_lines)
  VCF.info <- readLines(glue("Data/VEP_92_NSLC/NSLC_{Patient}.vcf"))[c(VCF.info_lines[1:length(VCF.info_lines)-1])] %>% gsub(".*ID=", "", .) %>% gsub(",.*", "", .)
  VEP.format <- readLines(glue("Data/VEP_92_NSLC/NSLC_{Patient}.vcf"))[VEP.format_line] %>% gsub(".*Format: ", "", .) %>% gsub("\">", "", .) %>% str_split(., "\\|") %>% unlist()
  VEP_data <- read.delim(glue("Data/VEP_92_NSLC/NSLC_{Patient}.vcf", sep = "\t", skip = VEP.header-1, header = TRUE))

}

# Data formatting
VEP.extra <- strsplit(VEP_data$INFO, ";\\s*(?=[^;]+$)", perl=TRUE)
VEP_data$INFO <- sapply(VEP.extra, "[[", 1)
VEP_data$EXTRA <- sapply(VEP.extra, "[[", 2)
VEP_data$EXTRA <- lapply(VEP_data$EXTRA, gsub, pattern = 'CSQ=', replacement='')
VEP_data.EXTRA_split_rows <- separate_rows(VEP_data, EXTRA, sep=",")
VEP_data.EXTRA_split_cols <- separate(VEP_data.EXTRA_split_rows, EXTRA, c(VEP.format), sep = "\\|")
assign(glue("VEP_data.NSLC_{Patient}"), VEP_data.EXTRA_split_cols)




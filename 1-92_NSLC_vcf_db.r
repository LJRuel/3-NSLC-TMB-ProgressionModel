if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("GenomicRanges")
#BiocManager::install("VariantAnnotation")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(BSgenome.Hsapiens.UCSC.hg19)
#txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#BSgenome <- BSgenome.Hsapiens.UCSC.hg19

library(tidyverse)
library(glue)
library(xlsx)
library(GenomicRanges)
#library(VariantAnnotation)


Patient <- "0001"

# Mutation types annotated by VEP are available here: https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html
region_types <- list()
region_types[["intergenic"]] <- c("downstream_gene_variant", "upstream_gene_variant", "intergenic_variant")
region_types[["intron"]] <- c("intron_variant", "non_coding_transcript_variant")
region_types[["regulatory"]] <- c("regulatory_region_variant", "TF_binding_site_variant")
region_types[["exon"]] <- c("non_coding_transcript_exon_variant", "missense_variant", "synonymous_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "coding_sequence_variant", "frameshift_variant", "stop_gained")
region_types[["splice"]] <- c("splice_region_variant", "splice_donor_variant", "splice_acceptor_variant", "splice_polypyrimidine_tract_variant", "splice_donor_region_variant", "splice_donor_5th_base_variant")

region_types

## Linux
if(Sys.info()['sysname'] == "Linux") {
  # Clinicopathological data - keeping only adenocarinomas
  clin_data <- read.xlsx("/home/l_ruelou01/NSLC_TMB-Model/Data/NCI_NeverSmoker_n92_20210812_TMB_drivers.xlsx", sheetIndex = 1)
  adeno_list <- clin_data$Patient_ID[clin_data$histology == "3"]
  
  VEP.raw <- readLines(glue("/mnt/sda2/TMB/Data/VEP_92_NSLC/NSLC_{Patient}.vcf"))
  VEP.header <- grep('##', VEP.raw, invert = TRUE, fixed = TRUE, )[1]
  VCF.info_lines <- grep('##INFO=', VEP.raw, fixed = TRUE)
  VEP.format_line <- last(VCF.info_lines)
  VCF.info <- VEP.raw[c(VCF.info_lines[1:length(VCF.info_lines)-1])] %>% gsub(".*ID=", "", .) %>% gsub(",.*", "", .)
  VEP.format <- VEP.raw[VEP.format_line] %>% gsub(".*Format: ", "", .) %>% gsub("\">", "", .) %>% str_split(., "\\|") %>% unlist()
  VEP_data <- read.delim(glue("/mnt/sda2/TMB/Data/VEP_92_NSLC/NSLC_{Patient}.vcf"), sep = "\t", skip = VEP.header-1, header = TRUE)
  # Removing MT chromosome variants
  VEP_data <- VEP_data[VEP_data$X.CHROM != "MT",]
  # Adding patient ID to the ID field
  VEP_data$ID <- paste0(glue("{Patient}:"), VEP_data$ID)
}
  
## Windows
if(Sys.info()['sysname'] == "Windows") {
  # Clinicopathological data - keeping only adenocarinomas
  clin_data <- read.xlsx("Data/NCI_NeverSmoker_n92_20210812_TMB_drivers.xlsx", sheetIndex = 1)
  adeno_list <- clin_data$Patient_ID[clin_data$histology == "3"]
  
  # Getting the description and data lines from vcf files
  VEP.raw <- readLines(glue("Data/VEP_92_NSLC/NSLC_{Patient}.vcf"))
  VEP.header <- grep('##', VEP.raw, invert = TRUE, fixed = TRUE, )[1]
  VCF.info_lines <- grep('##INFO=', VEP.raw, fixed = TRUE)
  VEP.format_line <- last(VCF.info_lines)
  VCF.info <- VEP.raw[c(VCF.info_lines[1:length(VCF.info_lines)-1])] %>% gsub(".*ID=", "", .) %>% gsub(",.*", "", .)
  VEP.format <- VEP.raw[VEP.format_line] %>% gsub(".*Format: ", "", .) %>% gsub("\">", "", .) %>% str_split(., "\\|") %>% unlist()
  VEP_data <- read.delim(glue("Data/VEP_92_NSLC/NSLC_{Patient}.vcf"), sep = "\t", skip = VEP.header-1, header = TRUE)
  
  # Removing MT chromosome variants
  VEP_data <- VEP_data[VEP_data$X.CHROM != "MT",]
  
  # Adding patient ID to the ID field
  VEP_data$ID <- paste0(glue("{Patient}:"), VEP_data$ID)
}

## Data formatting
VEP.extra <- strsplit(VEP_data$INFO, ";\\s*(?=[^;]+$)", perl=TRUE)
VEP_data$INFO <- sapply(VEP.extra, "[[", 1)
VEP_data$EXTRA <- sapply(VEP.extra, "[[", 2)
VEP_data$EXTRA <- lapply(VEP_data$EXTRA, gsub, pattern = 'CSQ=', replacement='')
VEP_data.EXTRA_split_rows <- separate_rows(VEP_data, EXTRA, sep=",")
VEP_data.EXTRA_split_cols <- separate(VEP_data.EXTRA_split_rows, EXTRA, c(VEP.format), sep = "\\|")

# The "PICK" field represents the suggested entry to keep for each variant according to VEP. See info here: http://useast.ensembl.org/info/docs/tools/vep/script/vep_other.html#most_severe_eg
VEP_data.pick <- VEP_data.EXTRA_split_cols[VEP_data.EXTRA_split_cols$PICK==1,]

# Data harmonizing: reporting the worst (first) consequence of variants when multiple consequences are annotated by VEP.
# Order of consequence severity: https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html
worst_consequence <- sapply(strsplit(VEP_data.pick$Consequence, "&", perl=TRUE), "[[", 1)

# Verifying that all consequences are registered in the available region types (region_types variable). 
# If a consequence type is not registered in the region_types variable, then a message prompts to add 
#   the missing consequence to the proper region type, according to https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html.
if(!all(unique(worst_consequence) %in% as.vector(unlist(region_types)))) {
  message("A missing region type has been detected. Please add the missing
          region in the propoer region_types variable (see script) according to 
          https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html.")
} else {
  VEP_data.pick$region_type <- unlist(sapply(1:nrow(VEP_data.pick), function(x) {names(region_types[grep(worst_consequence[x], region_types)])}))
  VEP_data.pick <- VEP_data.pick %>% relocate(region_type, .after=Consequence)
}

# Final data set for each patient
assign(glue("VEP_data.NSLC_{Patient}"), VEP_data.pick)


## TO DO : 
# - Comparer VEP avec l'annotation des GRanges (surtout pour les types de régions enregistrées).
# - Est-ce qu'on garde les synonymous variants?



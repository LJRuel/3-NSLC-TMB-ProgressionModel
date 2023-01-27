

## TODO
# - Repeat the TMB operations, but for each mutation type (SNVs, insertion, deletion, CNV(?), ...).
# - Find the respective total size of each region type.
# - Get CNV data -> need BAM files.
# - Make a graph of all variants with a Manhattan-like style plot (one for each chromosome?).
#   One color per region type.
# - Make a similar graph that combines all the patients.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("AnnotationHub")

library(GenomicRanges)
library(GenomicFeatures)
library(tidyverse)
library(data.table)
library(xlsx)
library(glue)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################################################################
######################### Data cleaning and formatting ######################### 
################################################################################

# Mutation types annotated by VEP are available here: https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html
region_types <- list()
region_types[["intergenic"]] <- c("downstream_gene_variant", 
                                  "upstream_gene_variant", 
                                  "intergenic_variant")

region_types[["intron"]] <- c("intron_variant", 
                              "non_coding_transcript_variant")

region_types[["regulatory"]] <- c("regulatory_region_variant", 
                                  "TF_binding_site_variant")

region_types[["exon"]] <- c("non_coding_transcript_exon_variant", 
                            "missense_variant", 
                            "synonymous_variant", 
                            "3_prime_UTR_variant", 
                            "5_prime_UTR_variant", 
                            "coding_sequence_variant", 
                            "frameshift_variant", 
                            "stop_gained")

region_types[["splice"]] <- c("splice_region_variant", 
                              "splice_donor_variant", 
                              "splice_acceptor_variant", 
                              "splice_polypyrimidine_tract_variant", 
                              "splice_donor_region_variant", 
                              "splice_donor_5th_base_variant")

# Initating the hash table containing the TMB data for patients.
VEP.TMB_records <- new.env(hash = TRUE)

# Clinicopathological data - keeping only adenocarinomas
clin_data <- as.data.table(read.xlsx("Data/NCI_NeverSmoker_n92_20210812_TMB_drivers.xlsx", sheetIndex = 1))
clin_data <- clin_data[histology == 3, ]

# Patients list
Patient <- "0001"

## Importing vcf data
if(glue("NSLC-{Patient}") %in% clin_data[, Patient_ID]) {
  VEP.raw <- readLines(glue("Data/VEP_92_NSLC/NSLC_{Patient}.vcf"))
  VEP.header <- grep('##', VEP.raw, invert = TRUE, fixed = TRUE, )[1]
  VCF.info_lines <- grep('##INFO=', VEP.raw, fixed = TRUE)
  VEP.format_line <- last(VCF.info_lines)
  VCF.info <- VEP.raw[c(VCF.info_lines[1:length(VCF.info_lines)-1])] %>% gsub(".*ID=", "", .) %>% gsub(",.*", "", .)
  VEP.format <- VEP.raw[VEP.format_line] %>% gsub(".*Format: ", "", .) %>% gsub("\">", "", .) %>% str_split(., "\\|") %>% unlist()
  VEP_data <- setDT(read.delim(glue("Data/VEP_92_NSLC/NSLC_{Patient}.vcf"), sep = "\t", skip = VEP.header-1, header = TRUE))
  setnames(VEP_data, "X.CHROM", "CHROM")
  
  # Removing MT chromosome variants
  VEP_data <- VEP_data[CHROM != "MT",]
  
  # Adding patient ID to the ID field
  VEP_data[, ID := paste0(glue("{Patient}:"), VEP_data[,ID])]
} else {
  message(glue("The patient {Patient} is excluded. Handling next patient."))
  break
}

## Data formatting
VEP_data[, EXTRA := strsplit(VEP_data[,INFO], ";\\s*(?=[^;]+$)", perl=TRUE)]
VEP_data[, INFO := sapply(VEP_data[,EXTRA], "[[", 1)]
VEP_data[, EXTRA := sapply(VEP_data[,EXTRA], "[[", 2)]
VEP_data[, EXTRA := lapply(VEP_data[,EXTRA], gsub, pattern = 'CSQ=', replacement='')]
VEP_data <- setDT(separate_rows(VEP_data, EXTRA, sep=","))
VEP_data <- setDT(separate(VEP_data, EXTRA, c(VEP.format), sep = "\\|"))

# The "PICK" field represents the suggested entry to keep for each variant according to VEP. See info here: http://useast.ensembl.org/info/docs/tools/vep/script/vep_other.html#most_severe_eg
VEP_data <- VEP_data[PICK==1,]

# Data harmonizing: reporting the worst (first) consequence of variants when multiple consequences are annotated by VEP.
# Order of consequence severity: https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html
worst_consequence <- sapply(strsplit(VEP_data[,Consequence], "&", perl=TRUE), "[[", 1)

# Verifying that all consequences are registered in the available region types (region_types variable). 
# If a consequence type is not registered in the region_types variable, then a message prompts to add 
#   the missing consequence to the proper region type, according to https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html.
if(!all(unique(worst_consequence) %in% as.vector(unlist(region_types)))) {
  missing_consequence <- unique(worst_consequence[!(worst_consequence %in% as.vector(unlist(region_types)))])
  stop("Missing region type detected: ", paste(missing_consequence, collapse = ", "), 
  ".\nPlease add the missing region in the proper 'region_types' variable (see script)",  
  "\naccording to https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html.")
  cat()
  
} else {
  VEP_data[, region_type := unlist(sapply(1:nrow(VEP_data), function(x) {names(region_types[grep(worst_consequence[x], region_types)])}))]
  # Bringing region_type column next to Consequence column
  setcolorder(VEP_data, colnames(VEP_data)[c(1:13, ncol(VEP_data), 14:(ncol(VEP_data)-1))])
}


# Adding mutation type
VEP_data <- VEP_data %>%
  mutate(mutation_type = case_when(nchar(REF) == nchar(ALT) ~ "SNV",
                                   nchar(REF) > nchar(ALT) ~ "deletion",
                                   nchar(REF) < nchar(ALT) ~ "insertion"))
# Bringing mutation_type column next to region_type column
setcolorder(VEP_data, colnames(VEP_data)[c(1:14, ncol(VEP_data), 15:(ncol(VEP_data)-1))])


# Final data set for each patient
#assign(glue("VEP_data.NSLC_{Patient}"), VEP_data)
#fwrite(VEP_data.NSLC_0001, glue("Data/VEP_92_NSLC/VEP_NSLC_{Patient}.csv"))

################################################################################
################################# Regions size #################################
################################################################################

# References: https://gist.github.com/crazyhottommy/4681a30700b2c0c1ee02cbc875e7c4e9
## make a txdb
GRCh37.txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
GRCh37.txdb <- keepSeqlevels(GRCh37.txdb, paste0("chr", c(1:22, "X", "Y")), pruning.mode = "coarse")

## exons
exons <- exonsBy(GRCh37.txdb, "gene")
exons2<- exonicParts(GRCh37.txdb, linked.to.single.gene.only = TRUE)

width(exons) %>% unlist() %>% sum()
width(exons2) %>% unlist() %>% sum()

## introns
introns <- intronicParts(GRCh37.txdb, linked.to.single.gene.only = TRUE)
width(introns) %>% sum()


################################################################################
########################## TMB according to regions ############################
################################################################################

# The GRCh37 whole genome size is 3,101,788,170 bases = 3,102 Mb. Reference: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/
# The WGS TMB contained in the clinical data (clin_data) has been calculated using a general size of 3000 Mb.
# WGS TMB without synonymous variants
VEP.WGS_TMB.no_syn <- round(nrow(VEP_data[Consequence != "synonymous_variant"]) / 3102, digits = 3)
# WGS TMB with synonymous variants
VEP.WGS_TMB.syn <- round(nrow(VEP_data) / 3102, digits = 3)

# TMBs depending on the region. They are calculated based on the whole genome size instead of their own region total cumulative size.
# This only allows to get the ratio of each region based TMB on the total WGS TMB.
# Same with WGS TMB, two version of th exon TMB include synonymous variants or not.
VEP.intergenic_TMB <- round(nrow(VEP_data[region_type == "intergenic"]) / 3102, digits = 3)
VEP.intron_TMB <- round(nrow(VEP_data[region_type == "intron"]) / 3102, digits = 3)
VEP.regulatory_TMB <- round(nrow(VEP_data[region_type == "regulatory"]) / 3102, digits = 3)
VEP.exon_TMB.syn <- round(nrow(VEP_data[region_type == "exon"]) / 3102, digits = 3)
VEP.exon_TMB.no_syn <- round(nrow(VEP_data[region_type == "exon" & Consequence != "synonymous_variant"]) / 3102, digits = 3)
VEP.splice_TMB <- round(nrow(VEP_data[region_type == "splice"]) / 3102, digits = 3)

VEP.TMBs_list.no_syn <- c(VEP.intergenic_TMB, VEP.intron_TMB, VEP.regulatory_TMB, VEP.exon_TMB.no_syn, VEP.splice_TMB)
VEP.TMBs_list.syn <- c(VEP.intergenic_TMB, VEP.intron_TMB, VEP.regulatory_TMB, VEP.exon_TMB.syn, VEP.splice_TMB)

# Making a table for WGS and region based TMB calculated above. Not including synonymous variants.
if(round(sum(VEP.TMBs_list.no_syn), digits = 3) != VEP.WGS_TMB.no_syn) {
  message("The WGS TMB with no synonymous_variant is incomplete. Check that the variant count from VEP_data is correct.")
} else {
  VEP.TMBs.no_syn <- setNames(
    list(
      setNames(c(VEP.WGS_TMB.no_syn, VEP.TMBs_list.no_syn), c("total_WGS_TMB", paste0(names(region_types), "_WGS_TMB"))),
      setNames(sapply(1:length(VEP.TMBs_list.no_syn), function(x) round(VEP.TMBs_list.no_syn[x]/VEP.WGS_TMB.no_syn, digits = 4)), paste0(names(region_types), "_WGS_TMB_ratio"))), 
    c("WGS_TMBs_no_synonymous", "WGS_TMBs_no_synonymous_ratios")
  )
}

# Making a table for WGS and region based TMB calculated above. Including synonymous variants.
if(round(sum(VEP.TMBs_list.syn), digits = 3) != VEP.WGS_TMB.syn) {
  message("The WGS TMB with no synonymous_variant is incomplete. Check that the variant count from VEP_data is correct.")
} else {
  VEP.TMBs.syn <- setNames(
    list(
      setNames(c(VEP.WGS_TMB.syn, VEP.TMBs_list.syn), c("total_WGS_TMB", paste0(names(region_types), "_WGS_TMB"))),
      setNames(sapply(1:length(VEP.TMBs_list.syn), function(x) round(VEP.TMBs_list.syn[x]/VEP.WGS_TMB.syn, digits = 4)), paste0(names(region_types), "_WGS_TMB_ratio"))), 
    c("WGS_TMBs_synonymous", "WGS_TMBs_synonymous_ratios")
  )
}

VEP.TMB_records[[glue("NSLC_{Patient}")]] = c(VEP.TMBs.no_syn, VEP.TMBs.syn)


################################################################################
####################### TMB according to mutation types ########################
################################################################################









################################################################################
################################### ARCHIVES ###################################
################################################################################

## Comparing GRanges annotation with VEP annotations
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("GenomicRanges")
#BiocManager::install("VariantAnnotation")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
BSgenome <- BSgenome.Hsapiens.UCSC.hg19

vcf <- readVcf(glue("/mnt/sda2/TMB/Data/92_patients_NSLC_filtered_VCFS/NSLC_{Patient}.vcf"), "hg19")
vcf <- keepSeqlevels(vcf, c(1:22, "X", "Y"), pruning.mode = "coarse")
seqlevels(vcf) <- paste0("chr", c(1:22, "X", "Y"))
vcf.range <- rowRanges(vcf)

all.variants <- locateVariants(vcf.range, txdb, AllVariants())
mcols(all.variants)[c("LOCATION", "PRECEDEID", "FOLLOWID")]
which(all.variants == trim(all.variants))
which(end(all.variants) > seqlengths(BSgenome)[as.character(seqnames(all.variants))])

# Did any coding variants match more than one gene?
splt.match <- split(mcols(all.variants)$GENEID, mcols(all.variants)$QUERYID)
table(sapply(splt.match, function(x) length(unique(x)) > 1))

# Summarize the number of variants by gene ID
splt.genes <- split(mcols(all.variants)$QUERYID, mcols(all.variants)$GENEID)
sapply(splt.genes, function(x) length(unique(x)))


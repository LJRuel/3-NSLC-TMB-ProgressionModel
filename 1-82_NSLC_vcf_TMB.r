
## TODO
# - Get CNV data -> need BAM files.

## Questions
# - Est-ce qu'on garde les synonymnous variants?
# - Est-ce qu'on garde les pseudogenes dans les variants annotes?

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("AnnotationHub")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
#BiocManager::install("biomaRt")
#BiocManager::install("ensembldb")
#BiocManager::install("EnsDb.Hsapiens.v75")

library(GenomicFeatures)
library(AnnotationHub)
library(tidyverse)
library(data.table)
library(xlsx)
library(glue)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################################################################
######################### Data cleaning and formatting ######################### 
################################################################################

# Mutation types annotated by VEP are available here: https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html
region_types <- list()
region_types[["exons"]] <- c("non_coding_transcript_exon_variant", 
                            "missense_variant", 
                            "synonymous_variant", 
                            "3_prime_UTR_variant", 
                            "5_prime_UTR_variant", 
                            "coding_sequence_variant", 
                            "frameshift_variant", 
                            "stop_gained",
                            "stop_lost",
                            "inframe_deletion",
                            "inframe_insertion",
                            "start_lost",
                            "stop_retained_variant",
                            "protein_altering_variant")

region_types[["splice"]] <- c("splice_region_variant", 
                              "splice_donor_variant", 
                              "splice_acceptor_variant", 
                              "splice_polypyrimidine_tract_variant", 
                              "splice_donor_region_variant", 
                              "splice_donor_5th_base_variant")

region_types[["introns"]] <- c("intron_variant", 
                              "non_coding_transcript_variant")

region_types[["regulatory"]] <- c("regulatory_region_variant", 
                                  "TF_binding_site_variant",
                                  "TFBS_ablation")

region_types[["intergenic"]] <- c("downstream_gene_variant", 
                                  "upstream_gene_variant", 
                                  "intergenic_variant",
                                  "mature_miRNA_variant")

# Initating the hash table containing the TMB data for patients.
VEP.TMB_records <- data.table()

# Clinicopathological data - keeping only adenocarinomas
clin_data <- as.data.table(read.xlsx("Data/NCI_NeverSmoker_n92_20210812_TMB_drivers.xlsx", sheetIndex = 1))
clin_data <- clin_data[histology == 3, ]

################################################################################
################################# Regions size #################################
################################################################################

# References: https://gist.github.com/crazyhottommy/4681a30700b2c0c1ee02cbc875e7c4e9
# The following code does not need to be run again. The relevent information is
# contained in "GRCh37.region_sizes" variable at the end of this section. 
# To run this section again, select the commented lines and remove the 
# comment marks with Ctrl+Shift+C

ah = AnnotationHub()
#AnnotationHub::query(ah, c('gtf', 'Homo_sapiens', 'GRCh37'))
GRCh37.gtf<- ah[['AH10684']]


## subset the gtf files for only protein_coding genes and lincRNAs
GRCh37.gtf<- GRCh37.gtf[GRCh37.gtf$gene_biotype %in% c('protein_coding', 'lincRNA')]
table(GRCh37.gtf$gene_biotype)

## make a txdb and keep conventional chromosomes
GRCh37.txdb <- makeTxDbFromGRanges(GRCh37.gtf)
GRCh37.txdb <- keepSeqlevels(GRCh37.txdb, c(1:22, 'X', 'Y'), pruning.mode = 'coarse')


## cds (not compiled in the total genome size because included in exons)
cds <- cdsBy(GRCh37.txdb, 'gene') %>% unlist() %>% unstrand() %>% GenomicRanges::reduce()
sum(width(cds))


## exons: 101,578,353 bp
exons <- exonicParts(GRCh37.txdb) %>% unstrand() %>% GenomicRanges::reduce() 
exons_size <- as.double(sum(width(exons)))
mcols(exons)$region <- "exons"


## introns: 1,345,484,609 bp
# The overlap with some exons is also removed by GenomicRanges::setdiff() function.
introns <- intronicParts(GRCh37.txdb) %>% unstrand() %>% GenomicRanges::setdiff(., exons)
introns_size <- as.double(sum(width(introns)))
mcols(introns)$region <- "introns"


## regulatory regions: 221,273,310 bp
# Gathering regulatory regions data from Ensembl
edb <- EnsDb.Hsapiens.v75
ensembl_reg = useEnsembl(biomart="regulation", dataset = "hsapiens_regulatory_feature", GRCh=37)
listAttributes(ensembl_reg)
all.regulatory <- as.data.table(getBM(attributes = c("chromosome_name", "bound_seq_region_start",
                                                     "bound_seq_region_end", "feature_type_name"),
                                      mart = ensembl_reg))

# Rearranging the regulatory regions output
setcolorder(all.regulatory, c("chromosome_name", "bound_seq_region_start", "bound_seq_region_end", "feature_type_name"))
all.regulatory <- all.regulatory[order(chromosome_name, bound_seq_region_start)]
regulatory <- GRanges(all.regulatory) %>% unstrand() %>% GenomicRanges::reduce()
regulatory <- keepSeqlevels(regulatory, c(1:22, "X", "Y"), pruning.mode = "coarse")
regulatory_rest <- c(exons, introns) %>% sort() %>% GenomicRanges::setdiff(regulatory, .)
regulatory_size <- as.double(sum(width(regulatory_rest)))
mcols(regulatory_rest)$region <- "regulatory"

## intergenic regions: 1,427,341,140 bp
chrom_grngs <- as(seqinfo(GRCh37.txdb), 'GRanges')
collapsed_tx <- GenomicRanges::reduce(transcripts(GRCh37.txdb))
strand(collapsed_tx) <- '*'
intergenic <- GenomicRanges::setdiff(chrom_grngs, collapsed_tx)
# Intergenic region size and removing overlapping regions with other annotated regions.
non_inter_ranges <- c(exons, introns, regulatory_rest) %>% sortSeqlevels() %>% sort()
intergenic_rest <- GenomicRanges::setdiff(intergenic, non_inter_ranges)
mcols(intergenic_rest)$region <- "intergenic"
intergenic_size <- as.double(sum(width(intergenic_rest)))


# *** Not needed. 5' and 3' UTR.
five_UTR <- fiveUTRsByTranscript(GRCh37.txdb) %>% unlist() %>% GenomicRanges::reduce()
sum(width(five_UTR))

three_UTR <- threeUTRsByTranscript(GRCh37.txdb) %>% unlist() %>% GenomicRanges::reduce()
sum(width(three_UTR))

## total genome size: 3095677412 bp
genome <- c(exons, introns, regulatory_rest, intergenic_rest) %>% sort()
genome_size <- sum(c(exons_size, introns_size, regulatory_size, intergenic_size))


## Final region sizes used for TMB calculation. 
GRCh37.region_sizes <- as.data.table(setNames(list(3095677412, 101578353, 1345484609, 221273310, 1427341140), 
                                              list("genome_size", "exons_size",
                                                    "introns_size", "regulatory_size",
                                                   "intergenic_size")))

## Values are in bp. TMB requires values to be in megabase, therefore a simple division by 10^6 is needed when calculating TMBs.
GRCh37.region_sizes <- round(GRCh37.region_sizes/1e6, digits = 2)


################################################################################
################################## Main loop ###################################
################################################################################

# Patients list
Patients_list <- gsub("NSLC-", "",sort(clin_data[, Patient_ID]))

for(Patient in Patients_list) {
  
  Patient = "0001"
  
  message(glue("Processing TMB values for patient {Patient}."))
  
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
    "\naccording to https://grch37.ensembl.org/info/genome/variation/prediction/predicted_data.html")
    
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
  
  # Adding end position of variant. Relevant for deletion.
  VEP_data <- VEP_data %>%
    mutate(END_POS = case_when(mutation_type == "SNV" ~ POS,
                               mutation_type == "insertion" ~ POS,
                               mutation_type == "deletion" ~ POS+nchar(REF))) %>% relocate(END_POS, .after = POS)
  
  # Bringing mutation_type column next to region_type column
  setcolorder(VEP_data, colnames(VEP_data)[c(1:14, ncol(VEP_data), 15:(ncol(VEP_data)-1))])
  
  
  # Final data set for each patient
  #fwrite(VEP_data, glue("Data/VEP_82_NSLC_TMB/VEP_NSLC-{Patient}.csv"))
 
  
  
  ##############################################################################
  ######################### TMB according to regions ###########################
  ##############################################################################
  
  # The GRCh37 whole genome size is 3,129,589,526 bases = 3,192 Mb as calculated above (see GRCh37.region_sizes variable).
  # The WGS TMB contained in the clinical data (clin_data) has been calculated using a general size of 3000 Mb.
  # WGS TMB without synonymous variants
  VEP.WGS_TMB.no_syn <- round(nrow(VEP_data[Consequence != "synonymous_variant"]) / GRCh37.region_sizes$genome_size, digits = 3)
  # WGS TMB with synonymous variants
  VEP.WGS_TMB.syn <- round(nrow(VEP_data) / GRCh37.region_sizes$genome_size, digits = 3)
  
  # TMBs depending on the region. They are calculated based on the whole genome size instead of their own region total cumulative size.
  # This only allows to get the ratio of each region based TMB on the total WGS TMB.
  # Same with WGS TMB, two version of th exon TMB include synonymous variants or not.
  
  ################
  ### Option 1 ###
  ################
  
  ## TMB calculated according to each region's size, e.g. intron TMB with intron variants and intron size.
  # Using this option limits the region types to exon, intron and intergenic. 
  # Regulatory are included in intergenic and splice are included in exon.
  VEP.exon_TMB.syn <- round(nrow(VEP_data[region_type == "exons" | region_type == "splice"]) / GRCh37.region_sizes$exons_size, digits = 3)
  VEP.exon_TMB.no_syn <- round(nrow(VEP_data[region_type == "exons" & Consequence != "synonymous_variant" | region_type == "splice"]) / GRCh37.region_sizes$exons_size, digits = 3)
  VEP.intron_TMB <- round(nrow(VEP_data[region_type == "introns"]) / GRCh37.region_sizes$introns_size, digits = 3)
  VEP.intergenic_TMB <- round(nrow(VEP_data[region_type == "intergenic"]) / GRCh37.region_sizes$intergenic_size, digits = 3)
  VEP.regulatory_TMB <- round(nrow(VEP_data[region_type == "regulatory"]) / GRCh37.region_sizes$regulatory_size, digits = 3)
  
  VEP.TMBs_list.no_syn <- list(VEP.exon_TMB.no_syn, VEP.intron_TMB, VEP.regulatory_TMB, VEP.intergenic_TMB)
  VEP.TMBs_list.syn <- list(VEP.exon_TMB.syn, VEP.intron_TMB, VEP.regulatory_TMB, VEP.intergenic_TMB)
  
  # Making a table for WGS and region based TMB calculated above. Not including synonymous variants.
  VEP.TMBs.no_syn <- as.data.table(setNames(c(VEP.WGS_TMB.no_syn, VEP.TMBs_list.no_syn),
                                               c("genome_TMB", 
                                                 paste0(names(region_types)[c(1,3:5)], "_TMB")))) # Splice variants are included in exons
  
  # Making a table for WGS and region based TMB calculated above. Including synonymous variants.
  VEP.TMBs.syn <- as.data.table(setNames(c(VEP.WGS_TMB.syn, VEP.TMBs_list.syn), 
                                           c("genome_TMB_with_synonymous", 
                                             paste0(names(region_types)[c(1,3:5)], "_TMB_with_synonymous")))) # Splice variants are included in exons
  
  ################
  ### Option 2 ###
  ################
  
  ### TMB calculated according to the whole genome size, e.g. intron TMB with intron variants and genome size.
  # This option is representing the different TMBs as a proportion of the whole genome TMB, i.e. their weights on the WGS TMB.
  VEP.WGS.exon_TMB.syn <- round(nrow(VEP_data[region_type == "exons"]) / GRCh37.region_sizes$genome_size, digits = 3)
  VEP.WGS.exon_TMB.no_syn <- round(nrow(VEP_data[region_type == "exons" & Consequence != "synonymous_variant"]) / GRCh37.region_sizes$genome_size, digits = 3)
  VEP.WGS.splice_TMB <- round(nrow(VEP_data[region_type == "splice"]) / GRCh37.region_sizes$genome_size, digits = 3)
  VEP.WGS.intergenic_TMB <- round(nrow(VEP_data[region_type == "intergenic"]) / GRCh37.region_sizes$genome_size, digits = 3)
  VEP.WGS.regulatory_TMB <- round(nrow(VEP_data[region_type == "regulatory"]) / GRCh37.region_sizes$genome_size, digits = 3)
  VEP.WGS.intron_TMB <- round(nrow(VEP_data[region_type == "introns"]) / GRCh37.region_sizes$genome_size, digits = 3)
  
  VEP.WGS.TMBs_list.no_syn <- list(VEP.WGS.exon_TMB.no_syn, VEP.WGS.splice_TMB, VEP.WGS.intron_TMB, VEP.WGS.regulatory_TMB, VEP.WGS.intergenic_TMB)
  VEP.WGS.TMBs_list.syn <- list(VEP.WGS.exon_TMB.syn, VEP.WGS.splice_TMB, VEP.WGS.intron_TMB, VEP.WGS.regulatory_TMB, VEP.WGS.intergenic_TMB)
  
  # Making a table for WGS and region based TMB calculated above. Not including synonymous variants.
  VEP.WGS.TMBs.no_syn <- as.data.table(setNames(c(VEP.WGS_TMB.no_syn, VEP.WGS.TMBs_list.no_syn),
                                                  c("genome_TMB", paste0(names(region_types), "_WGS_TMB"))))
  
  # Making a table for WGS and region based TMB calculated above. Including synonymous variants.
  VEP.WGS.TMBs.syn <- as.data.table(setNames(c(VEP.WGS_TMB.syn, VEP.WGS.TMBs_list.syn), 
                                               c("genome_TMB", paste0(names(region_types), "_WGS_TMB_with_synonymous"))))
  
  # Making a table for the region specific TMB to WGS TMB ratio. Not including synonymous variants.
  VEP.WGS.TMBs_ratio.no_syn <- setNames(round(sapply(1:length(VEP.WGS.TMBs_list.no_syn),
                                                     function(x) VEP.WGS.TMBs_list.no_syn[[x]]/VEP.WGS_TMB.no_syn), digits = 3),
                                        c(paste0(names(region_types), "_TMB_ratio")))
  
  # Making a table for the region specific TMB to WGS TMB ratio. Including synonymous variants.
  VEP.WGS.TMBs_ratio.syn <- setNames(c(round(sapply(1:length(VEP.WGS.TMBs_list.syn),
                                                     function(x) VEP.WGS.TMBs_list.syn[[x]]/VEP.WGS_TMB.syn), digits = 3)),
                                        c(paste0(names(region_types), "_TMB_ratio_syn")))
  
  
  all.TMB.regions <- list(unlist(c(VEP.TMBs.no_syn,
                                   VEP.TMBs.syn[, .(genome_TMB_with_synonymous,
                                                    exons_TMB_with_synonymous)])),
                          unlist(c(VEP.WGS.TMBs.no_syn,
                                   VEP.WGS.TMBs.syn[, .(exons_WGS_TMB_with_synonymous)])),
                          VEP.WGS.TMBs_ratio.no_syn,
                          VEP.WGS.TMBs_ratio.syn)
  
  all.TMB.regions <- setNames(all.TMB.regions, c("TMB_per_region", 
                                                 "WGS_TMB_per_region", 
                                                 "WGS_TMB_per_region_ratio", 
                                                 "WGS_TMB_per_region_ratio_with_synonymous"))
  
  ## TODO
  ## Change order in lists
  ## Put a description of every list type found in all.TMB.regions.
  setcolorder(all.TMB.regions, colnames(all.TMB.regions)[c(1, 6, 2, 7, 3:5, 8, ncol(all.TMB.regions), 9:(ncol(all.TMB.regions)-1))])
  
  ##############################################################################
  ###################### TMB according to mutation types #######################
  ##############################################################################
  
  VEP.SNV_TMB.no_syn <- round(nrow(VEP_data[mutation_type == "SNV" & Consequence != "synonymous_variant"]) / GRCh37.region_sizes$genome_size, digits = 3)
  VEP.SNV_TMB.syn <- round(nrow(VEP_data[mutation_type == "SNV"]) / GRCh37.region_sizes$genome_size, digits = 3)
  
  VEP.del_TMB.no_syn <- round(nrow(VEP_data[mutation_type == "deletion" & Consequence != "synonymous_variant"]) / GRCh37.region_sizes$genome_size, digits = 3)
  VEP.del_TMB.syn <- round(nrow(VEP_data[mutation_type == "deletion"]) / GRCh37.region_sizes$genome_size, digits = 3)
  
  VEP.ins_TMB.no_syn <- round(nrow(VEP_data[mutation_type == "insertion" & Consequence != "synonymous_variant"]) / GRCh37.region_sizes$genome_size, digits = 3)
  VEP.ins_TMB.syn <- round(nrow(VEP_data[mutation_type == "insertion"]) / GRCh37.region_sizes$genome_size, digits = 3)
  
  for(mutation in unique(VEP_data[, mutation_type])) {
    exons <- round(nrow(VEP_data[(region_type == "exons" | region_type == "splice") & mutation_type == mutation & Consequence != "synonymous_variant"]) / GRCh37.region_sizes$exons_size, digits = 3)
    assign(glue("exons_{mutation}_TMB"), exons)
    exons_syn <- round(nrow(VEP_data[(region_type == "exons" | region_type == "splice") & mutation_type == mutation]) / GRCh37.region_sizes$exons_size, digits = 3)
    assign(glue("exons_syn_{mutation}_TMB"), exons_syn)
    introns <- round(nrow(VEP_data[region_type == "introns" & mutation_type == mutation]) / GRCh37.region_sizes$introns_size, digits = 3)
    assign(glue("introns_{mutation}_TMB"), introns)
    intergenic <- round(nrow(VEP_data[(region_type == "intergenic" | region_type == "regulatory") & mutation_type == mutation]) / GRCh37.region_sizes$intergenic_size, digits = 3)
    assign(glue("intergenic_{mutation}_TMB"), intergenic)
  }
  
  all.TMB.mutations <- data.table(exons_SNV_TMB, exons_deletion_TMB, exons_insertion_TMB,
                                  exons_syn_SNV_TMB, exons_syn_deletion_TMB, exons_syn_insertion_TMB,
                                  introns_SNV_TMB, introns_deletion_TMB, introns_insertion_TMB,
                                  intergenic_SNV_TMB, intergenic_deletion_TMB, intergenic_insertion_TMB)
  
  
  ##############################################################################
  ########################## TMB data consolidation ############################
  ##############################################################################
  
  # Consolidation of TMB by region type and TMB by mutation type.
  
  # To access data: VEP.TMB_records[["NSLC_PatientID"]]$TMB_by_xx$yy
  #   where xx is either "TMB_by_region" or "TMB_by_mutation" and yy is an inner 
  #   xx list variable (see all.TMB.regions and all.TMB.mutations variables).
  #   E.g. VEP.TMB_records[["NSLC_0001"]]$TMB_by_region$exons_TMB
  
  VEP.TMB_records <- rbindlist(list(VEP.TMB_records, c(setNames(list(glue("NSLC-{Patient}")), "Patient_ID"), all.TMB.regions, all.TMB.mutations)))
  fwrite(VEP.TMB_records, glue("Data/VEP_82_NSLC_TMB/VEP_NSLC-{Patient}_TMB.csv"))
} # End of main loop

################################################################################
################################### ARCHIVES ###################################
################################################################################

# ## Comparing GRanges annotation with VEP annotations
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# #BiocManager::install("GenomicRanges")
# #BiocManager::install("VariantAnnotation")
# #BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# #BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
# 
# library(VariantAnnotation)
# library(BSgenome.Hsapiens.UCSC.hg19)
# BSgenome <- BSgenome.Hsapiens.UCSC.hg19
# 
# vcf <- readVcf(glue("/mnt/sda2/TMB/Data/92_patients_NSLC_filtered_VCFS/NSLC_{Patient}.vcf"), "hg19")
# vcf <- keepSeqlevels(vcf, c(1:22, "X", "Y"), pruning.mode = "coarse")
# seqlevels(vcf) <- paste0("chr", c(1:22, "X", "Y"))
# vcf.range <- rowRanges(vcf)
# 
# all.variants <- locateVariants(vcf.range, txdb, AllVariants())
# mcols(all.variants)[c("LOCATION", "PRECEDEID", "FOLLOWID")]
# which(all.variants == trim(all.variants))
# which(end(all.variants) > seqlengths(BSgenome)[as.character(seqnames(all.variants))])
# 
# # Did any coding variants match more than one gene?
# splt.match <- split(mcols(all.variants)$GENEID, mcols(all.variants)$QUERYID)
# table(sapply(splt.match, function(x) length(unique(x)) > 1))
# 
# # Summarize the number of variants by gene ID
# splt.genes <- split(mcols(all.variants)$QUERYID, mcols(all.variants)$GENEID)
# sapply(splt.genes, function(x) length(unique(x)))
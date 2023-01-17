if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = '3.16')
BiocManager::install("GenomicRanges")
BiocManager::install("VariantAnnotation")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

library(GenomicRanges)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(stringr)
library(tidyr)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
BSgenome <- BSgenome.Hsapiens.UCSC.hg19

# Reference : https://bioconductor.org/packages/devel/bioc/vignettes/VariantAnnotation/inst/doc/VariantAnnotation.pdf

# Linux
VEP.first_line <- grep('##', readLines("/mnt/sda2/TMB/Data/VEP_92_NSLC/NSLC_0001.vep"), invert = TRUE, fixed = TRUE, )[1]
VEP_data <- read.delim("/mnt/sda2/TMB/Data/VEP_92_NSLC/NSLC_0001.vep", sep = "\t", skip = VEP.first_line-1, header = TRUE)
# Windows
#VEP.first_line <- grep('##', readLines("../VEP_92_NSLC/NSLC_0001.vcf"), invert = TRUE, fixed = TRUE, )[1]
#VEP_data <- read.delim("../VEP_92_NSLC/NSLC_0001.vcf", sep = "\t", skip = VEP.first_line-1, header = TRUE)

max_extra_col_names <- which.max(str_count(VEP_data$Extra, ";"))
extra_col_names <- gsub("=.*", "", str_split_1(VEP_data[max_extra_col_names,ncol(VEP_data)], ";"))

# Now the col names have been set, split all the data rows for each new col name.
for (i in nrow(VEP_data)) {
  
  
}

B = sapply(str_split(VEP_data$Extra, ";"), gsub, pattern =".*=", replacement = "")
B <- str_split_fixed(VEP_data$Extra, ";", n=2)

# Match the result with the column name.
separate(VEP_data, Extra, c(extra_col_names), sep = ";")
gsub(".*=", "", str_split_1(VEP_data$Extra[1], ";"))


VCF <- keepSeqlevels(VCF, c(1:22, "X", "Y"), pruning.mode = "coarse")
seqlevels(VCF) <- paste0("chr", c(1:22, "X", "Y"))
VCF.range <- rowRanges(VCF)

all.variants <- locateVariants(VCF.range, txdb, AllVariants())

which(all.variants == trim(all.variants))
which(end(all.variants) > seqlengths(BSgenome)[as.character(seqnames(all.variants))])
all.variants[all.variants$QUERYID == 153]


# Did any coding variants match more than one gene?
splt <- split(mcols(all.variants)$GENEID, mcols(all.variants)$QUERYID)
table(sapply(splt, function(x) length(unique(x)) > 1))

# Summarize the number of variants by gene ID
splt <- split(mcols(all.variants)$QUERYID, mcols(all.variants)$GENEID)
sapply(splt, function(x) length(unique(x)))

?`trim,GenomicRanges-method`

# ? faire : 
# Comparer annotations avec celles de VEP.
# 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicRanges")
BiocManager::install("VariantAnnotation")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")


library(GenomicRanges)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
BSgenome <- BSgenome.Hsapiens.UCSC.hg19

# Reference : https://bioconductor.org/packages/devel/bioc/vignettes/VariantAnnotation/inst/doc/VariantAnnotation.pdf

VCF <- readVcf("../92_patients_NSLC_filtered_VCFS/NSLC_0001.vcf", "hg19")
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

# À faire : 
# Comparer annotations avec celles de VEP.
# 

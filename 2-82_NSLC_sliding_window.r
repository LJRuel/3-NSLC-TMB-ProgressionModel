
library(GenomicFeatures)
library(tidyverse)
library(data.table)
library(glue)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

Patient="0001"

VEP_data <- fread(glue("Data/VEP_82_NSLC_TMB/VEP_NSLC-{Patient}.csv"))

#https://www.google.com/search?q=r+genomic+sliding+window&oq=r+genomic+sliding+window&aqs=edge..69i57j0i30i546i625.11434j0j1&sourceid=chrome&ie=UTF-8


library(SummarizedExperiment)
#library(tidyverse)
library(dplyr)
library(diffUTR)
library(edgeR)
library(DESeq2)

gtfPATH = "./data/gencode.v38.annotation.gtf"

# since we lost rowdata while merging 2500+ bam files, we add it back.
bins <- prepareBins(gtfPATH, stranded=FALSE)
saveRDS(bins, file="./outputs/diffUTR_bins_gencode_v38.RDS")

hist()

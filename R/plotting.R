library(SummarizedExperiment)
library(tidyverse)
library(dplyr)
library(diffUTR)
library(edgeR)


gtfPATH = "./data/gencode.v38.annotation.gtf"
bins <- prepareBins(gtfPATH, stranded=FALSE)



hist(log2(width(ranges(bins))), main='log2BinWidth')

sum(width(ranges(bins)))

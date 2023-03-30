library(SummarizedExperiment)
#library(tidyverse)
library(dplyr)
library(diffUTR)
library(edgeR)
library(DESeq2)

gtfPATH = "./data/gencode.v38.annotation.gtf"

exp <- readRDS("outputs/exonExpression_hdl.RDS") # explore.R generates this file

# since we lost rowdata while merging 2500+ bam files, we add it back.
bins <- prepareBins(gtfPATH, stranded=FALSE)
exp_se <- SummarizedExperiment(assays=list(counts=exp),
                              rowData=bins)
# Look into bin width and filter
MinBinWidth <- 45 # it was chosen so that the number of the remaining bins < 1M
row_width_filter <- width(ranges(bins)) > MinBinWidth
exp_se <- exp_se[row_width_filter,]

# remove outlying bins. (bins with too many counts)
maxAvgCount = 1000000
myRowMeans <- assays(exp_se)$counts %>% rowMeans(.)
myRowMeans_sorted <- sort(myRowMeans, decreasing=T)
myRowMeans_sorted[1:20]
row_excessive_filter <- myRowMeans < maxAvgCount
exp_se <- exp_se[row_excessive_filter,]


# Filter for lowly expressed exons.
min_num_samples = floor(ncol(exp_se) * 0.015)
expression_level_filter <- rowSums(cpm(assays(exp_se)$counts > 3)) >= min_num_samples
exp_passing <- exp_se[expression_level_filter, ] # no need sample filter as we are using samples used for TWAS

# counts assay is currently a data.frame, need to convert into data.matrix
assays(exp_passing)$counts <- data.matrix(assays(exp_passing)$counts)

# VST transformation
# convert RangedSummarizedExperiment object to DESeqDataSet
dds <- DESeqDataSet(exp_passing, design=~1)
dds <- estimateSizeFactors(dds)
saveRDS(dds, file="./outputs/exonExpression_hdl_preVST.RDS")

#vst_obj <- vst(dds, blind=T)

# output the
#saveRDS(vst_obj, file="./outputs/exonExpression_hdl_preCovariatesAdj.RDS")



# 10 PCs

# follow sandeep's pipeline.

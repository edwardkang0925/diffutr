library(SummarizedExperiment)
library(tidyverse)
library(dplyr)
library(diffUTR)
library(edgeR)


rds <- readRDS("outputs/exonExpression_hdl_preCovariatesAdj.RDS")
rds_preVST <- readRDS("/Users/test/projects/diffutr/outputs/exonExpression_hdl_preVST.RDS")
gtfPATH = "./data/gencode.v38.annotation.gtf"
bins <- prepareBins(gtfPATH, stranded=FALSE)

# Hot fix [reproduce the filter I used]
exp <- readRDS("outputs/exonExpression_hdl.RDS") # explore.R generates this file
exp_se <- SummarizedExperiment(assays=list(counts=exp),rowData=bins)
MinBinWidth <- 45 # it was chosen so that the number of the remaining bins < 1M
row_width_filter <- width(ranges(bins)) > MinBinWidth
exp_se <- exp_se[row_width_filter,]

maxAvgCount = 1000000
myRowMeans <- assays(exp_se)$counts %>% rowMeans(.)
row_excessive_filter <- myRowMeans < maxAvgCount
exp_se <- exp_se[row_excessive_filter,]


min_num_samples = floor(ncol(exp_se) * 0.015)
expression_level_filter <- rowSums(cpm(assays(exp_se)$counts > 3)) >= min_num_samples
exp_se <- exp_se[expression_level_filter,]


bins <- bins[row_width_filter,]
bins <- bins[row_excessive_filter,]
bins <- bins[expression_level_filter,]

# append bin size as rowData
bin_width <- width(ranges(bins))
rowData(rds)[, 'binWidth'] <- bin_width
rowData(rds)[,'avgCount'] <- assays(rds)[[1]] %>% rowMeans(.)
# bin size vs average count across subjects
plot(rowData(rds)[, 'binWidth'],  rowData(rds)[,'avgCount'], xlab='binWidth', ylab='avgNormCount', main = "post VST")

# Compare with preVST version of counts
rowData(rds_preVST)[, 'binWidth'] <- bin_width
rowData(rds_preVST)[,'avgCount'] <- assays(rds_preVST)$counts %>% rowMeans(.)
plot(rowData(rds_preVST)[, 'binWidth'],  rowData(rds_preVST)[,'avgCount'], xlab='binWidth', ylab='avgCount', main="pre VST")

# log(binwidth)
plot(log2(rowData(rds)[, 'binWidth']),  rowData(rds)[,'avgCount'], xlab='logBinWidth', ylab='avgNormCount', main = "log post VST")

# histogram of bin width distribution after filtering.
hist(bin_width, main="binWidth")

hist(log2(bin_width), main='logBinWidth')

# check the large bin sizes
largeBins <- bins[which(bin_width > 400000),]
counts <- assays(rds)[[1]]
counts[which(bin_width>400000),] %>% rowMeans(.)

.libPaths(c(.libPaths(), "/ref/mblab/software/r-envs/llfs/4.1"))  # package dependency.
setwd("/scratch/mblab/edwardkang/llfs_diffUTR")

################## Take care above lines when env issue is resolved currently using WOO's rstudio sbatch#######

############# current situation: R4.2 diffUTR to generate initial count file of 2000+bam files.
############# processing the massive data requires cluster, so here I am with R4.1

library(SummarizedExperiment)
library(tidyverse)
library(dplyr)
library(edgeR)
library(DESeq2)


#specify dependency files
combinedCount_rds = 'combined/diffUTR_count_combined.RDS'
phenotypeDataPATH = 'data/ADJUSTED_HEART_DISEASE_RELATED_TRAITS_FINAL_ROUND/adjusted_fhshdl.csv'
whatdatallPATH = 'data/whatdatall.csv'
qualimapPATH = 'data/chase/llfs_qualimap_20230323.csv'
outputPATH = 'outputs/'
binsPATH = 'data/diffUTR_bins_gencode_v38.RDS'
plotPATH = "outputs/plot_local/"

#
phenotype_df <- read.csv(phenotypeDataPATH)
whatdatall_df <- read.csv(whatdatallPATH)
qualimap_df <- read.csv(qualimapPATH)
rds <- readRDS(combinedCount_rds)
rds <- as.data.frame(rds)


# first from phenotype dataframe, extract subject id, use whatdatall to map to sample id to filter columns of rds
phenotype_subjectIDs <- phenotype_df[,1] # select subjectID the second column is the phenotype value

mappedSampleIDs = c()
unmappedSubjectIDs = c()
for (subjectID in phenotype_subjectIDs){
  sampleID = unique(whatdatall_df[which(whatdatall_df$subject == subjectID),]$id) # had to use 'which' to avoid mysterious NA rows. It happens when NA exists in the column
  if (length(sampleID) == 1){
    mappedSampleIDs = c(mappedSampleIDs, sampleID)
  }else{
    print("Unwanted case: subjectID -> sampleID not unique")
    print(paste0(subjectID, " -> ", sampleID))
    unmappedSubjectIDs <- c(unmappedSubjectIDs, subjectID)
  }
}

# using qualimap file, choose which of the duplicates to drop based on the percent_intergenic. (keep the lowest one)
qualimap_df_related <- qualimap_df[qualimap_df$id %in% mappedSampleIDs,]
qualimap_df_visit1_related <- qualimap_df_related[grepl('1', qualimap_df_related$visit),] # drop other than visit 1
qualimap_df_visit1_related_sorted <- qualimap_df_visit1_related[order(qualimap_df_visit1_related$percent_intergenic, decreasing=F),]
qualimap_df_visit1_related_sorted_unique <- qualimap_df_visit1_related_sorted %>% distinct(id, .keep_all = T)

# first handle the colnames to separate out samples. (visit 1 or 2 + whether in phenotype data)
filesToKeep = c()
for (bamfile in colnames(rds)){
  # Hard coded based on the filename format of diffUTR runs on 2561 bamfiles
  if (!grepl('pool',bamfile, fixed=T)){ # remove differently formatted bamfiles from the downstream analysis
    sampleID_visitcode <- strsplit(bamfile, split='[.]')[[1]][5]

    sampleID <- str_extract(sampleID_visitcode, '^\\d+')

    visitCode <- str_extract(sampleID_visitcode, "visit_\\d|visit\\d|vist\\d|vist_\\d")
    if(grepl('1', visitCode) & sampleID %in% mappedSampleIDs){ # visit 1 sample with phenotype
      visitCodeToUse = qualimap_df_visit1_related_sorted_unique[qualimap_df_visit1_related_sorted_unique$id == sampleID,c('visit')]
      if (visitCode == visitCodeToUse){
        filesToKeep = c(filesToKeep, bamfile)
      }
    }
    if (is.na(visitCode)){
      print(bamfile)
    }
  }
}

# drop unrelavent bamfiles. relavence determined by subjectID, visitcode pair for a trait file.
rds <- rds[, colnames(rds) %in% filesToKeep]


# rename the colnames. matching sandeep's format
subjects = c()
for (bamfile in colnames(rds)){
  sampleID_visitcode <- strsplit(bamfile, split='[.]')[[1]][5]

  sampleID <- str_extract(sampleID_visitcode, '^\\d+')
  subjectID <- unique(whatdatall_df[which(whatdatall_df$id == sampleID),]$subject)
  subject <- paste0('s', subjectID, '_v1')
  subjects = c(subjects, subject)
}
colnames(rds) <- subjects # rename colname from bamfilename to s<subjectID>_v<visitcode>

#saveRDS(rds, file=paste0(outputPATH, "exonExpression_hdl.RDS"))

# Re-combining bins information to exon count `rds`
bins <- readRDS(binsPATH)
exonCounts <- SummarizedExperiment(assays=list(counts=rds), rowRanges=bins)


# Filter for lowly expressed exons.
min_num_samples = floor(ncol(exonCounts) * 0.015)
expression_level_filter <- rowSums(cpm(assays(exonCounts)$counts > 3)) >= min_num_samples
#sum(expression_level_filter)
exonCounts_expressionFiltered <- exonCounts[expression_level_filter, ]

# filtering on bins based on width
minBinWidths <- c(1,3,5,10,15,20,25,30,35,40)
findProperBinWidth_hist <- function(minW, maxW, exonCounts_rse){
  rowWidthFilter <- ( width(ranges(rowRanges(exonCounts_rse))) >= minW &
                        width(ranges(rowRanges(exonCounts_rse))) <= maxW )
  return(sum(rowWidthFilter))
}
# for histogram, minbinW search while fixing maxW as 1000. Currently, graphics error
nrb <- lapply(minBinWidths,findProperBinWidth_hist, maxW=1000, exonCounts_rse=exonCounts_expressionFiltered )
data_hist <- data.frame(
  cbind(minBinWidths, nrb)
)
data_hist # console output for picking minW. continue with minW = 3 #being
saveRDS(data_hist, file=paste0(plotPATH, "minBinWidth_nrb.RDS"))

# actual filtering with bin Width begins.
minW = 10
maxW = 1000
rowWidthFilter <- ( width(ranges(rowRanges(exonCounts_expressionFiltered))) >= minW &
                      width(ranges(rowRanges(exonCounts_expressionFiltered))) <= maxW )
exonCounts_expressionFiltered <- exonCounts_expressionFiltered[rowWidthFilter, ]



# counts assay is currently a data.frame, need to convert into data.matrix
assays(exonCounts_expressionFiltered)$counts <- data.matrix(assays(exonCounts_expressionFiltered)$counts)

# VST transformation
# convert RangedSummarizedExperiment object to DESeqDataSet
dds <- DESeqDataSet(exonCounts_expressionFiltered, design=~1)
dds <- estimateSizeFactors(dds)
vst_obj <- vst(dds, blind=T)
assayNames(vst_obj) <- c("counts")
saveRDS(vst_obj, file=paste0(outputPATH, "vst_20230404_dds.RDS")) #output

# pc file generate
vst_count <- assay(vst_obj)
vst_count <- t(vst_count)
results <-  prcomp(vst_count, center = TRUE, scale = FALSE)
#variance explained
var_explained <- results$sdev^2 / sum(results$sdev^2) * 100
print(var_explained)

# output pc file
write.csv(results$x, paste0(outPATH,"pcs_20230404.csv"), quote = FALSE)










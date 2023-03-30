library(SummarizedExperiment)
library(tidyverse)
library(dplyr)


inputpath = '/Users/test/projects/diffutr/outputs/combined_allBams/diffUTR_count_combined.RDS'
phenotypeDataPATH = '/Users/test/projects/diffutr/data/ADJUSTED_HEART_DISEASE_RELATED_TRAITS_FINAL_ROUND/adjusted_fhshdl.csv'
whatdatallPATH = '/Users/test/projects/diffutr/data/whatdatall.csv'
qualimapPATH = '/Users/test/projects/diffutr/data/chase/llfs_qualimap_20230323.csv'

phenotype_df <- read.csv(phenotypeDataPATH)
whatdatall_df <- read.csv(whatdatallPATH)
qualimap_df <- read.csv(qualimapPATH)
rds <- readRDS(inputpath)
rds <- as.data.frame(rds)

phenotype_subjectIDs <- phenotype_df[,1]


# first from phenotype dataframe, extract subject id, use whatdatall to map to sample id to filter columns of rds
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

# drop unrelavent bamfiles.
rds <- rds[, colnames(rds) %in% filesToKeep]
saveRDS(rds, file=paste0("./outputs/", "exonExpression_hdl.RDS"))

# preprocessing (handling rownames) + VST



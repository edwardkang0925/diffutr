library(SummarizedExperiment)
library(dplyr)
library(stringr)

whatdatallPATH = '/Users/test/projects/diffutr/data/whatdatall.csv'
pcPATH = '/Users/test/projects/diffutr/outputs/pcs_exonExpression.csv'
vstPATH = '/Users/test/projects/diffutr/outputs/vst_normalized.csv'

whatdatall_df <- read.csv(whatdatallPATH)
pc_df <- read.csv(pcPATH)
#vst_df <- read.csv(vstPATH)

subjects = c()
for (bamfile in pc_df$X){
  sampleID_visitcode <- strsplit(bamfile, split='[.]')[[1]][5]

  sampleID <- str_extract(sampleID_visitcode, '^\\d+')
  subjectID <- unique(whatdatall_df[which(whatdatall_df$id == sampleID),]$subject)
  subject <- paste0('s', subjectID, '_v1')
  subjects = c(subjects, subject)
}

pc_df$X <- subjects
colnames(vst_df) <- subjects

write.csv(pc_df, "outputs/pcs_exonExpression_renamed.csv", quote=FALSE)
write.csv(vst_df, "outputs/vst_normalized_renamed.csv", quote = FALSE)

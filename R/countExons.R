library(diffUTR)
library(SummarizedExperiment)
library(optparse)
library(foreach)
library(tidyverse)

bamfiles = c('data/10367282_visit_1.markdup.sorted.bam',
             'data/10367282_visit_2.markdup.sorted.bam',
             'data/10512482_visit_1.markdup.sorted.bam',
             'data/10512482_visit_2.markdup.sorted.bam',
             'data/10698322_visit_1.markdup.sorted.bam',
             'data/10698322_visit_2.markdup.sorted.bam',
             'data/11138864_visit_1.markdup.sorted.bam',
             'data/11138864_visit_2.markdup.sorted.bam',
             'data/11289212_visit_1.markdup.sorted.bam',
             'data/11289212_visit_2.markdup.sorted.bam',
             'data/45213600_visit1.markdup.sorted.bam')
gtfFile = './data/gencode.v38.annotation.gtf'

whatdatall <- read.csv("data/whatdatall.csv")


bins <- prepareBins(gtfFile, stranded=FALSE)

rse <- countFeatures(bamfiles, bins, strandSpecific=0, isPairedEnd=TRUE,  GTF.featureType='exon', GTF.attrType='exon_id', nthreads=6)

# extract sampleID from the colnames (bamfile names) and map to subjectID
subject_id = c()
for (sample_id in str_extract(colnames(rse), "^\\d+") ){
  subject_id = c(subject_id, subset(whatdatall, id == sample_id)$subject[1])
}
subject_id


# extract visit
str_extract(colnames(rse), 'visit_\\d|visit\\d|vist\\d')


#save the SE object from bamfile
saveRDS(rse, file=paste0("./outputs/merged11Bam_counts_droppedDup.RDS"))




















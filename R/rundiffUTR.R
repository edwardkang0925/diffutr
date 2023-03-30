#!/usr/bin/env Rscript

library(diffUTR)
library(SummarizedExperiment)
library(optparse)
library(foreach)

TEST = TRUE

parser <- OptionParser()
#parser <- add_option(parser, c("--test"), action="store_true",
#                     default=FALSE, help="Activate test run")
parser <- add_option(parser, c("--gtf"), type="character",
                     default="",
                     help="path to gtf file",
                     metavar="'/path/to/file.gtf'")
parser <- add_option(parser, c("--bam"), type="character",
                     default="",
                     help="path to bam file",
                     metavar="'/path/to/file.bam'")
parser <- add_option(parser, c("--output_prefix"), type="character",
                     default="./outputs/",
                     help="output location",
                     metavar="/path/to/output/dir")
parser <- add_option(parser, c("--jobid"), type="integer",
                     default="1",
                     help="Slurm array job id",
                     metavar="some integer N")
# main -------------------------------------------------------------------------

opt = if(TEST){
  parse_args(parser,
             args = c(
               "--gtf=./data/gencode.v38.annotation.gtf",
               "--bam=./data/10367282_visit_2.markdup.sorted.bam"))
} else{
  parse_args(parser)
}

# check cmd line arguments
x = foreach(
  i = names(opt)
) %do% {
  input_value=opt[[i]]
  if(input_value==''){
    stop(sprintf("ARGUMENT --%s IS REQUIRED",i))
  }else if(i %in% c('dna','rna')){
    if(!file.exists(input_value)){
      stop(sprintf("FILE %s DOES NOT EXIST",input_value))
    }
  }
}

# llfs reads are unstranded, so no need to use stranded=TRUE
bins <- prepareBins(opt$gtf, stranded=FALSE)

# reading a single bam file
bamfiles <- c(opt$bam)

# strandSpecific=0 for unstranded data. GTF.attrType defines which feature level to groupby
rse <- countFeatures(bamfiles, bins, strandSpecific=0, isPairedEnd=TRUE,  GTF.featureType='exon', GTF.attrType='exon_id')

# output SE
saveRDS(rse, file=paste0(opt$output_prefix, "diffUTR_count_", opt$jobid, ".RDS"))

## EPA

count_assay <- assays(rse)$counts

count_assay['ENSG00000001626.16.99',]
c(rowData(rse)[1,'exon_id'][[1]], rowData(rse)[2,'exon_id'][[1]],rowData(rse)[3,'exon_id'][[1]])

exon_list = c()
for (binName in rownames(rse)){
  exons <- rowData(rse)[binName,'exon_id'][[1]]
  exon_list = c(exon_list, exon)
}

for (binName in rownames(rse)){
  exons <- rowData(rse)[binName,'exon_id'][[1]]
  for (exon in exons){
    count <- count_assay[binName,]
    print(paste(exon, " ", count))
    # currently there are duplicated exon id, sum up by exon_id
  }
}






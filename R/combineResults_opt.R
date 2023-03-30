#!/usr/bin/env Rscript

library(SummarizedExperiment)
library(optparse)

parser <- OptionParser()
#parser <- add_option(parser, c("--test"), action="store_true",
#                     default=FALSE, help="Activate test run")
parser <- add_option(parser, c("--input_prefix"), type="character",
                     default="./outputs/",
                     help="location of RDS files to combine",
                     metavar="/path/to/dir/with/RDSfiles")

parser <- add_option(parser, c("--output_prefix"), type="character",
                     default="./outputs/",
                     help="output location",
                     metavar="/path/to/output/dir")


# main -------------------------------------------------------------------------

opt = parse_args(parser)

filenames <- list.files(opt$input_prefix, pattern="*.RDS", full.names=TRUE)

combined_obj <- readRDS(filenames[1])
combined_obj <- assays(combined_obj)$counts

for (filename in filenames[c(2:length(filenames))]){
  print(filename)
  new_obj <- readRDS(filename)
  new_obj <- assays(new_obj)$counts
  combined_obj <- cbind(combined_obj, new_obj)
}

saveRDS(combined_obj, file=paste0(opt$output_prefix, "RSubRead_count_combined", ".RDS"))

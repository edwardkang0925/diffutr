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
rdsobjs <- lapply(filenames, readRDS) # FIXME: probably lthis is computationally horrible.
combined_objs <- do.call(cbind, rdsobjs)

saveRDS(combined_objs, file=paste0(opt$output_prefix, "RSubRead_count_combined", ".RDS"))

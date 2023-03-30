library(Rsubread)

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
               "--bam=./data/45213600_visit1.markdup.sorted.bam"))
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



exon_counts_bothendsmapped_igdup <- featureCounts(files=opt$bam, annot.ext=opt$gtf,
                                                  isGTFAnnotationFile=TRUE, GTF.featureType='exon',
                                                  GTF.attrType="exon_id",
                                                  useMetaFeatures=FALSE, allowMultiOverlap=TRUE,
                                                  strandSpecific = 0, isPairedEnd = TRUE, nthreads=1,
                                                  requireBothEndsMapped=TRUE, ignoreDup = TRUE)
#drop duplicated rows
exon_counts_bothendsmapped_igdup$counts <- unique(exon_counts_bothendsmapped_igdup$counts)
#save only the counts table
saveRDS(exon_counts_bothendsmapped_igdup$counts, file=paste0(opt$output_prefix, "rSubRead_count_", opt$jobid, ".RDS"))




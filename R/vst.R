library(DESeq2)


# read saved RDS obj and outputs its assay as csv
dds <- readRDS("outputs/exonExpression_hdl_preCovariatesAdj.RDS")

vst_expr = assays(dds)
write.csv(vst_expr, "outputs/vst_normalized.csv", quote = FALSE)

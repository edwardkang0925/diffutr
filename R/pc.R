library(ggplot2)
library(DESeq2)

vst <- readRDS("outputs/exonExpression_hdl_preCovariatesAdj.RDS")

vst_expr <- assay(vst)
# transpose
vst_expr <- t(vst_expr)

# pc
results <- prcomp(vst_expr, center = TRUE, scale = FALSE)

#variance explained
var_explained <- results$sdev^2 / sum(results$sdev^2) * 100
print(var_explained)

#Write PCs into a matrix
write.csv(results$x, "outputs/pcs_exonExpression.csv", quote = FALSE)

# args <- commandArgs(trailingOnly = TRUE)
library(tidyverse)
library(MASS)
library(dplyr)
library(here)
library(stringr)
library(fastDummies)
# library(vsn)
library(hexbin)
library(ggplot2)
library(reshape2)

# misc_file_path <- args[1]
# print(misc_file_path)
# output_file_path <- args[2]
# print(output_file_path)
# phenotype_file_path <- args[3]
# print(phenotype_file_path)
# phenotype <- args[4]
# print(phenotype)
# number_of_plates <- args[5]
# print(number_of_plates)
# visit_type <- args[6]
# print(visit_type)

# manually defined parameters for testing purpose
misc_file_path="data/MISC_FINAL_ROUND/"
output_file_path="outputs/"
phenotype_file_path="data/ADJUSTED_HEART_DISEASE_RELATED_TRAITS_FINAL_ROUND/"
phenotype="adjusted_fhshdl"
number_of_plates=28
visit_type="visit_1"

blood_all_filepath <- paste(misc_file_path, "blood_all.csv", sep = "")
pcs_filepath <- paste(output_file_path, "pcs_exonExpression_renamed.csv", sep = "")
metadata_filepath <- paste(misc_file_path, "sample_meta_20220608.csv", sep = "")
subject_age_fc_filepath <- paste(misc_file_path, "subject_fc_sex_age_revised_no_blood_cancer.csv", sep = "")
phenotype_filepath <- paste(phenotype_file_path, phenotype,".csv", sep = "")
exon_expression_filepath <- paste(output_file_path, "vst_normalized_renamed.csv", sep = "")

#reading and pre-processing_files

blood_all_df <- read.table(file = blood_all_filepath,  header = TRUE, sep = ",")
print(dim(blood_all_df))
pcs_df <- read.table(file = pcs_filepath, header = TRUE, sep = ",")
print(dim(pcs_df))
metadata_df <- read.table(file = metadata_filepath, header = TRUE, sep = ",")
print(dim(metadata_df))
subject_age_fc_df <- read.table(file = subject_age_fc_filepath, header = TRUE, sep = ",")
print(dim(subject_age_fc_df))
phenotype_df <- read.table(file = phenotype_filepath, header = TRUE, sep = ",")
print(dim(phenotype_df))
exon_expression_df <- read.table(file = exon_expression_filepath, header = TRUE, sep = ",")
print(dim(exon_expression_df))
#remove a row with fc null input
subject_age_fc_df <- subject_age_fc_df[subject_age_fc_df$fc != "", ]

#parse pcs_filepath
#EDIT: VISIT
pcs_df <- pcs_df[, (names(pcs_df) %in% c("X", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))]
pcs_df <- pcs_df %>%
  filter(str_detect(X, "v1"))
#parse subject in pcs_df
pcs_df <- pcs_df %>%
  separate(X, into = c("subjecttemp", "visit"), sep = "_")
pcs_df <- pcs_df %>%
  separate(subjecttemp, into = c("rem", "subject"), sep = "s")
pcs_df$rem <- NULL
pcs_df$visit <- NULL

#parse metadata_df
metadata_df <- dummy_cols(metadata_df, select_columns = "plate", remove_first_dummy = TRUE, remove_selected_columns = TRUE)
plates <- list()
for (i in c(2:number_of_plates)){
  plates <- append(plates, paste("plate_", as.character(i), sep = ""))
}

#EDIT:VISIT
metadata_df <- metadata_df %>%
  filter(str_detect(subject_count_headers, "v1"))
covariates_metadata_df <- append(plates, c("subject", "percent_intergenic"))
metadata_df <- metadata_df[, (names(metadata_df) %in% covariates_metadata_df)]
metadata_df <- metadata_df[metadata_df$subject %in% pcs_df$subject,]

#parse subject_age_fc_df
subject_age_fc_df <- dummy_cols(subject_age_fc_df, select_columns = "fc", remove_first_dummy = TRUE, remove_selected_columns = TRUE)
subject_age_fc_df <- subject_age_fc_df[subject_age_fc_df$visitcode == 1,]
#parse blood_all_df
blood_all_df <- blood_all_df %>%
  filter(str_detect(visitcode, "Visit 1"))
blood_all_df <- blood_all_df[, (names(blood_all_df) %in% c("subject", "amon", "aneu", "rbc","wbc", "plt"))]

#change subject id to numeric
pcs_df$subject <- as.numeric(as.character(pcs_df$subject))
metadata_df$subject <- as.numeric(as.character(metadata_df$subject))
blood_all_df$subject <- as.numeric(as.character(blood_all_df$subject))
subject_age_fc_df$subject <- as.numeric(as.character(subject_age_fc_df$subject))
phenotype_df$subject <- as.numeric(as.character(phenotype_df$subject))

#input to phenotype adjustments. We only want subjects that are in subject_age_fc_df, blood_all_df, metadata_df, pcs_df
# temp_list <- list(blood_all_df, metadata_df, pcs_df, subject_age_fc_df)
# temp_df <- temp_list %>% purrr::reduce(full_join, by = "subject")
# temp_df[temp_df==""]<-NA
# temp_df <- temp_df %>% drop_na()
# write.csv(temp_df, paste(output_file_path, "input_phenotype_adjustment_",visit_type,".csv", sep = ""), row.names = FALSE, quote = FALSE)

#combine dfs
combined_df_list <- list(phenotype_df,subject_age_fc_df,blood_all_df, metadata_df, pcs_df)
combined_df <- combined_df_list %>% purrr::reduce(full_join, by = "subject")
combined_df[combined_df==""]<-NA
combined_df <- combined_df %>% drop_na()

#change blood count to numeric and add age squared
combined_df$plt <- as.numeric(as.character(combined_df$plt))
combined_df$rbc <- as.numeric(as.character(combined_df$rbc))
combined_df$wbc <- as.numeric(as.character(combined_df$wbc))
rownames(combined_df) <- combined_df$subject
combined_df$subject <- NULL
combined_df <- combined_df[ order(as.numeric(row.names(combined_df))), ]
combined_df$age_2 <- combined_df$age ^ 2
print("Do first and second match? The first is the phenotype df. The second is the combined df.")
print(dim(phenotype_df))
print(dim(combined_df))

# <----------------- combined df dimension should be same as the dimension of phenotype df if Lisa and Sandeep's data processing is accurate.--------------->

#remove phenotype from combined df
combined_df <- combined_df[-c(1)]

# < --------------------------- RESIDUAL GENERATION ------------------------>
# parse gene expression data
exon_expression_df <- t(exon_expression_df)
colnames(exon_expression_df) <- exon_expression_df[1,]
exon_expression_df <- exon_expression_df[-c(1),]
#EDIT: VISIT
exon_expression_df <- exon_expression_df[grepl("v1", rownames(exon_expression_df)),]

split_level_1 <- strsplit(rownames(exon_expression_df), "_")
split_level_1_grep_list <- lapply(split_level_1, `[[`, 1)
subject_ids <- as.numeric(lapply(split_level_1_grep_list, function(x) gsub("s", " ", x, fixed=TRUE)))
rownames(exon_expression_df) <- subject_ids

#make sure residuals were calculated for samples for which we have genotype pcs, phenotype, and other covariates
exon_expression_df <- exon_expression_df[rownames(exon_expression_df) %in% rownames(combined_df),]
exon_expression_df <- exon_expression_df[ order(as.numeric(row.names(exon_expression_df))), ]
combined_df <- combined_df[ order(as.numeric(row.names(combined_df))), ]
combined_df$visitcode <- NULL

print("The dimension of the gene expression matrix and the combined df should match.")
print(dim(exon_expression_df))

residuals_df = data.frame(matrix(ncol = dim(exon_expression_df)[2] , nrow = dim(exon_expression_df)[1]))
rownames(residuals_df) <- rownames(combined_df)
colnames(residuals_df) <- colnames(exon_expression_df)
index <- 1

#check if residuals are normal. Create list of p values and w values for ben shapiro test.
residuals_normality_shapiro_p_values <- vector(mode='list', length=dim(exon_expression_df)[2])
residuals_normality_shapiro_w_values <- vector(mode='list', length=dim(exon_expression_df)[2])
base_covariates_list = c("age", "sex", "fc_DK", "fc_NY", "fc_PT", "amon", "aneu", "plt", "rbc", "wbc", "percent_intergenic", "age_2")
additional_covariates_list = colnames(combined_df)[12:48]
temp_exon_expression_df <-  exon_expression_df[,1:2]
#HERE
count = 0
#residual adjustment using stepwise regression
for (values in colnames(exon_expression_df)){
  #new code start
  count = count + 1
  if(count %% 2 == 0){
    print(count)
  }
  combined_df["gene"] <- as.numeric(exon_expression_df[,values])
  fixed_effects <- as.formula(paste("gene", paste(c(base_covariates_list, additional_covariates_list), collapse="+"), sep="~"))
  model <- lm(fixed_effects, data=combined_df)
  model_step <- step(model, scope=list(upper = as.formula(paste("~", paste(c(base_covariates_list, additional_covariates_list), collapse="+"))),
                                       lower = as.formula(paste("~", paste(base_covariates_list, collapse="+")))), trace=FALSE)
  temp_residuals <- resid(model_step)
  residuals_df[,index] <- temp_residuals
  #new code end
  combined_df["gene"] <- NULL
  residuals_normality_shapiro_p_values[[index]] <- shapiro.test(temp_residuals)$p.value
  residuals_normality_shapiro_w_values[[index]] <- shapiro.test(temp_residuals)$W
  index <- index + 1
}

print(summary(model_step))
negative_log_residuals_normality_shapiro_p_values <- as.numeric(lapply(residuals_normality_shapiro_p_values, function(x) -log10(x)))

# <--------- GENERATE MMAP INPUT -------------- >
residuals_df$EGO <- rownames(residuals_df)

rownames(phenotype_df) <- phenotype_df$subject
phenotype_df$subject <- NULL
phenotype_df$EGO <- rownames(phenotype_df)

phenotype_df$EGO <- as.numeric(as.character(phenotype_df$EGO))
residuals_df$EGO <- as.numeric(as.character(residuals_df$EGO))

residuals_plus_covariates_list <- list(phenotype_df, residuals_df)
residuals_plus_covariates <- residuals_plus_covariates_list %>% purrr::reduce(full_join, by = "EGO")
residuals_plus_covariates[residuals_plus_covariates==""]<-NA
residuals_plus_covariates <- residuals_plus_covariates %>% drop_na()
residuals_plus_covariates <- residuals_plus_covariates %>%
  relocate(EGO)

#write twas residuals/MMAP input in the output directory
print(dim(residuals_plus_covariates))
write.csv(residuals_plus_covariates, paste(output_file_path,"twas_input_",visit_type,"_",phenotype,".csv", sep = ""), row.names = FALSE, quote = FALSE)

#histogram for the distribution of p-values along with the number of genes that passed the normality test.
# jpeg(file = paste(output_file_path, "Figures/", "twas_log10ShapiroTest_",visit_type, phenotype, ".jpeg", sep = ""))
normality_test_passed_genes <- sum(negative_log_residuals_normality_shapiro_p_values < 1.30102999566)
print(normality_test_passed_genes)
# hist(negative_log_residuals_normality_shapiro_p_values,
#      main = paste("Shapiro Test - ", as.character(normality_test_passed_genes),"genes"),
#      xlab="-log10p",
#      col="darkmagenta",
#      breaks=50,
#      ylim=c(0,7000))
# dev.off()

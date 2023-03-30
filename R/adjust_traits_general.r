
###########################################################
#           User inputs                             
###########################################################

args = commandArgs(trailingOnly=TRUE)

master <- read.csv(args[1]) # csv containing covariate values for everyone
trait_covar <- read.csv(args[2]) # csv containing the names of the covariates for each trait (no genetic PCs), dimension = trait x covariates, see example for what the headers should be

trait_input <- args[3] # This must be the same as the column label in the master file, and must be the same as the column name in the traits_covar file

output_dir <- args[4]

dir.create(file.path(output_dir), showWarnings = FALSE)
dir.create(file.path(output_dir, "stepwise_models"), showWarnings = FALSE)
dir.create(file.path(output_dir, "residual_hist"), showWarnings = FALSE)
dir.create(file.path(output_dir, "residual_qq"), showWarnings = FALSE)
dir.create(file.path(output_dir, "raw_value_hist"), showWarnings = FALSE)
dir.create(file.path(output_dir, "raw_value_qq"), showWarnings = FALSE)
print(args[2])
print(trait_covar$trait)


adjust_trait <- function(trait, df_trait_covar, df_master){

    # expression covariates
    expr_covar <- c('amon', 'aneu', 'plt', 'rbc', 'wbc', 'percent_intergenic', 'ePC1', 'ePC2', 'ePC3')

    # parsing the csv to get a list of covariates for the current trait
    covar_list <- as.character(df_trait_covar[df_trait_covar$trait == trait, ])
    covar_list_lower <- covar_list[nzchar(covar_list)]
    covar_list_lower <- covar_list_lower[covar_list_lower != 'fc']

    covar_list_upper <- c(covar_list_lower, paste0('PC', seq(1, 10))) # append PC1-PC10 to the covariates list 

    # keep only samples (rows) without any NAs  
    samples_complete <- df_master[complete.cases(df_master[, c(covar_list_upper, expr_covar)]), covar_list_upper]
    samples_complete_id <- df_master[complete.cases(df_master[, c(covar_list_upper, expr_covar)]), 'subject']
    
    # output raw value histograms and normal qq plots for complete samples 
    png(paste0(output_dir, "/raw_value_hist/raw_value_hist_", trait, ".png"))
    hist(samples_complete[, trait], main=paste0("raw ", trait))
    dev.off()

    png(paste0(output_dir, "/raw_value_qq/raw_value_qq_", trait, ".png"))
    min_y <- ceiling(min(samples_complete[, trait])) + 1
    observed <- sort(samples_complete[, trait])
    qqnorm(observed, main=paste0("raw ", trait))
    qqline(observed)
    dev.off()
    
    # creates a formula for the minimum model used in the stepwise regression, the formula looks like this: trait ~ age + age^2 + sex + DK + BC + NY
    fixed_effects <- as.formula(paste(trait, paste(covar_list_upper[2:length(covar_list_upper)], collapse="+"), sep="~"))

    model <- lm(fixed_effects, data=samples_complete)

    # all covariates in the fixed_effects formula will always be adjusted, while PC1-PC10 is added in each step
    model_step <- step(model, scope=list(upper = as.formula(paste("~", paste(covar_list_upper[2:length(covar_list_upper)], collapse="+"))), 
                                lower = as.formula(paste("~", paste(covar_list_lower[2:length(covar_list_lower)], collapse="+")))), trace=FALSE)
    print(summary(model_step))

    residuals <- resid(model_step) # get resuduals after stepwise
    
    # output residual histograms and normal qq plots
    png(paste0(output_dir, "/residual_hist/residual_hist_", trait, ".png"))
    hist(residuals, main=paste0("residual ", trait))
    dev.off()
    png(paste0(output_dir, "/residual_qq/residual_qq_", trait, ".png"))
    min_y <- ceiling(min(residuals)) + 1
    observed <- sort(residuals)
    qqnorm(observed, main=paste0("residual ", trait))
    qqline(observed)
    dev.off()

    alpha <- 0.05

    # shapiro test for normality, if p <= 0.05, the residual distribution is significantly different from normal, so we inverse normal transform them
    residual_shapiro <- shapiro.test(residuals)
    normal <- TRUE
    if(residual_shapiro$p.value <= 0.05){
        # not normal, inverse normal transform residuals
        # source: https://www.biostars.org/p/80597/
        print(paste0(trait, " trait not normal, INT"))
        residuals <- qnorm((rank(residuals,na.last="keep")-0.5)/sum(!is.na(residuals)))
        png(paste0(output_dir, "/residual_hist/residual_hist_int_", trait, ".png"))
        hist(residuals,main=paste0("residual INT ", trait))
        dev.off()
        normal <- FALSE
    }

    
    df_residuals <- data.frame(subject=samples_complete_id)
    df_residuals[, trait] <- residuals
    
    # saving the included variables in stepwise regression to csv. 
    final_params <- variable.names(model_step)
    final_params <- final_params[2:length(final_params)] # exclude intercept

    # save shapiro test summary
    if(normal){
        normality_check <- rep(NA, length(final_params))
    }
    else{
        print(residual_shapiro$p.values)
        normality_check <- c(residual_shapiro$p.value, rep(NA, length(final_params)-1))
   
    }

    df_model_meta <- data.frame(final_params=final_params, normality_check=normality_check)

    # final residuals, inverse normal transformed if it failed the shapiro test
    write.csv(df_residuals, paste0(output_dir, "/adjusted_", trait, ".csv"), row.names=FALSE) 
    # stepwise regression model summary
    write.csv(df_model_meta, paste0(output_dir, "/stepwise_models/stepwise_model_", trait, ".csv"), row.names=FALSE)
}

###########################################################
#           Main                             
###########################################################



if(trait_input == 'all'){
    for(t in trait_covar$trait){
        print(paste0("Trait to adjust : ", t))

        adjust_trait(t, trait_covar, master)
    }
} else{
    print(paste0("Trait to adjust : ", trait_input))

    adjust_trait(trait_input, trait_covar, master)
}

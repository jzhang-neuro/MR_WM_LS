# List of risk factors to analyze
risk_factors <- c("Hypertension", "Type2Diabetes", "Lipoprotein_Metabolism")

# List of all IDP traits (assuming you have 288 dMRI traits files)
idp_files <- list.files(pattern = "exp_dat_IDP_.*\\.csv")

# Initialize dataframes to store final results
significant_results <- data.frame()
heterogeneity_results <- data.frame()
pleiotropy_results <- data.frame()
directionality_results <- data.frame()

# Set significance threshold (adjust this if necessary)
significance_threshold <- 0.05 / (length(risk_factors) * length(idp_files))  # Bonferroni correction

# Loop through each risk factor
for (risk_factor in risk_factors) {
  
  # Load preprocessed data for the current risk factor
  risk_factor_file <- paste0("exp_dat_", risk_factor, ".csv")
  risk_dat <- read.csv(risk_factor_file)
  
  # Loop through each IDP trait file
  for (i in seq_along(idp_files)) {
    idp_file <- idp_files[i]
    
    # Load preprocessed data for the current IDP trait
    idp_dat <- read.csv(idp_file)
    
    # Extract the IDP name from the file name
    idp_name <- sub(".*exp_dat_IDP_(.*)\\.csv", "\\1", idp_file)
    
    # Harmonize exposure (risk factor) and outcome (IDP trait) data
    dat <- harmonise_data(
      exposure_dat = risk_dat,
      outcome_dat = idp_dat
    )
    
    # Filter the data for valid instruments (where mr_keep == TRUE)
    valid_dat <- dat[dat$mr_keep == TRUE, ]
    
    # Check if there are enough valid SNPs after filtering
    if (is.null(valid_dat) || nrow(valid_dat) < 2) {
      cat("Not enough SNPs in valid instruments for MR analysis between", risk_factor, "and", idp_name, "\n")
      next
    }
    
    # Perform MR analysis
    res <- mr(valid_dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_mode", "mr_weighted_median"))
    
    
    # Extract IVW method p-value
    ivw_pval <- res[res$method == "Inverse variance weighted", "pval"]
    
    # Perform MR-RAPS analysis
    mr_raps_result <- mr.raps(valid_dat$beta.exposure, valid_dat$beta.outcome, valid_dat$se.exposure, valid_dat$se.outcome)
    
    # Organize MR-RAPS results
    mr_raps_df <- data.frame(
      id.exposure = unique(valid_dat$id.exposure),
      id.outcome = unique(valid_dat$id.outcome),
      outcome = unique(valid_dat$outcome),
      exposure = unique(valid_dat$exposure),
      method = "mr_raps",
      nsnp = nrow(valid_dat),
      b = mr_raps_result$beta.hat,
      se = mr_raps_result$beta.se,
      pval = mr_raps_result$beta.p.value
    )
    
    # Combine MR and MR-RAPS results
    res <- rbind(res, mr_raps_df)
    
    # Perform MR-PRESSO outlier detection if enough instruments are available
    if (nrow(valid_dat) >= 2) {
      tryCatch({
        presso_result <- run_mr_presso(dat = valid_dat)
        
        # Create dataframes for MR-PRESSO results
        mr_presso_df <- data.frame(
          analysis = c("Raw", "Outlier-corrected"),
          causal_estimate = presso_result[[1]]$`Main MR results`$`Causal Estimate`,
          sd = presso_result[[1]]$`Main MR results`$Sd,
          t_stat = presso_result[[1]]$`Main MR results`$`T-stat`,
          p_value = presso_result[[1]]$`Main MR results`$`P-value`
        )
        
        global_test_df <- data.frame(
          rss_obs = presso_result[[1]]$`MR-PRESSO results`$`Global Test`$RSSobs,
          p_value = presso_result[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
        )
        
        outlier_test_df <- presso_result[[1]]$`MR-PRESSO results`$`Outlier Test`
        
        distortion_test_df <- data.frame(
          outliers_indices = presso_result[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`,
          distortion_coefficient = presso_result[[1]]$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`,
          p_value = presso_result[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue
        )
        
        # Combine all MR-PRESSO results into a list
        result_list <- list(
          MR_PRESSO_results = mr_presso_df,
          MR_PRESSO_global_test = global_test_df,
          MR_PRESSO_outlier_test = outlier_test_df,
          MR_PRESSO_distortion_test = distortion_test_df
        )
        
        # Save each dataframe in the result list to a CSV file
        file_prefix <- gsub(" ", "_", idp_name)
        risk_factor_prefix <- gsub(" ", "_", risk_factor)  # Ensure risk factor name has no spaces
        
        for (name in names(result_list)) {
          write.csv(result_list[[name]], 
                    file = paste0(i, "_", name, "_", file_prefix, "_", risk_factor_prefix, ".csv"), 
                    row.names = FALSE)
        }
        
        # Extract Outlier-corrected results
        main_mr_results <- presso_result[[1]]$`Main MR results`
        outlier_corrected <- main_mr_results[main_mr_results$`MR Analysis` == 'Outlier-corrected', ]
        
        causal_estimate <- outlier_corrected$`Causal Estimate`
        sd <- outlier_corrected$Sd
        p_value <- outlier_corrected$`P-value`
        
        # Extract p-values from the Outlier Test
        outlier_test <- presso_result[[1]]$`MR-PRESSO results`$`Outlier Test`
        outlier_test$Pvalue <- as.numeric(gsub("[^0-9.]", "", outlier_test$Pvalue))
        
        # Calculate the number of outliers with Pvalue < 0.05
        nsnp_outlier <- sum(outlier_test$Pvalue < 0.05, na.rm = TRUE)
        
        # Calculate the effective number of instruments
        nsnp_effective <- nrow(valid_dat) - nsnp_outlier
        
        # Create MR-PRESSO Outlier-corrected results dataframe
        mr_presso_outlier_corrected <- data.frame(
          id.exposure = unique(valid_dat$id.exposure),
          id.outcome = unique(valid_dat$id.outcome),
          outcome = unique(valid_dat$outcome),
          exposure = unique(valid_dat$exposure),
          method = "MR-PRESSO_Outlier_corrected",
          nsnp = nsnp_effective,
          b = causal_estimate,
          se = sd,
          pval = p_value
        )
        
        # Add MR-PRESSO Outlier-corrected results to res dataframe
        res <- rbind(res, mr_presso_outlier_corrected)
        
      }, error = function(e) {
        message("Error in MR-PRESSO analysis: ", e$message)
      })
    } else {
      message("Not enough instrumental variables for MR-PRESSO analysis. Skipping this step.")
    }
    
    # Generate odds ratios
    or_df <- generate_odds_ratios(res)
    
    # Merge results and odds ratios
    result <- merge(res, or_df, by = c("id.exposure", "id.outcome", "outcome", "exposure", "method", "nsnp", "b", "se", "pval"))
    
    # Calculate overall R^2 and F-statistic
    overall_R2 <- sum(valid_dat$R2)
    overall_R2.easy <- sum(valid_dat$R2.easy)
    
    k <- nrow(valid_dat)  # Number of instruments
    overall_F_statistic <- (overall_R2 * (N - k - 1)) / (k * (1 - overall_R2))
    
    # Add overall R^2 and F-statistic to result
    result$overall_R2 <- overall_R2
    result$overall_F_statistic <- overall_F_statistic
    
    # Save MR results to a CSV file
    output_file <- paste0("res_", risk_factor, "_on_", idp_name, ".csv")
    write.csv(result, file = output_file, row.names = FALSE)
    
    # Check if any "Inverse variance weighted" method p-values are significant
    if (any(res$method == "Inverse variance weighted" & res$pval < significance_threshold)) {
      res$risk_factor <- risk_factor  # Add risk factor name column
      res$idp_trait <- idp_name  # Add IDP trait column
      res$i <- i  # Add loop index column
      significant_results <- rbind(significant_results, res)
    }
    
    # Visualization: scatter plot
    p1 <- mr_scatter_plot(res, valid_dat)
    ggsave(p1[[1]], file = paste0("scatter_", risk_factor, "_on_", idp_name, ".pdf"), width = 7, height = 7)
    
    # Single SNP analysis: forest plot
    res_single <- mr_singlesnp(valid_dat)
    p2 <- mr_forest_plot(res_single)
    ggsave(p2[[1]], file = paste0("forest_", risk_factor, "_on_", idp_name, ".pdf"), width = 7, height = 7)
    
    # Leave-one-out analysis: leave-one-out plot
    res_loo <- mr_leaveoneout(valid_dat)
    p3 <- mr_leaveoneout_plot(res_loo)
    ggsave(p3[[1]], file = paste0("loo_", risk_factor, "_on_", idp_name, ".pdf"), width = 7, height = 7)
    
    # Funnel plot
    p4 <- mr_funnel_plot(res_single)
    ggsave(p4[[1]], file = paste0("funnel_", risk_factor, "_on_", idp_name, ".pdf"), width = 7, height = 7)
    
    # Additional analyses:
    # Heterogeneity test
    het <- mr_heterogeneity(valid_dat)
    het$risk_factor <- risk_factor
    het$idp_trait <- idp_name
    het$i <- i  # Add loop index column
    heterogeneity_results <- rbind(heterogeneity_results, het)
    
    # Pleiotropy test
    pleio <- mr_pleiotropy_test(valid_dat)
    pleio$risk_factor <- risk_factor
    pleio$idp_trait <- idp_name
    pleio$i <- i  # Add loop index column
    pleiotropy_results <- rbind(pleiotropy_results, pleio)
    
    # Add loop index column, risk factor name, ivw_pval to Steiger directionality test
    steiger_df$risk_factor <- risk_factor
    steiger_df$idp_trait <- idp_name
    steiger_df$ivw_pval <- ivw_pval
    steiger_df$i <- i
    directionality_results <- rbind(directionality_results, steiger_df)
  }
}

# Save significant results to a CSV file
write.csv(significant_results, "significant_results.csv", row.names = FALSE)
write.csv(heterogeneity_results, "heterogeneity_results.csv", row.names = FALSE)
write.csv(pleiotropy_results, "pleiotropy_results.csv", row.names = FALSE)
write.csv(directionality_results, "directionality_results.csv", row.names = FALSE)

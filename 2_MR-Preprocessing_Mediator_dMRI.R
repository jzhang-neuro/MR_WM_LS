# Load required libraries
library(vroom)       # Efficient reading of large data files
library(TwoSampleMR) # Tools for two-sample Mendelian randomization

# Load custom function file
source("calculate_R2_F.R")  # Custom function to calculate R² and F-statistic

# Read the IDP trait data file
idp_traits <- read.csv("IDP_dMRI_trait.csv")

# Set filtering threshold and output directory for filtered weak instruments
Ffilter <- 10
output_Ffilter_path <- file.path(getwd(), "exp_dat_eaf_Ffilter_dMRI")

# Loop through each row in the IDP_dMRI_trait.csv file
for (i in 1:nrow(idp_traits)) {
  idp_name <- idp_traits[i, "IDP_dMRI_Trait"]
  tsv_url <- idp_traits[i, "ftp"]
  
  # Define the temporary file name for the downloaded data
  temp_file <- paste0(i, "_", gsub(" ", "_", idp_name), "_", idp_traits[i, "accessionId"], ".tsv")
  
  # Download and read the tsv file
  download.file(tsv_url, destfile = temp_file)
  IDP <- vroom(temp_file, col_names = TRUE)
  
  # Add samplesize to the IDP data
  current_samplesize <- idp_traits$samplesize[i]
  IDP$samplesize <- current_samplesize
  
  # Format the data for MR analysis
  exp_dat1 <- format_data(IDP, type = "exposure",
                          snp_col = "variant_id",
                          beta_col = "beta",
                          se_col = "standard_error",
                          effect_allele_col = "effect_allele",
                          other_allele_col = "other_allele",
                          eaf_col = "effect_allele_frequency",
                          pval_col = "p_value",
                          samplesize_col = "samplesize")
  
  # Filter the data based on p-value (threshold: 5e-8)
  exp_dat2 <- exp_dat1[exp_dat1$pval.exposure < 5e-8, ]
  
  # Check if exp_dat2 has fewer than 2 rows or is NULL
  if (is.null(exp_dat2) || nrow(exp_dat2) < 2) {
    cat("No significant SNPs found in exp_dat2 for clumping for IDP:", idp_name, "\n")
    next  # Skip to the next iteration
  }
  
  # Clump the data to ensure independence of SNPs
  exp_dat3 <- clump_data(exp_dat2)
  
  # Check if exp_dat3 has fewer than 2 rows or is NULL
  if (is.null(exp_dat3) || nrow(exp_dat3) < 2) {
    cat("Not enough SNPs in exp_dat3 for MR analysis for", idp_name, "\n")
    next
  }
  
  # Save the filtered data to a CSV file
  write.csv(exp_dat3, file = paste0("exp_dat_", gsub(" ", "_", idp_name), ".csv"))
  
  # Get the sample size for the current trait
  N <- idp_traits$samplesize[i]
  
  # Calculate R² and F-statistics using custom functions
  exp_dat3 <- calculate_R2_F(exp_dat3, N)
  
  # Filter out weak instruments based on the F-statistic threshold
  filtered_data <- exp_dat3[exp_dat3$F > Ffilter, ]
  
  # Construct the output file name and path, adding a "_Ffiltered" suffix
  output_csv_name <- paste0("exp_dat_Ffilter_", gsub(" ", "_", idp_name), ".csv")
  output_csv_path <- file.path(output_Ffilter_path, output_csv_name)
  
  # Save the filtered data to a new CSV file
  write.csv(filtered_data, output_csv_path, row.names = FALSE)
  
  # Extract outcome data (Example: Lacunar stroke)
  out_dat <- extract_outcome_data(
    snps = filtered_data$SNP,
    outcomes = 'ebi-a-GCST90014123'
  )
  
  # Check if out_dat is empty or has fewer than 2 rows
  if (is.null(out_dat) || nrow(out_dat) < 2) {
    cat("Not enough SNPs in out_dat for MR analysis for", idp_name, "\n")
    next
  }
  
  # Harmonize the exposure and outcome data
  dat <- harmonise_data(
    exposure_dat = filtered_data,
    outcome_dat = out_dat
  )
  
  # Perform Steiger directionality test
  steiger_result <- directionality_test(dat)
  
  # Save Steiger test results to a dataframe
  steiger_df <- data.frame(
    id.exposure = unique(dat$id.exposure),
    id.outcome = unique(dat$id.outcome),
    outcome = unique(dat$outcome),
    exposure = unique(dat$exposure),
    snp_r2.exposure = steiger_result$snp_r2.exposure,
    snp_r2.outcome = steiger_result$snp_r2.outcome,
    correct_direction = steiger_result$correct_causal_direction,
    steiger_pval = steiger_result$steiger_pval
  )
  
  # Check if the directionality is correct, if not, skip the analysis
  if (!steiger_result$correct) {
    cat("Directionality test failed for", idp_name, "\n")
    next
  }
  
  # Save the harmonized data to a CSV file
  write.csv(dat, file = paste0(i, "_dat_dMRI_on_lacunar_", gsub(" ", "_", idp_name), ".csv"))
}

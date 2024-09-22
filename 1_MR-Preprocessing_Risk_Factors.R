# MR Exposure Data Preprocessing for Three Categories of Risk factors
# 1. FinnGen datasets (Local Method to Read); 2. Published Supplementary Materials (Local Method to Read); 3. Online Method to Extract Preprocessed Data
# Required libraries
library(vroom)
library(dplyr)
library(TwoSampleMR)

# Preprocessing I. Read, Format, clumping, and Save
# 1. Local Method to Read FinnGen Data
# Initialize a list of FinnGen datasets and corresponding sample sizes
finngen_datasets <- list(
  list(file = "finngen_R10_I9_HYPTENS.gz", samplesize = 412113, name = "Hypertension"),
  list(file = "finngen_R10_BMI_IRN.gz", samplesize = 290820, name = "BMI"),
  list(file = "finngen_R10_E4_LIPOPROT.gz", samplesize = 395429, name = "Lipoprotein_Metabolism"),
  list(file = "finngen_R10_E4_HYPERLIPNAS.gz", samplesize = 365889, name = "Hyperlipidaemia"),
  list(file = "finngen_R10_G6_SLEEPAPNO.gz", samplesize = 410385, name = "SleepApnea")
)

# Read FinnGen datasets
for (data_info in finngen_datasets) {
  data <- vroom(data_info$file, col_names = TRUE)
  
  # Add the sample size column
  data <- data %>% mutate(samplesize = data_info$samplesize)
  
  # Format exposure data
  exp_dat1 <- format_data(data, type = "exposure", 
                          snp_col = "rsids", beta_col = "beta", se_col = "sebeta", 
                          effect_allele_col = "alt", other_allele_col = "ref", 
                          eaf_col = "af_alt", pval_col = "pval", samplesize_col = "samplesize")
  
  # Filter data with p-value < 5e-8 and apply clumping
  exp_dat2 <- exp_dat1[exp_dat1$pval.exposure < 5e-8,]
  exp_dat3 <- clump_data(exp_dat2)
  
  # Save the processed data
  write.csv(exp_dat3, file = paste0("exp_dat_", data_info$name, "_FinnGen.csv"), row.names = FALSE)
}

# 2. Local Method to Read Data from Published Supplementary Materials
# Hyperhomocysteinemia (2013, PubMedID 23824729), Insomnia (2019, PubMed ID 30804565), Type 2 diabetes (2018, PubMed ID 30297969)
pubmed_datasets <- list(
  list(file = "HHcy_SNP_table.csv", samplesize = 44147, name = "Hyperhomocysteinemia", 
       snp_col = "SNP", beta_col = "Beta", se_col = "SE", effect_allele_col = "Effect_allele", 
       other_allele_col = "other_allele", eaf_col = "EAF", pval_col = "p_meta"),
  
  list(file = "NG_S6_insomnia_lead_snp.csv", samplesize = 1331010, name = "Insomnia", 
       snp_col = "rsID", beta_col = "log_OR", se_col = "SE", effect_allele_col = "A1", 
       other_allele_col = "A2", eaf_col = "EAF", pval_col = "P", additional_mutation = "log(OR)"),
  
  list(file = "Diabetes_NatGenet_2018.csv", name = "Type2Diabetes", 
       custom_process = TRUE)  # This dataset requires custom processing
)

for (data_info in pubmed_datasets) {
  if (data_info$custom_process) {
    # Custom processing for Diabetes dataset
    diabetes_data <- read.csv(data_info$file)
    diabetes_data <- diabetes_data %>%
      mutate(
        rsids = Index.variant, beta = log(OR), 
        sebeta = (log(OR_UCI) - log(OR_LCI)) / (2 * 1.96), 
        alt = Risk.allele, ref = Other.allele, af_alt = RAF / 100, 
        pval = as.numeric(sub("x10-", "e-", p.value)), 
        samplesize = Cases + Controls
      ) %>%
      select(rsids, beta, sebeta, alt, ref, af_alt, pval, samplesize)
    write.csv(diabetes_data, file = "Diabetes_NatGenet_2018_preprocessed.csv")
    exp_dat1 <- format_data(diabetes_data, type = "exposure", snp_col = "rsids", 
                            beta_col = "beta", se_col = "sebeta", effect_allele_col = "alt", 
                            other_allele_col = "ref", eaf_col = "af_alt", pval_col = "pval", 
                            samplesize_col = "samplesize")
  } else {
    # General process for other datasets
    data <- read.csv(data_info$file)
    
    if (!is.null(data_info$additional_mutation)) {
      # Apply mutation for Insomnia dataset (log(OR))
      data <- data %>% mutate(log_OR = log(OR))
    }
    
    data <- data %>% mutate(samplesize = data_info$samplesize)
    
    exp_dat1 <- format_data(data, type = "exposure", 
                            snp_col = data_info$snp_col, beta_col = data_info$beta_col, 
                            se_col = data_info$se_col, effect_allele_col = data_info$effect_allele_col, 
                            other_allele_col = data_info$other_allele_col, eaf_col = data_info$eaf_col, 
                            pval_col = data_info$pval_col, samplesize_col = "samplesize")
  }
  
  exp_dat2 <- exp_dat1[exp_dat1$pval.exposure < 5e-8,]
  exp_dat3 <- clump_data(exp_dat2)
  
  # Save the processed data
  write.csv(exp_dat3, file = paste0("exp_dat_", data_info$name, "_PubMed.csv"), row.names = FALSE)
}

# 3. Online Method to Extract Preprocessed Data
online_datasets <- list(
  list(outcome = 'ieu-b-25', name = "Cigarettes_per_day"),
  list(outcome = 'ieu-b-73', name = "Alcohol_per_week")
)

for (data_info in online_datasets) {
  exp_dat3 <- extract_instruments(outcomes = data_info$outcome)
  
  # Save the processed data
  write.csv(exp_dat3, file = paste0("exp_dat_", data_info$name, "_online.csv"), row.names = FALSE)
}

## Preprocessing II. Initialize datasets for risk factors
# Load custom function file
source("calculate_R2_F.R")  # Custom function to calculate R² and F-statistic

# Set filtering threshold and output directory for filtered weak instruments
Ffilter <- 10
output_Ffilter_path <- file.path(getwd(), "exp_dat_eaf_Ffilter_Risk_factors")

# List all files that match the pattern 'exp_dat_*.csv'
files <- list.files(pattern = "exp_dat_.*\\.csv")

# Initialize a list to store the file paths, names, sources, and sample sizes for each risk factor
risk_factor_datasets <- list()

# Loop through the files and extract information
for (file in files) {
  
  # Split the filename into components
  # Assuming the format is exp_dat_{name}_{source}.csv
  file_info <- strsplit(file, "_")[[1]]
  
  # Extract the risk factor name and source (FinnGen, PubMed, Online)
  risk_factor_name <- file_info[3]  # e.g., Hypertension, BMI, etc.
  # risk_factor_source <- gsub("\\.csv", "", file_info[4])  # e.g., FinnGen, PubMed
  
  # Read the data from the CSV file
  exp_data <- read.csv(file)
  
  # Extract sample size if the column exists
  if ("samplesize" %in% colnames(exp_data)) {
    samplesize <- unique(exp_data$samplesize)
  } else {
    samplesize <- NA  # Set to NA if the sample size is unknown
  }
  
  # Add the risk factor information to the risk_factor_datasets list
  risk_factor_datasets[[length(risk_factor_datasets) + 1]] <- list(
    file = file,
    name = risk_factor_name,
    source = risk_factor_source,
    samplesize = samplesize
  )
}

# Each element in risk_factor_datasets now contains:
# - file: Path to the file
# - name: Risk factor name
# - source: Data source (FinnGen, PubMed, Online)
# - samplesize: Sample size (if available)

## Preprocessing III. Calculate R² and F-statistics, Filter out weak instruments, Extract outcome data, Harmonize, Steiger test,Save
# Iterate over all risk factor datasets, similar to the IDPs handling process
  # Load preprocessed data
  for (data_info in risk_factor_datasets) {
    # Load the preprocessed CSV file
    exp_dat3 <- read.csv(data_info$file)
    
    # Check if exp_dat3 has fewer than 2 rows or is NULL
    if (is.null(exp_dat3) || nrow(exp_dat3) < 2) {
      cat("Not enough SNPs in exp_dat3 for MR analysis for", data_info$name, "\n")
      next
    }
    
    # Save the preprocessed and filtered exposure data to a CSV file
    write.csv(exp_dat3, file = paste0("exp_dat_", gsub(" ", "_", data_info$name), ".csv"))
    
    # Get the sample size for the current trait
    N <- data_info$samplesize
    
    # Calculate R² and F-statistics using custom functions
    exp_dat3 <- calculate_R2_F(exp_dat3, N)
    
    # Filter out weak instruments based on the F-statistic threshold
    filtered_data <- exp_dat3[exp_dat3$F > Ffilter, ]
    
    # Construct the output file name and path, adding a "_Ffiltered" suffix
    output_csv_name <- paste0("exp_dat_Ffilter_", gsub(" ", "_", data_info$name), ".csv")
    output_csv_path <- file.path(output_Ffilter_path, output_csv_name)
    
    # Save the filtered data to a new CSV file
    write.csv(filtered_data, output_csv_path, row.names = FALSE)
    
    # Extract outcome data (Lacunar stroke)
    out_dat <- extract_outcome_data(
      snps = filtered_data$SNP,
      outcomes = 'ebi-a-GCST90014123'
    )
    
    # Check if out_dat is empty or has fewer than 2 rows
    if (is.null(out_dat) || nrow(out_dat) < 2) {
      cat("Not enough SNPs in out_dat for MR analysis for", data_info$name, "\n")
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
      cat("Directionality test failed for", data_info$name, "\n")
      next
    }
    
    # Save the harmonized data to a CSV file
    write.csv(dat, file = paste0(data_info$name, "_dat_riskfactor_on_lacunar_", gsub(" ", "_", data_info$name), ".csv"))
  }
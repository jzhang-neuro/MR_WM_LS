# Set path
setwd("/Custom_Path/...")

# Function for mediation analysis
mediation_analysis <- function(directA_file, directB_file, total_effect_file, output_file) {
  
  # Read Direct A: Exposure (Hypertension) -> Mediator
  res_mi <- subset(read.csv(directA_file), method == "Inverse variance weighted")
  
  # Read Direct B: Mediator -> Outcome (Lacunar Stroke)
  res_mi1 <- subset(read.csv(directB_file), method == "Inverse variance weighted")
  
  # Read total effect: Exposure (Hypertension) -> Outcome (Lacunar Stroke)
  res_mi2 <- subset(read.csv(total_effect_file), method == "Inverse variance weighted")
  
  # Calculate mediation effect
  midbeta = res_mi$b * res_mi1$b
  
  # Calculate mediation proportion
  mid_of_sum = midbeta / res_mi2$b
  
  # Calculate standard error for the mediation effect
  S = sqrt(res_mi$b^2 * res_mi1$se^2 + res_mi1$b^2 * res_mi$se^2)
  
  # Z-statistic for mediation effect
  Z = midbeta / S
  
  # P-value for mediation effect
  P = 2 * pnorm(q = abs(Z), lower.tail = FALSE)
  
  # Calculate 95% CI for total effect
  total_effect_beta = res_mi2$b
  total_effect_se = res_mi2$se
  total_effect_lower = total_effect_beta - 1.96 * total_effect_se
  total_effect_upper = total_effect_beta + 1.96 * total_effect_se
  
  # Calculate 95% CI for Direct A
  directA_beta = res_mi$b
  directA_se = res_mi$se
  directA_lower = directA_beta - 1.96 * directA_se
  directA_upper = directA_beta + 1.96 * directA_se
  
  # Calculate 95% CI for Direct B
  directB_beta = res_mi1$b
  directB_se = res_mi1$se
  directB_lower = directB_beta - 1.96 * directB_se
  directB_upper = directB_beta + 1.96 * directB_se
  
  # Calculate 95% CI for mediation effect
  mediation_effect_lower = midbeta - 1.96 * S
  mediation_effect_upper = midbeta + 1.96 * S
  
  # Standard error of the mediation proportion
  se_mid_of_sum = sqrt((1 / total_effect_beta)^2 * S^2 + (midbeta / total_effect_beta^2)^2 * total_effect_se^2)
  
  # Calculate 95% CI for mediation proportion
  mid_of_sum_lower = mid_of_sum - 1.96 * se_mid_of_sum
  mid_of_sum_upper = mid_of_sum + 1.96 * se_mid_of_sum
  
  # Create result dataframe
  results <- data.frame(
    Effect = c("Total effect", "Direct A", "Direct B", "Mediation effect", "Mediation Proportion"),
    Beta = c(total_effect_beta, directA_beta, directB_beta, midbeta, mid_of_sum),
    Lower_95_CI = c(total_effect_lower, directA_lower, directB_lower, mediation_effect_lower, mid_of_sum_lower),
    Upper_95_CI = c(total_effect_upper, directA_upper, directB_upper, mediation_effect_upper, mid_of_sum_upper),
    Mediation_Proportion = c(NA, NA, NA, mid_of_sum, mid_of_sum),
    Z_Statistic = c(NA, NA, NA, Z, NA),
    P_Value = c(NA, NA, NA, P, NA)
  )
  
  # Save results to CSV file
  write.csv(results, file = output_file, row.names = FALSE)
}

# Perform mediation analysis for the three pairs
# Group 1: L_ALIC_MD
mediation_analysis(
  directA_file = "res_Hypertension_on_MD_ALIC_L_pos_0001.csv",  # Direct A: Hypertension -> L_ALIC_MD
  directB_file = "res_MD_ALIC_L_on_Lacunar.csv",  # Direct B: L_ALIC_MD -> Lacunar
  total_effect_file = "res_Hypertension_on_lacunar.csv",  # Total effect: Hypertension -> Lacunar Stroke
  output_file = "res_mediation_Hypertension_LALIC-MD_on_Lacunar.csv"  # Output results
)

# Group 2: R_ALIC_MD
mediation_analysis(
  directA_file = "res_Hypertension_on_MD_ALIC_R_pos_0001.csv",  # Direct A: Hypertension -> R_ALIC_MD
  directB_file = "res_MD_ALIC_R_on_Lacunar.csv",  # Direct B: R_ALIC_MD -> Lacunar
  total_effect_file = "res_Hypertension_on_lacunar.csv",  # Total effect: Hypertension -> Lacunar Stroke
  output_file = "res_mediation_Hypertension_RALIC-MD_on_Lacunar.csv"  # Output results
)

# Group 3: R_ALIC_ISOVF
mediation_analysis(
  directA_file = "res_Hypertension_on_ISOVF_ALIC_R_pos_0001.csv",  # Direct A: Hypertension -> R_ALIC_ISOVF
  directB_file = "res_ISOVF_ALIC_R_on_Lacunar.csv",  # Direct B: R_ALIC_ISOVF -> Lacunar
  total_effect_file = "res_Hypertension_on_lacunar.csv",  # Total effect: Hypertension -> Lacunar Stroke
  output_file = "res_mediation_Hypertension_RALIC-ISOVF_on_Lacunar.csv"  # Output results
)

# Define a function to calculate R^2 and F statistics, with N as an input parameter
calculate_R2_F <- function(data, N) {
  
  # Calculate R^2 based on exposure beta, effect allele frequency (eaf), and sample size
  data <- transform(data, R2 = (2 * (beta.exposure^2) * eaf.exposure * (1 - eaf.exposure)) / 
                      (2 * (beta.exposure^2) * eaf.exposure * (1 - eaf.exposure) + 
                         2 * (se.exposure^2) * N * eaf.exposure * (1 - eaf.exposure)))
  
  # Calculate the F-statistic based on R^2 and sample size
  data <- transform(data, F = (R2 * (N - 2)) / (1 - R2))
  
  # Return the data with added R^2 and F columns
  return(data)
}

# Mendelian Randomization Analysis of Risk Factors and dMRI Traits

This repository contains R code for conducting Mendelian Randomization (MR) analyses of various risk factors and their relationships with dMRI traits, specifically focusing on neurological diseases. The analyses include preprocessing, MR analyses, and mediation analysis.

## Contents

1. **R Scripts**
   - **1_MR-Preprocessing_Risk_Factors.R**: Code for preprocessing risk factor data in the MR pipeline.
   - **2_MR-Preprocessing_Mediator_dMRI.R**: Code for preprocessing dMRI mediator data.
   - **3-MR-Analysis_Risk_Factors_on_Lacunar.R**: Code for MR analysis of risk factors on neurological diseases.
   - **4-MR-Analysis_dMRL_on_Lacunar.R**: Code for MR analysis of dMRI traits on neurological diseases.
   - **5-MR-Analysis_Risk_Factors_on_dMRI.R**: Code for MR analysis of risk factors on dMRI traits.
   - **6-MR-TwoStep_Mediation.R**: Code for conducting two-step mediation analyses involving risk factors and dMRI traits.

2. **Function File**
   - **calculate_R2_F.R**: A function for calculating R² and F statistics used in the analyses.

3. **Data File**
   - **IDP_dMRI_trait.csv**: A list of IDP dMRI traits used in the analyses.

## Getting Started

### Prerequisites
To run the scripts in this repository, ensure you have R installed on your computer along with the following R packages:
- `vroom`
- `dplyr`
- `ggplot2`
- `MendelianRandomization`

You can install the required packages using the following command in R:

```r
install.packages(c("vroom", "dplyr", "ggplot2", "MendelianRandomization"))

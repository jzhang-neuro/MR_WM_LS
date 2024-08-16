library(vroom)
library(TwoSampleMR)
library(ggplot2)
library(dplyr)
library(mr.raps)
library(tidyverse)

# function
source("/Users/neurotide/Downloads/R_MR_Gene/calculate_R2_F.R")
source("/Users/neurotide/Downloads/R_MR_Gene/calc_easy_R2_F.R")

# 创建一个空列表，用于保存结果
results_list <- list()

# 创建一个空的数据框，用于保存所有MR阳性结果
significant_results <- data.frame()

# 创建空的数据框，用于保存heterogeneity和pleiotropy测试结果
heterogeneity_results <- data.frame()
pleiotropy_results <- data.frame()

# 创建空的数据框，用于保存Steiger方向性检验结果
directionality_results <- data.frame()

# 正式读取IDP_dMRI_trait.csv文件
idp_traits <- read.csv("IDP_dMRI_trait.csv")

# 读取CSV文件，告诉read.csv没有列名
exp_dat_WM_549 <- read.csv("exp_dat_WM_549.csv", header = FALSE)

# 提取文件路径
file_paths <- exp_dat_WM_549[[1]]

# 设置过滤阈值和弱工具变量过滤后的输出文件夹路径
Ffilter <- 10

output_Ffilter_path <- file.path(getwd(), "exp_dat_eaf_Ffilter_WM_549")

# 遍历IDP_dMRI_trait.csv文件中的每一行
for (i in 1:nrow(idp_traits)) {
  idp_name <- idp_traits[i, "IDP_dMRI_Trait"]
  
  # 构造文件名
  csv_file <- paste0("exp_dat_", gsub(" ", "_", idp_name), ".csv")
  
  # 指定子文件夹路径
  subfolder_path <- "exp_dat_eaf_WM_549"
  
  # 构建CSV文件的完整路径
  csv_file_path <- file.path(subfolder_path, csv_file)
  
  # 检查文件是否存在于目录中
  if (csv_file %in% file_paths) {
    # 读取CSV文件
    exp_dat3 <- read.csv(csv_file_path, header = TRUE)
  } else {
    message(paste("No exp_dat for", idp_name))
    next  # 跳过这个循环
  }
  
  # 获取样本量
  N <- exp_dat3[1, "samplesize.exposure"]
  
  # 计算R^2和F值
  exp_dat3 <- calculate_R2_F(exp_dat3)
  
  # 计算简易版R^2和F值
  exp_dat3 <- calc_easy_R2_F(exp_dat3)
  
  # 过滤弱工具变量
  filtered_data <- exp_dat3[exp_dat3$F > Ffilter, ]
  
  # 构建输出文件路径，添加 "_Ffiltered" 后缀
  output_csv_name <- paste0("exp_dat_Ffilter_", gsub(" ", "_", idp_name), ".csv")
  output_csv_path <- file.path(output_Ffilter_path, output_csv_name)
  
  # 写入过滤后的数据到新的CSV文件
  write.csv(filtered_data, output_csv_path, row.names = FALSE)
  
  # 提取结果数据 Outcome
  out_dat <- extract_outcome_data(
    snps = filtered_data$SNP,
    outcomes = 'Outcome_ID')
  
  # 检查 out_dat 是否为空或行数少于2
  if (is.null(out_dat) || nrow(out_dat) < 2) {
    cat("Not enough SNPs in out_dat for MR analysis for", idp_name, "\n")
    next
  }
  # 协调数据(做MR分析时：暴露结局有EAF，用action=2;若无EAF，可用action=1作为主分析，action=3作为删除回文SNP的敏感性分析)
  dat <- harmonise_data(
    exposure_dat = filtered_data, 
    outcome_dat = out_dat,
    action = 1
  )
  
  # 进行Steiger方向性检验
  steiger_result <- directionality_test(dat)
  
  # 将Steiger检验结果保存到数据框
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
  
  # 检查方向性是否正确，如果不正确则跳过该分析
  if (!steiger_result$correct) {
    cat("Directionality test failed for", idp_name, "\n")
    next
  }
  
  # 保存协调后的数据
  write.csv(dat, file = paste0(i, "_dat_lacunar_", gsub(" ", "_", idp_name), ".csv"))
  
  # 筛选有效的工具变量
  valid_dat <- dat[dat$mr_keep == TRUE, ]
 
   # 判断工具变量个数是否足够
  if (is.null(valid_dat) || nrow(valid_dat) < 2) {
    cat("Not enough SNPs in valid instruments for MR analysis for", idp_name, "\n")
    next
  }
 
  # 执行MR分析
  res <- mr(valid_dat, method_list = c("mr_ivw", "mr_ivw_fe", "mr_egger_regression", "mr_weighted_mode", "mr_weighted_median"))
  
  # 提取IVW方法的 p 值
  ivw_pval <- res[res$method == "Inverse variance weighted", "pval"]
  
  # 计算mr.raps
  mr_raps_result <- mr.raps(valid_dat$beta.exposure, valid_dat$beta.outcome, valid_dat$se.exposure, valid_dat$se.outcome)
  
  # 整理mr.raps结果
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
  
  # 合并res和mr_raps的结果
  res <- rbind(res, mr_raps_df)
  
  # 执行MR-PRESSO异常值检测偏倚
  # 判断工具变量个数是否足够
  if (nrow(valid_dat) >= 4) {
    # 尝试运行 MR-PRESSO 分析
    tryCatch({
      presso_result <- run_mr_presso(dat = valid_dat)
      
      # 创建保存结果的数据框
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
      
      # 将所有部分的结果合并到一个列表中
      result_list <- list(
        MR_PRESSO_results = mr_presso_df,
        MR_PRESSO_global_test = global_test_df,
        MR_PRESSO_outlier_test = outlier_test_df,
        MR_PRESSO_distortion_test = distortion_test_df
      )
      
      # 将列表中的每个数据框写入到CSV文件中
      file_prefix <- gsub(" ", "_", idp_name)
      for (name in names(result_list)) {
        write.csv(result_list[[name]], file = paste0(i, "_", name, "_", file_prefix, ".csv"), row.names = FALSE)
      }
      
      # 提取Outlier-corrected的结果
      main_mr_results <- presso_result[[1]]$`Main MR results`
      outlier_corrected <- main_mr_results[main_mr_results$`MR Analysis` == 'Outlier-corrected', ]
      
      causal_estimate <- outlier_corrected$`Causal Estimate`
      sd <- outlier_corrected$Sd
      p_value <- outlier_corrected$`P-value`
      
      # 提取Outlier Test中的Pvalue并转换为数值
      outlier_test <- presso_result[[1]]$`MR-PRESSO results`$`Outlier Test`
      outlier_test$Pvalue <- as.numeric(gsub("[^0-9.]", "", outlier_test$Pvalue))
      
      # 计算Pvalue < 0.05的行数
      nsnp_outlier <- sum(outlier_test$Pvalue < 0.05, na.rm = TRUE)
      
      # 计算有效的nsnp
      nsnp_effective <- nrow(valid_dat) - nsnp_outlier
      
      # 创建MR-PRESSO Outlier-corrected的结果数据框
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
      
      # 将MR-PRESSO Outlier-corrected的结果整合到res数据框中
      res <- rbind(res, mr_presso_outlier_corrected)
      
    }, error = function(e) {
      message("Error in MR-PRESSO analysis: ", e$message)
    })
  } else {
    message("Not enough instrumental variables for MR-PRESSO analysis. Skipping this step.")
  }
  
  # 生成odds ratios
  or_df <- generate_odds_ratios(res)
  
  # 合并结果和 odds ratios 数据框
  result <- merge(res, or_df, by = c("id.exposure", "id.outcome", "outcome", "exposure", "method", "nsnp", "b", "se", "pval"))
  
  # 计算总体的R^2
  overall_R2 <- sum(valid_dat$R2)
  overall_R2.easy <- sum(valid_dat$R2.easy)
  
  # 工具变量的数量是有效工具变量的数量
  k <- nrow(valid_dat)
  
  # 计算总体的F值
  overall_F_statistic <- (overall_R2 * (N - k - 1)) / (k * (1 - overall_R2))
  overall_F.easy_statistic <- (overall_R2.easy * (N - k - 1)) / (k * (1 - overall_R2.easy))
  
  # 添加总体的R^2和F统计量到结果
  result$overall_R2 <- overall_R2
  result$overall_F_statistic <- overall_F_statistic
  result$overall_R2.easy <- overall_R2.easy
  result$overall_F.easy_statistic <- overall_F.easy_statistic
  
  # 保存结果
  # p_value <- min(result$pval)  # 获取最小 p 值
  suffix <- ifelse(ivw_pval <= 0.0001, "_pos_0001",
                   ifelse(ivw_pval <= 0.00025, "_pos_00025",
                          ifelse(ivw_pval <= 0.001, "_pos_001",
                                 ifelse(ivw_pval <= 0.005, "_pos_005",
                                        ifelse(ivw_pval <= 0.05, "_pos_05", "")))))  # 根据 p 值范围设置后缀
  
  write.csv(result, file = paste0(i, "_res_", gsub(" ", "_", idp_name), suffix, ".csv"), row.names = FALSE)
  
  # 检查是否存在 "mr_ivw" 这个方法的 p 值小于 0.05
  if (any(res$method == "Inverse variance weighted" & res$pval < 0.05)) {
    result$IDP <- idp_name  # 添加IDP名称列
    result$i <- i  # 添加循环编号列
    significant_results <- rbind(significant_results, result)
  }
  
  # 可视化
  p1 <- mr_scatter_plot(res, valid_dat)
  ggsave(p1[[1]], file = paste0(i, "_", gsub(" ", "_", idp_name), "_scatter.pdf"), width = 7, height = 7)
  
  res_single <- mr_singlesnp(valid_dat)
  p2 <- mr_forest_plot(res_single)
  ggsave(p2[[1]], file = paste0(i, "_", gsub(" ", "_", idp_name), "_forest.pdf"), width = 7, height = 7)
  
  res_loo <- mr_leaveoneout(valid_dat)
  p3 <- mr_leaveoneout_plot(res_loo)
  ggsave(p3[[1]], file = paste0(i, "_", gsub(" ", "_", idp_name), "_loo.pdf"), width = 7, height = 7)
  
  p4 <- mr_funnel_plot(res_single)
  ggsave(p4[[1]], file = paste0(i, "_", gsub(" ", "_", idp_name), "_funnel.pdf"), width = 7, height = 7)
  
  # 额外的分析步骤：
  # 存储heterogeneity结果并添加IDP名称列
  het <- mr_heterogeneity(valid_dat)
  het$ivw_pval <- ivw_pval
  het$IDP <- idp_name
  het$i <- i  # 添加循环编号列
  heterogeneity_results <- rbind(heterogeneity_results, het)
  
  # 存储pleiotropy测试结果并添加IDP名称列
  pleio <- mr_pleiotropy_test(valid_dat)
  pleio$ivw_pval <- ivw_pval
  pleio$IDP <- idp_name
  pleio$i <- i  # 添加循环编号列
  pleiotropy_results <- rbind(pleiotropy_results, pleio)
  
  # 将循环编号列、IDP、ivw_pval加到steiger_df
  steiger_df$ivw_pval <- ivw_pval
  steiger_df$IDP <- idp_name
  steiger_df$i <- i
  # 将本循环的Steiger结果追加到directionality_results数据框
  directionality_results <- rbind(directionality_results, steiger_df)
}

# 保存heterogeneity, pleiotropy, Steiger方向性检验结果为CSV文件
write.csv(heterogeneity_results, file = "heterogeneity_results.csv", row.names = FALSE)
write.csv(pleiotropy_results, file = "pleiotropy_results.csv", row.names = FALSE)
write.csv(directionality_results, file = "directionality_results.csv", row.names = FALSE)
# 保存所有阳性结果到一个文件
write.csv(significant_results, file = "significant_results.csv", row.names = FALSE)

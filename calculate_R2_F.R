# 定义计算R^2和F值的函数
calculate_R2_F <- function(data) {
  # 计算R^2
  data <- transform(data, R2 = (2 * (beta.exposure^2) * eaf.exposure * (1 - eaf.exposure)) / 
                      (2 * (beta.exposure^2) * eaf.exposure * (1 - eaf.exposure) + 
                         2 * (se.exposure^2) * samplesize.exposure * eaf.exposure * (1 - eaf.exposure)))
  
  # 计算F值
  N <- data$samplesize.exposure[1]
  data <- transform(data, F = (R2 * (N - 2)) / (1 - R2))
  
  return(data)
}
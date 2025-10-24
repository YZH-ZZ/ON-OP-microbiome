# 安装和加载必要的包
if (!require("MMUPHin")) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("MMUPHin")
}

if (!require("tidyverse")) install.packages("tidyverse")
if (!require("vegan")) install.packages("vegan")
if (!require("ape")) install.packages("ape")
if (!require("patchwork")) install.packages("patchwork")
if (!require("ggrepel")) install.packages("ggrepel")

library(MMUPHin)
library(tidyverse)
library(vegan)
library(ape)
library(patchwork)
library(ggrepel)

# 设置随机种子
set.seed(123)

# 1. 读取数据 ----------------------------------------------------
cat("Reading microbial composition data...\n")
microbial_data <- read.table("123.txt", header = TRUE, sep = "\t", 
                             stringsAsFactors = FALSE, check.names = FALSE)
cat("Microbial data dimensions:", dim(microbial_data), "\n")

cat("Reading batch information data...\n")
batch_info <- read.table("234.txt", header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)
cat("Batch info dimensions:", dim(batch_info), "\n")

# 2. 数据预处理 --------------------------------------------------
# 提取样本名、队列信息和微生物特征
sample_names <- microbial_data[, 1]
cohort_info <- microbial_data[, 2]
microbial_features <- microbial_data[, -(1:2)]

# 转换为矩阵
feature_matrix <- as.matrix(microbial_features)
rownames(feature_matrix) <- sample_names

# 检查并处理零方差特征
zero_variance_features <- apply(feature_matrix, 2, var) == 0
cat("Removing", sum(zero_variance_features), "zero-variance features\n")
feature_matrix <- feature_matrix[, !zero_variance_features]

# 添加伪计数并转换为相对丰度
pseudocount <- min(feature_matrix[feature_matrix > 0]) / 2
feature_matrix_pseudo <- feature_matrix + pseudocount
feature_matrix_rel <- feature_matrix_pseudo / rowSums(feature_matrix_pseudo)

# 准备样本元数据
sample_metadata <- data.frame(
  sample_id = sample_names,
  cohort = cohort_info
) %>% 
  left_join(batch_info, by = c("cohort" = "Cohort")) %>%
  mutate(sample_id = as.character(sample_id))

# 设置元数据行名
rownames(sample_metadata) <- sample_metadata$sample_id

# 3. 使用MMUPHin进行批次校正 ------------------------------------
cat("\nPerforming batch correction with MMUPHin...\n")

batch_variables <- c("Bacth1", "Bacth2", "Bacth3")
mmuphin_results <- list()

for (batch_var in batch_variables) {
  cat("\nCorrecting for batch variable:", batch_var, "\n")
  
  current_batch <- sample_metadata[[batch_var]]
  batch_counts <- table(current_batch)
  cat("Batch distribution:\n")
  print(batch_counts)
  
  if (length(unique(current_batch)) < 2) {
    cat("Skipping", batch_var, "- insufficient batch variation\n")
    next
  }
  
  # 准备MMUPHin输入数据
  feature_abd_mmuphin <- t(feature_matrix_rel)
  common_samples <- intersect(colnames(feature_abd_mmuphin), rownames(sample_metadata))
  
  if (length(common_samples) == 0) {
    cat("ERROR: No common samples found!\n")
    next
  }
  
  feature_abd_subset <- feature_abd_mmuphin[, common_samples, drop = FALSE]
  metadata_subset <- sample_metadata[common_samples, , drop = FALSE]
  
  tryCatch({
    fit <- adjust_batch(
      feature_abd = feature_abd_subset,
      batch = batch_var,
      data = metadata_subset,
      control = list(verbose = FALSE)
    )
    
    corrected_abd <- t(fit$feature_abd_adj)
    mmuphin_results[[batch_var]] <- corrected_abd
    
    cat("Successfully corrected for", batch_var, "\n")
    
  }, error = function(e) {
    cat("Error in MMUPHin correction:", e$message, "\n")
  })
}

# 4. 评估批次校正效果 -------------------------------------------
if (length(mmuphin_results) > 0) {
  successful_correction <- names(mmuphin_results)[1]
  final_corrected <- mmuphin_results[[successful_correction]]
  
  cat("\nUsing correction from:", successful_correction, "\n")
  
  # 保存校正后的数据
  corrected_df <- data.frame(
    Sample_ID = rownames(final_corrected),
    Cohort = sample_metadata$cohort[match(rownames(final_corrected), sample_metadata$sample_id)],
    final_corrected,
    check.names = FALSE
  )
  
  write.table(corrected_df, "microbial_data_MMUPHin_corrected.txt", 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # 5. PCoA分析和可视化（仅显示队列，不显示batch）---------------
  pcoa_analysis <- function(data_mat, title, batch_var) {
    # 计算Bray-Curtis距离
    bray_dist <- vegdist(data_mat, method = "bray")
    
    # PCoA分析
    pcoa_result <- pcoa(bray_dist)
    variance_explained <- pcoa_result$values$Relative_eig * 100
    
    # PERMANOVA检验 - 批次效应
    batch_groups <- sample_metadata[[batch_var]][match(rownames(data_mat), sample_metadata$sample_id)]
    permanova_batch <- adonis2(bray_dist ~ batch_groups)
    batch_r2 <- permanova_batch$R2[1]
    batch_pvalue <- permanova_batch$`Pr(>F)`[1]
    
    # PERMANOVA检验 - 队列效应（生物学差异）
    cohort_groups <- sample_metadata$cohort[match(rownames(data_mat), sample_metadata$sample_id)]
    permanova_cohort <- adonis2(bray_dist ~ cohort_groups)
    cohort_r2 <- permanova_cohort$R2[1]
    cohort_pvalue <- permanova_cohort$`Pr(>F)`[1]
    
    # 创建绘图数据（仅使用队列信息）
    pcoa_df <- data.frame(
      PCo1 = pcoa_result$vectors[, 1],
      PCo2 = pcoa_result$vectors[, 2],
      Cohort = cohort_groups,  # 只使用队列信息
      Sample = rownames(data_mat)
    )
    
    # 过滤坐标轴范围
    
    # 计算每个队列的中心点（用于统计信息显示）
    cohort_centers <- pcoa_df %>%
      group_by(Cohort) %>%
      summarise(
        Center_PCo1 = mean(PCo1),
        Center_PCo2 = mean(PCo2)
      )
    
    # 创建PCoA图（仅用颜色区分队列）
    p <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2, color = Cohort)) +
      geom_point(alpha = 0.7, size = 2) +
      stat_ellipse(level = 0.68, alpha = 0.3, linewidth = 0.8)  +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.5) +
      labs(
        title = title,
        subtitle = paste0(
          "Cohort effect: R² = ", round(cohort_r2, 3), ", p = ", 
          ifelse(cohort_pvalue < 0.001, "< 0.001", format.pval(cohort_pvalue, digits = 3))
        ),
        x = paste0("PCo1 (", round(variance_explained[1], 1), "%)"),
        y = paste0("PCo2 (", round(variance_explained[2], 1), "%)"),
        color = "Cohort"
      ) +
      theme_test() +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, face = "bold"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()
      ) +
      scale_color_brewer(palette = "Set1") +
      guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
    
    return(list(
      plot = p,
      permanova_batch = list(r2 = batch_r2, pvalue = batch_pvalue),
      permanova_cohort = list(r2 = cohort_r2, pvalue = cohort_pvalue),
      filtered_data = pcoa_df,
      cohort_centers = cohort_centers
    ))
  }
  
  # 6. 生成校正前后PCoA图 --------------------------------------
  cat("\nGenerating PCoA plots...\n")
  
  before_results <- pcoa_analysis(
    feature_matrix_rel, 
    "Before Batch Correction",
    successful_correction
  )
  
  after_results <- pcoa_analysis(
    final_corrected,
    "After Batch Correction", 
    successful_correction
  )
  
  # 7. 创建统计指标对比图 --------------------------------------
  # 批次效应指标
  batch_metrics <- data.frame(
    Condition = c("Before", "After"),
    Batch_R2 = c(before_results$permanova_batch$r2, after_results$permanova_batch$r2),
    Batch_Pvalue = c(before_results$permanova_batch$pvalue, after_results$permanova_batch$pvalue)
  )
  
  batch_plot <- ggplot(batch_metrics, aes(x = Condition, y = Batch_R2, fill = Condition)) +
    geom_col(alpha = 0.7, width = 0.6) +
    geom_text(aes(label = paste0("R² = ", round(Batch_R2, 3), 
                                 "\np = ", 
                                 ifelse(Batch_Pvalue < 0.001, "< 0.001", 
                                        format.pval(Batch_Pvalue, digits = 3)))), 
              vjust = -0.3, size = 3.5, lineheight = 0.8) +
    labs(title = "Batch Effect Reduction",
         y = "PERMANOVA R² (Batch Effect)",
         x = "") +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    ) +
    scale_fill_manual(values = c("Before" = "#E41A1C", "After" = "#377EB8")) +
    ylim(0, max(batch_metrics$Batch_R2) * 1.3)
  
  # 队列效应指标（生物学信号）
  cohort_metrics <- data.frame(
    Condition = c("Before", "After"),
    Cohort_R2 = c(before_results$permanova_cohort$r2, after_results$permanova_cohort$r2),
    Cohort_Pvalue = c(before_results$permanova_cohort$pvalue, after_results$permanova_cohort$pvalue)
  )
  
  cohort_plot <- ggplot(cohort_metrics, aes(x = Condition, y = Cohort_R2, fill = Condition)) +
    geom_col(alpha = 0.7, width = 0.6) +
    geom_text(aes(label = paste0("R² = ", round(Cohort_R2, 3), 
                                 "\np = ", 
                                 ifelse(Cohort_Pvalue < 0.001, "< 0.001", 
                                        format.pval(Cohort_Pvalue, digits = 3)))), 
              vjust = -0.3, size = 3.5, lineheight = 0.8) +
    labs(title = "Cohort Effect Preservation",
         y = "PERMANOVA R² (Cohort Effect)", 
         x = "") +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11)
    ) +
    scale_fill_manual(values = c("Before" = "#4DAF4A", "After" = "#984EA3")) +
    ylim(0, max(cohort_metrics$Cohort_R2) * 1.3)
  
  # 8. 组合所有图形 --------------------------------------------
  # 第一行：PCoA图
  pcoa_row <- before_results$plot + after_results$plot
  
  # 第二行：统计指标
  metrics_row <- batch_plot + cohort_plot
  
  # 组合所有图形
  final_plot <- pcoa_row / metrics_row + 
    plot_layout(heights = c(2, 1))
  
  ggsave("batch_correction_cohort_focus.png", final_plot, 
         width = 16, height = 12, dpi = 300)
  
  # 9. 输出详细统计结果 ----------------------------------------
  cat("\n=== BATCH CORRECTION EVALUATION ===\n")
  
  cat("\nBATCH EFFECT (should decrease after correction):\n")
  cat("Before - R²:", round(before_results$permanova_batch$r2, 4), 
      "p =", ifelse(before_results$permanova_batch$pvalue < 0.001, "< 0.001", 
                    format.pval(before_results$permanova_batch$pvalue, digits = 4)), "\n")
  cat("After  - R²:", round(after_results$permanova_batch$r2, 4),
      "p =", ifelse(after_results$permanova_batch$pvalue < 0.001, "< 0.001", 
                    format.pval(after_results$permanova_batch$pvalue, digits = 4)), "\n")
  
  batch_r2_reduction <- before_results$permanova_batch$r2 - after_results$permanova_batch$r2
  batch_improvement_pct <- ifelse(before_results$permanova_batch$r2 > 0, 
                                  round(batch_r2_reduction/before_results$permanova_batch$r2 * 100, 1), 0)
  
  cat("Batch effect reduction: R² decreased by", round(batch_r2_reduction, 4),
      paste0("(", batch_improvement_pct, "% improvement)\n"))
  
  cat("\nCOHORT EFFECT (should be preserved after correction):\n")
  cat("Before - R²:", round(before_results$permanova_cohort$r2, 4),
      "p =", ifelse(before_results$permanova_cohort$pvalue < 0.001, "< 0.001", 
                    format.pval(before_results$permanova_cohort$pvalue, digits = 4)), "\n")
  cat("After  - R²:", round(after_results$permanova_cohort$r2, 4),
      "p =", ifelse(after_results$permanova_cohort$pvalue < 0.001, "< 0.001", 
                    format.pval(after_results$permanova_cohort$pvalue, digits = 4)), "\n")
  
  cohort_r2_change <- after_results$permanova_cohort$r2 - before_results$permanova_cohort$r2
  cat("Cohort effect change: R²", ifelse(cohort_r2_change >= 0, "increased", "decreased"), 
      "by", round(abs(cohort_r2_change), 4), "\n")
  
  # 10. 评估校正效果 -------------------------------------------
  cat("\n=== CORRECTION EFFECTIVENESS ASSESSMENT ===\n")
  
  if (batch_improvement_pct >= 70 && abs(cohort_r2_change) <= 0.05) {
    cat("✓ EXCELLENT: Strong batch reduction with preserved biology\n")
  } else if (batch_improvement_pct >= 50 && abs(cohort_r2_change) <= 0.08) {
    cat("✓ GOOD: Effective batch correction\n") 
  } else if (batch_improvement_pct >= 30) {
    cat("○ MODERATE: Some batch reduction achieved\n")
  } else {
    cat("✗ MINIMAL: Limited batch correction effect\n")
  }
  
  # 11. 保存统计结果表格 ---------------------------------------
  results_summary <- data.frame(
    Metric = c("Batch_R2_Before", "Batch_R2_After", "Batch_Pvalue_Before", "Batch_Pvalue_After",
               "Cohort_R2_Before", "Cohort_R2_After", "Cohort_Pvalue_Before", "Cohort_Pvalue_After",
               "Batch_R2_Reduction", "Batch_Improvement_Percent", "Cohort_R2_Change"),
    Value = c(before_results$permanova_batch$r2, after_results$permanova_batch$r2,
              before_results$permanova_batch$pvalue, after_results$permanova_batch$pvalue,
              before_results$permanova_cohort$r2, after_results$permanova_cohort$r2,
              before_results$permanova_cohort$pvalue, after_results$permanova_cohort$pvalue,
              batch_r2_reduction, batch_improvement_pct, cohort_r2_change)
  )
  
  write.table(results_summary, "batch_correction_statistical_summary.txt", 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("\nResults saved:\n")
  cat("- microbial_data_MMUPHin_corrected.txt\n")
  cat("- batch_correction_cohort_focus.png\n")
  cat("- batch_correction_statistical_summary.txt\n")
  
} else {
  cat("No successful batch correction was performed.\n")
}

cat("\n=== PROCESSING COMPLETED ===\n")
# 加载必要的包
library(randomForest)
library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 读取数据
data <- read.table("123.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 检查数据结构
cat("数据维度:", dim(data), "\n")
cat("分组情况:", unique(data$group), "\n")

# 数据预处理
data$group <- as.factor(data$group)

# 只保留指定的四个队列
target_cohorts <- c("PRJNA1039567", "PRJNA565497", "PRJNA631117", "PRJNA796857")
data <- data[data$Cohort %in% target_cohorts, ]
data$Cohort <- factor(data$Cohort, levels = target_cohorts)

cat("筛选后的队列:", unique(data$Cohort), "\n")
cat("筛选后数据维度:", dim(data), "\n")

# 函数：执行二分类留一队列交叉验证
binary_loocv_analysis <- function(data, class1, class2, analysis_name) {
  # 筛选指定类别的数据
  binary_data <- data[data$group %in% c(class1, class2), ]
  binary_data$group <- droplevels(binary_data$group)
  
  cat("\n=== 正在进行", analysis_name, "分析 ===\n")
  cat("类别:", class1, "vs", class2, "\n")
  cat("总样本数:", nrow(binary_data), "\n")
  
  # 获取包含指定类别的队列，且每个队列必须同时包含两个类别
  cohort_summary <- binary_data %>%
    group_by(Cohort) %>%
    summarize(
      n_samples = n(),
      n_class1 = sum(group == class1),
      n_class2 = sum(group == class2),
      has_both_classes = n_distinct(group) == 2
    )
  
  # 只选择同时包含两个类别的队列
  valid_cohorts <- cohort_summary$Cohort[cohort_summary$has_both_classes]
  
  cat("有效队列（同时包含", class1, "和", class2, "）:", paste(valid_cohorts, collapse = ", "), "\n")
  
  if (length(valid_cohorts) == 0) {
    cat("警告: 没有队列同时包含", class1, "和", class2, "，跳过分析\n")
    return(NULL)
  }
  
  # 存储结果
  roc_results <- list()
  importance_results <- list()
  
  # 对每个有效队列进行留一验证
  for (i in 1:length(valid_cohorts)) {
    validation_cohort <- valid_cohorts[i]
    discovery_cohorts <- valid_cohorts[-i]
    
    cat("\n第", i, "次验证 - 验证队列:", validation_cohort, "\n")
    cat("发现队列:", paste(discovery_cohorts, collapse = ", "), "\n")
    
    # 分割数据
    discovery_data <- binary_data[binary_data$Cohort %in% discovery_cohorts, ]
    validation_data <- binary_data[binary_data$Cohort == validation_cohort, ]
    
    # 检查验证集是否包含两个类别
    validation_counts <- table(validation_data$group)
    cat("验证集类别分布:", validation_counts, "\n")
    
    if (length(validation_counts) < 2) {
      cat("警告: 验证队列", validation_cohort, "只包含一个类别，跳过\n")
      next
    }
    
    # 检查发现集样本数
    discovery_counts <- table(discovery_data$group)
    cat("发现集类别分布:", discovery_counts, "\n")
    
    if (nrow(discovery_data) < 10) {
      cat("警告: 发现集样本数不足，跳过\n")
      next
    }
    
    # 移除样本ID和队列列
    discovery_features <- discovery_data[, -c(1, 2)]
    validation_features <- validation_data[, -c(1, 2)]
    
    # 设置平衡抽样大小
    min_class_size <- min(discovery_counts)
    samp_size <- rep(min_class_size, length(discovery_counts))
    
    # 训练随机森林模型
    rf_model <- tryCatch({
      randomForest(
        group ~ .,
        data = discovery_features,
        ntree = 500,
        importance = TRUE,
        strata = discovery_features$group,
        sampsize = samp_size
      )
    }, error = function(e) {
      cat("模型训练错误:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(rf_model)) {
      next
    }
    
    # 在验证集上进行预测
    validation_predictions <- predict(rf_model, validation_features, type = "prob")
    
    # 计算AUROC（以class2为正类）
    roc_obj <- tryCatch({
      roc(validation_features$group, validation_predictions[, as.character(class2)])
    }, error = function(e) {
      cat("ROC计算错误:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(roc_obj)) {
      next
    }
    
    roc_ci <- ci.auc(roc_obj)
    
    # 存储ROC结果
    roc_results[[i]] <- list(
      cohort = validation_cohort,
      roc_object = roc_obj,
      auc = roc_obj$auc,
      ci_lower = roc_ci[1],
      ci_upper = roc_ci[3],
      discovery_counts = discovery_counts,
      validation_counts = validation_counts,
      analysis_name = analysis_name
    )
    
    # 获取特征重要性
    imp <- importance(rf_model)
    importance_results[[i]] <- data.frame(
      cohort = validation_cohort,
      analysis = analysis_name,
      feature = rownames(imp),
      importance = imp[, "MeanDecreaseAccuracy"],
      stringsAsFactors = FALSE
    )
    
    cat("AUROC:", round(roc_obj$auc, 3), "\n")
  }
  
  return(list(roc_results = roc_results, importance_results = importance_results))
}

# 函数：创建ROC数据框
create_roc_dataframe <- function(roc_results) {
  roc_data <- data.frame()
  
  for (i in 1:length(roc_results)) {
    if (!is.null(roc_results[[i]])) {
      result <- roc_results[[i]]
      
      # 安全地提取ROC曲线坐标
      tryCatch({
        # 获取ROC曲线的坐标点
        roc_coords_list <- coords(result$roc_object, transpose = FALSE)
        
        if (!is.null(roc_coords_list) && length(roc_coords_list$sensitivity) > 0) {
          # 创建坐标数据框
          coords_df <- data.frame(
            Cohort = result$cohort,
            Analysis = result$analysis_name,
            Sensitivity = roc_coords_list$sensitivity,
            Specificity = roc_coords_list$specificity,
            FPR = 1 - roc_coords_list$specificity,
            stringsAsFactors = FALSE
          )
          
          # 添加AUC和CI信息（每行都添加相同的值）
          coords_df$AUC <- result$auc
          coords_df$CI_lower <- result$ci_lower
          coords_df$CI_upper <- result$ci_upper
          
          roc_data <- rbind(roc_data, coords_df)
        }
      }, error = function(e) {
        cat("提取ROC坐标时出错:", e$message, "\n")
      })
    }
  }
  
  return(roc_data)
}

# 函数：创建综合ROC图
create_combined_roc_plot <- function(nc_on_results, nc_op_results) {
  # 收集所有ROC数据
  all_roc_data <- data.frame()
  
  if (!is.null(nc_on_results) && length(nc_on_results$roc_results) > 0) {
    nc_on_data <- create_roc_dataframe(nc_on_results$roc_results)
    all_roc_data <- rbind(all_roc_data, nc_on_data)
  }
  
  if (!is.null(nc_op_results) && length(nc_op_results$roc_results) > 0) {
    nc_op_data <- create_roc_dataframe(nc_op_results$roc_results)
    all_roc_data <- rbind(all_roc_data, nc_op_data)
  }
  
  if (nrow(all_roc_data) == 0) {
    cat("没有可用的ROC数据\n")
    return(NULL)
  }
  
  # 创建标注数据（每个队列和分析的唯一AUC值）
  annotation_data <- all_roc_data %>%
    group_by(Cohort, Analysis) %>%
    summarize(
      AUC = first(AUC),
      CI_lower = first(CI_lower),
      CI_upper = first(CI_upper),
      .groups = "drop"
    ) %>%
    mutate(
      Label = sprintf("%s\nAUC: %.3f (%.3f-%.3f)", 
                      Cohort, AUC, CI_lower, CI_upper),
      x_pos = 0.65,
      y_pos = 0.35
    )
  
  # 确保所有四个队列都在图中显示，即使某些队列在某些分析中没有结果
  all_cohorts_analysis <- expand.grid(
    Cohort = target_cohorts,
    Analysis = unique(all_roc_data$Analysis)
  )
  
  annotation_data <- annotation_data %>%
    right_join(all_cohorts_analysis, by = c("Cohort", "Analysis")) %>%
    mutate(
      Label = ifelse(is.na(AUC), 
                     sprintf("%s\nNo data", Cohort),
                     sprintf("%s\nAUC: %.3f (%.3f-%.3f)", Cohort, AUC, CI_lower, CI_upper)),
      x_pos = 0.65,
      y_pos = 0.35
    )
  
  # 创建综合ROC图
  p_roc <- ggplot(all_roc_data, aes(x = FPR, y = Sensitivity, color = Cohort, group = Cohort)) +
    geom_line(size = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray", size = 0.5) +
    geom_label(data = annotation_data, 
               aes(x = x_pos, y = y_pos, label = Label, color = Cohort),
               size = 3, show.legend = FALSE, alpha = 0.8,
               inherit.aes = FALSE) +
    facet_wrap(~ Analysis, ncol = 2) +
    labs(
      title = "ROC曲线分析",
      x = "假阳性率 (1 - 特异性)",
      y = "真阳性率 (灵敏度)",
      color = "验证队列"
    ) +
    theme_test() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    ) +
    scale_color_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"))  # 固定颜色
  
  return(p_roc)
}

# 函数：创建综合特征重要性图
create_combined_importance_plot <- function(nc_on_results, nc_op_results, top_n = 5) {
  all_imp_data <- data.frame()
  
  if (!is.null(nc_on_results) && length(nc_on_results$importance_results) > 0) {
    nc_on_imp <- do.call(rbind, nc_on_results$importance_results)
    all_imp_data <- rbind(all_imp_data, nc_on_imp)
  }
  
  if (!is.null(nc_op_results) && length(nc_op_results$importance_results) > 0) {
    nc_op_imp <- do.call(rbind, nc_op_results$importance_results)
    all_imp_data <- rbind(all_imp_data, nc_op_imp)
  }
  
  if (nrow(all_imp_data) == 0) {
    cat("没有可用的重要性数据\n")
    return(NULL)
  }
  
  # 计算平均重要性，只取前5个
  avg_importance <- all_imp_data %>%
    group_by(feature, analysis) %>%
    summarize(mean_importance = mean(importance), .groups = "drop") %>%
    group_by(analysis) %>%
    arrange(desc(mean_importance)) %>%  # 从大到小排序
    slice_head(n = top_n) %>%
    ungroup()
  
  # 创建综合特征重要性图
  p_imp <- ggplot(avg_importance, aes(x = mean_importance, y = reorder(feature, mean_importance), fill = analysis)) +
    geom_col(alpha = 0.8) +
    facet_wrap(~ analysis, scales = "free_y", ncol = 2) +
    labs(
      title = paste("Top", top_n, "特征重要性"),
      x = "平均重要性得分",
      y = "微生物特征",
      fill = "分析类型"
    ) +
    theme_test() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    ) +
    scale_fill_brewer(palette = "Set1")
  
  return(p_imp)
}

# 执行NC vs ON分析
nc_on_results <- binary_loocv_analysis(data, "NC", "ON", "NC_vs_ON")

# 执行NC vs OP分析
nc_op_results <- binary_loocv_analysis(data, "NC", "OP", "NC_vs_OP")

# 创建综合图
cat("\n=== 创建综合图 ===\n")

# 创建综合ROC图
combined_roc_plot <- create_combined_roc_plot(nc_on_results, nc_op_results)

# 创建综合特征重要性图
combined_imp_plot <- create_combined_importance_plot(nc_on_results, nc_op_results, top_n = 5)

# 使用patchwork组合图片
if (!is.null(combined_roc_plot) && !is.null(combined_imp_plot)) {
  # 垂直组合图片
  final_plot <- combined_roc_plot / combined_imp_plot +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(face = 'bold', size = 14))
  
  # 保存综合图
  ggsave("Combined_Analysis_Results.png", final_plot, 
         width = 14, height = 16, dpi = 300, bg = "white")
  cat("综合图已保存为 Combined_Analysis_Results.png\n")
  
  # 显示图片
  print(final_plot)
  
} else if (!is.null(combined_roc_plot)) {
  # 只有ROC图
  ggsave("Combined_ROC_Analysis.png", combined_roc_plot, 
         width = 12, height = 8, dpi = 300, bg = "white")
  cat("ROC综合图已保存为 Combined_ROC_Analysis.png\n")
  print(combined_roc_plot)
  
} else if (!is.null(combined_imp_plot)) {
  # 只有特征重要性图
  ggsave("Combined_Importance_Analysis.png", combined_imp_plot, 
         width = 10, height = 8, dpi = 300, bg = "white")
  cat("特征重要性综合图已保存为 Combined_Importance_Analysis.png\n")
  print(combined_imp_plot)
}

# 输出详细的AUC结果
cat("\n=== 详细AUC结果汇总 ===\n")

# NC vs ON 结果汇总
if (!is.null(nc_on_results) && length(nc_on_results$roc_results) > 0) {
  cat("\n*** NC vs ON 结果 ***\n")
  for (i in 1:length(nc_on_results$roc_results)) {
    if (!is.null(nc_on_results$roc_results[[i]])) {
      result <- nc_on_results$roc_results[[i]]
      cat("队列:", result$cohort, 
          " | AUC:", round(result$auc, 3),
          " | 95%CI: [", round(result$ci_lower, 3), "-", round(result$ci_upper, 3), "]",
          " | 发现集样本:", paste(names(result$discovery_counts), result$discovery_counts, collapse = "/"),
          " | 验证集样本:", paste(names(result$validation_counts), result$validation_counts, collapse = "/"),
          "\n")
    }
  }
} else {
  cat("\n*** NC vs ON: 无有效结果 ***\n")
}

# NC vs OP 结果汇总
if (!is.null(nc_op_results) && length(nc_op_results$roc_results) > 0) {
  cat("\n*** NC vs OP 结果 ***\n")
  for (i in 1:length(nc_op_results$roc_results)) {
    if (!is.null(nc_op_results$roc_results[[i]])) {
      result <- nc_op_results$roc_results[[i]]
      cat("队列:", result$cohort, 
          " | AUC:", round(result$auc, 3),
          " | 95%CI: [", round(result$ci_lower, 3), "-", round(result$ci_upper, 3), "]",
          " | 发现集样本:", paste(names(result$discovery_counts), result$discovery_counts, collapse = "/"),
          " | 验证集样本:", paste(names(result$validation_counts), result$validation_counts, collapse = "/"),
          "\n")
    }
  }
} else {
  cat("\n*** NC vs OP: 无有效结果 ***\n")
}

# 计算平均AUC
calculate_mean_auc <- function(roc_results) {
  if (is.null(roc_results)) return(NA)
  aucs <- sapply(roc_results, function(x) if(!is.null(x)) x$auc else NA)
  aucs <- aucs[!is.na(aucs)]
  if (length(aucs) > 0) {
    return(round(mean(aucs), 3))
  } else {
    return(NA)
  }
}

cat("\n=== 平均AUC ===\n")
cat("NC vs ON 平均AUC:", calculate_mean_auc(nc_on_results$roc_results), "\n")
cat("NC vs OP 平均AUC:", calculate_mean_auc(nc_op_results$roc_results), "\n")

# 保存所有结果
save(nc_on_results, nc_op_results, file = "binary_classification_results.RData")

cat("\n分析完成！\n")
cat("生成的文件:\n")
if (!is.null(combined_roc_plot) && !is.null(combined_imp_plot)) {
  cat("- Combined_Analysis_Results.png: 综合分析结果图\n")
} else if (!is.null(combined_roc_plot)) {
  cat("- Combined_ROC_Analysis.png: ROC综合分析图\n")
} else if (!is.null(combined_imp_plot)) {
  cat("- Combined_Importance_Analysis.png: 特征重要性综合图\n")
}
cat("- binary_classification_results.RData: 完整结果数据\n")
# 加载必要的包
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("vegan")) install.packages("vegan")
if (!require("ape")) install.packages("ape")
if (!require("patchwork")) install.packages("patchwork")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("rstatix")) install.packages("rstatix")

library(tidyverse)
library(vegan)
library(ape)
library(patchwork)
library(ggrepel)
library(rstatix)

# 设置随机种子
set.seed(123)

# 1. 读取数据 ----------------------------------------------------
cat("Reading microbial composition data...\n")
microbial_data <- read.table("123.txt", header = TRUE, sep = "\t", 
                             stringsAsFactors = FALSE, check.names = FALSE)
cat("Data dimensions:", dim(microbial_data), "\n")

# 2. 数据预处理 --------------------------------------------------
# 提取样本信息
sample_info <- microbial_data[, 1:3]
colnames(sample_info) <- c("Sample_ID", "Cohort", "Group")

# 提取微生物特征
microbial_features <- microbial_data[, -(1:3)]
feature_matrix <- as.matrix(microbial_features)
rownames(feature_matrix) <- sample_info$Sample_ID

# 转换为相对丰度
feature_matrix_rel <- feature_matrix / rowSums(feature_matrix)

# 3. 图1: 平均相对丰度前十微生物堆叠柱状图 ----------------------
cat("\n1. Creating stacked bar plot for top 10 microbes...\n")

# 计算每个微生物的平均相对丰度
mean_abundance <- apply(feature_matrix_rel, 2, mean)
top10_microbes <- names(sort(mean_abundance, decreasing = TRUE))[1:10]

# 准备堆叠柱状图数据
stack_data <- as.data.frame(feature_matrix_rel) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(sample_info, by = "Sample_ID") %>%
  select(Sample_ID, Group, all_of(top10_microbes)) %>%
  pivot_longer(cols = all_of(top10_microbes), 
               names_to = "Microbe", 
               values_to = "Abundance") %>%
  group_by(Group, Microbe) %>%
  summarise(Mean_Abundance = mean(Abundance), .groups = "drop")

# 计算Others
others_data <- stack_data %>%
  group_by(Group) %>%
  summarise(Group_Abundance = sum(Mean_Abundance)) %>%
  mutate(Others = 1 - Group_Abundance) %>%
  select(Group, Others) %>%
  pivot_longer(cols = Others, names_to = "Microbe", values_to = "Mean_Abundance")

# 合并数据
final_stack_data <- bind_rows(stack_data, others_data)

# 简化微生物名称（取属水平）
final_stack_data$Microbe_short <- gsub(".*;g__", "", final_stack_data$Microbe)
final_stack_data$Microbe_short <- ifelse(final_stack_data$Microbe_short == "", 
                                         "Unclassified", final_stack_data$Microbe_short)

# 按照平均丰度从大到小排序
microbe_order <- final_stack_data %>%
  filter(Microbe != "Others") %>%
  group_by(Microbe_short) %>%
  summarise(Total_Abundance = sum(Mean_Abundance)) %>%
  arrange(desc(Total_Abundance)) %>%
  pull(Microbe_short)

# 添加Others到最后
microbe_order <- c(microbe_order, "Others")

# 设置因子顺序
final_stack_data$Group <- factor(final_stack_data$Group, levels = c("NC", "ON", "OP"))
final_stack_data$Microbe_short <- factor(final_stack_data$Microbe_short, levels = microbe_order)

# 创建堆叠柱状图
plot1 <- ggplot(final_stack_data, aes(x = Mean_Abundance, y = Group, fill = Microbe_short)) +
  geom_col(position = "fill", width = 0.7) +
  labs(title = "Top 10 Microbial Genera by Average Relative Abundance",
       x = "Relative Abundance",
       y = "Group",
       fill = "Microbial Genus") +
  scale_x_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set3") +
  theme_test() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# 4. 图2: Alpha多样性分析 --------------------------------------
cat("2. Calculating alpha diversity and creating boxplots...\n")

# 计算Alpha多样性
alpha_diversity <- data.frame(
  Sample_ID = rownames(feature_matrix_rel),
  Shannon = diversity(feature_matrix_rel, index = "shannon"),
  Simpson = diversity(feature_matrix_rel, index = "simpson")
) %>% left_join(sample_info, by = "Sample_ID")

# 转换成长格式用于绘图
alpha_long <- alpha_diversity %>%
  pivot_longer(cols = c(Shannon, Simpson), 
               names_to = "Index", 
               values_to = "Value")

# 在每个队列内进行统计检验的函数
perform_cohort_alpha_stats <- function(data, cohort_name) {
  cohort_data <- data %>% filter(Cohort == cohort_name)
  unique_groups <- unique(cohort_data$Group)
  
  if (length(unique_groups) < 2) {
    return(NULL)  # 如果只有一组，不进行检验
  }
  
  results <- data.frame()
  
  for (index_type in c("Shannon", "Simpson")) {
    index_data <- cohort_data %>% filter(Index == index_type)
    
    if (length(unique_groups) == 2) {
      # Wilcoxon检验
      test_result <- wilcox_test(index_data, Value ~ Group)
      p_val <- test_result$p
      test_type <- "Wilcoxon"
    } else {
      # Kruskal-Wallis检验
      test_result <- kruskal_test(index_data, Value ~ Group)
      p_val <- test_result$p
      test_type <- "K-W"
    }
    
    results <- rbind(results, data.frame(
      Cohort = cohort_name,
      Index = index_type,
      p_value = p_val,
      test_type = test_type
    ))
  }
  
  return(results)
}

# 对每个队列执行统计检验
all_cohorts <- unique(alpha_long$Cohort)
stats_results <- data.frame()

for (cohort in all_cohorts) {
  cohort_stats <- perform_cohort_alpha_stats(alpha_long, cohort)
  if (!is.null(cohort_stats)) {
    stats_results <- rbind(stats_results, cohort_stats)
  }
}

# 创建Alpha多样性图
plot2 <- ggplot(alpha_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 0.8) +
  facet_grid(Index ~ Cohort, scales = "free_y") +
  labs(title = "Alpha Diversity by Group and Cohort",
       x = "Group",
       y = "Diversity Index Value",
       fill = "Group") +
  scale_fill_brewer(palette = "Set2") +
  theme_test() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  # 添加统计注释
  geom_text(data = stats_results,
            aes(x = 1.5, y = Inf, 
                label = paste0(test_type, "\np = ", round(p_value, 4))),
            vjust = 1.5, size = 3, inherit.aes = FALSE)

# 5. 图3: PCoA分析 ---------------------------------------------
cat("3. Performing PCoA analysis...\n")

# 计算Bray-Curtis距离
bray_dist <- vegdist(feature_matrix_rel, method = "bray")

# PCoA分析
pcoa_result <- pcoa(bray_dist)
variance_explained <- pcoa_result$values$Relative_eig * 100

# PERMANOVA检验
adonis_result <- adonis2(bray_dist ~ Group, data = sample_info)
adonis_r2 <- adonis_result$R2[1]
adonis_pvalue <- adonis_result$`Pr(>F)`[1]

# 准备PCoA绘图数据
pcoa_df <- data.frame(
  PCo1 = pcoa_result$vectors[, 1],
  PCo2 = pcoa_result$vectors[, 2],
  Group = sample_info$Group
)

# 创建PCoA图
plot3 <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2, color = Group)) +
  geom_point(alpha = 0.7, size = 2) +
  stat_ellipse(level = 0.68, alpha = 0.3, linewidth = 0.8) +
  labs(title = "PCoA Plot of Microbial Communities",
       subtitle = paste0("PERMANOVA: R² = ", round(adonis_r2, 3), 
                         ", p = ", ifelse(adonis_pvalue < 0.001, "< 0.001", 
                                          round(adonis_pvalue, 4))),
       x = paste0("PCo1 (", round(variance_explained[1], 1), "%)"),
       y = paste0("PCo2 (", round(variance_explained[2], 1), "%)"),
       color = "Group") +
  theme_test() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_color_brewer(palette = "Set1")

# 6. 图4: 差异微生物分析火山图 --------------------------------
cat("4. Performing pairwise differential abundance analysis...\n")

#  pairwise Wilcoxon检验函数
pairwise_wilcox_test_microbes <- function(data, group1, group2) {
  results <- data.frame(Microbe = character(),
                        p_value = numeric(),
                        log2FC = numeric(),
                        stringsAsFactors = FALSE)
  
  # 筛选两组数据
  group1_data <- data[sample_info$Group == group1, ]
  group2_data <- data[sample_info$Group == group2, ]
  
  for (microbe in colnames(data)) {
    # 只对在至少一组中有丰度的微生物进行检验
    if (sum(group1_data[, microbe] > 0) > 1 | sum(group2_data[, microbe] > 0) > 1) {
      test_result <- wilcox.test(group1_data[, microbe], group2_data[, microbe])
      
      # 计算log2 Fold Change
      mean1 <- mean(group1_data[, microbe])
      mean2 <- mean(group2_data[, microbe])
      log2FC <- log2((mean2 + 1e-10) / (mean1 + 1e-10))  # 添加伪计数避免除零
      
      results <- rbind(results, data.frame(
        Microbe = microbe,
        p_value = test_result$p.value,
        log2FC = log2FC,
        Comparison = paste0(group2, "_vs_", group1)
      ))
    }
  }
  
  results$p_adjust <- p.adjust(results$p_value, method = "fdr")
  return(results)
}

# 执行三组比较
comparisons <- list(
  c("NC", "ON"),
  c("NC", "OP"), 
  c("ON", "OP")
)

all_pairwise_results <- data.frame()

for (comp in comparisons) {
  comp_results <- pairwise_wilcox_test_microbes(as.data.frame(feature_matrix_rel), comp[1], comp[2])
  all_pairwise_results <- rbind(all_pairwise_results, comp_results)
}

# 简化微生物名称
all_pairwise_results$Microbe_short <- gsub(".*;g__", "", all_pairwise_results$Microbe)
all_pairwise_results$Microbe_short <- ifelse(all_pairwise_results$Microbe_short == "", 
                                             "Unclassified", all_pairwise_results$Microbe_short)

# 标记显著性
all_pairwise_results$Significance <- case_when(
  all_pairwise_results$p_adjust < 0.05 & abs(all_pairwise_results$log2FC) > 1 ~ "Significant (FDR<0.05, |FC|>2)",
  all_pairwise_results$p_adjust < 0.05 ~ "Significant (FDR<0.05)",
  TRUE ~ "Not significant"
)

# 创建火山图
plot4 <- ggplot(all_pairwise_results, aes(x = log2FC, y = -log10(p_adjust), color = Significance)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_text_repel(
    data = subset(all_pairwise_results, Significance != "Not significant"),
    aes(label = Microbe_short),
    size = 2.5,
    max.overlaps = 15,
    box.padding = 0.3
  ) +
  facet_wrap(~ Comparison, ncol = 2) +
  labs(title = "Pairwise Differential Abundance Analysis",
       x = "log2(Fold Change)",
       y = "-log10(Adjusted P-value)",
       color = "Significance") +
  scale_color_manual(values = c("Significant (FDR<0.05, |FC|>2)" = "red", 
                                "Significant (FDR<0.05)" = "orange", 
                                "Not significant" = "grey60")) +
  theme_test() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 7. 保存所有图形 ----------------------------------------------
cat("5. Saving all plots...\n")

# 单独保存每个图
ggsave("stacked_barplot_top10.png", plot1, width = 10, height = 6, dpi = 300)
ggsave("alpha_diversity.png", plot2, width = 12, height = 8, dpi = 300)
ggsave("pcoa_plot.png", plot3, width = 8, height = 6, dpi = 300)
ggsave("volcano_plot.png", plot4, width = 12, height = 8, dpi = 300)

# 组合图
combined_plot <- (plot1 | plot3) / (plot2 | plot4) + 
  plot_annotation(tag_levels = 'A') +
  plot_layout(heights = c(1, 1.5))

ggsave("combined_analysis.png", combined_plot, width = 16, height = 12, dpi = 300)

# 8. 输出统计结果 ----------------------------------------------
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("1. Top 10 microbes identified and visualized (sorted by abundance)\n")
cat("2. Alpha diversity statistics by cohort:\n")
print(stats_results)
cat("3. PCoA PERMANOVA results:\n")
cat("   R² =", round(adonis_r2, 4), "p =", 
    ifelse(adonis_pvalue < 0.001, "< 0.001", round(adonis_pvalue, 4)), "\n")
cat("4. Pairwise differential abundance analysis:\n")
for (comp in comparisons) {
  comp_name <- paste0(comp[2], "_vs_", comp[1])
  comp_data <- all_pairwise_results %>% filter(Comparison == comp_name)
  cat("   ", comp_name, ": Significant microbes (FDR < 0.05) =", 
      sum(comp_data$Significance != "Not significant"), "\n")
}

# 保存显著差异微生物列表
significant_microbes <- all_pairwise_results %>%
  filter(Significance != "Not significant") %>%
  select(Comparison, Microbe, Microbe_short, log2FC, p_value, p_adjust, Significance)

write.table(significant_microbes, "significant_microbes_pairwise.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nResults saved:\n")
cat("- stacked_barplot_top10.png\n")
cat("- alpha_diversity.png\n") 
cat("- pcoa_plot.png\n")
cat("- volcano_plot.png\n")
cat("- combined_analysis.png\n")
cat("- significant_microbes_pairwise.txt\n")

cat("\n=== ANALYSIS COMPLETED ===\n")
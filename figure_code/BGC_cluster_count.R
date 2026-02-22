#load package  加载需要的包
library(dplyr)
library(ggplot2)
library(eoffice)
library(ggplot2)
library(tidyr)
library(ggsci)

#Set working directory & read data  设置路径，读取数据
setwd("D:/metagenome/result/2025_data_analysis/Plot")
df <-  read.csv("../all_BGC_result.csv",check.names = FALSE)
bgc_class <- read.csv("../new_compound_BGC_cluster.csv")

# ------------------------- Step 1: Define columns 定义列范围 -------------------------
# Define the range of BGC abundance columns 定义 BGC 丰度所在的列范围
bgc_cols <- 2:(ncol(df)-7)
# Define the classification column (environment type) 定义分类列（环境类型）
class_col <- "Category"   



# ------------------------- Step 2: Convert abundance matrix to long format 转换为长表 -------------------------
df_long <- df %>%
  pivot_longer(cols = all_of(bgc_cols), 
               names_to = "BGC",       # Column name for BGC names 新列存放 BGC 名称
               values_to = "Abundance" # Column name for abundance 新列存放丰度
  )

# ------------------------- Step 3: Filter BGCs with abundance > 0 筛选丰度大于 0 的 BGC -------------------------
df_detected <- df_long %>%
  filter(as.numeric(Abundance) > 0)


# ------------------------- Step 4: Join with BGC category annotation 加入 BGC 分类信息 -------------------------
df_detected <- df_detected %>%
  left_join(bgc_class, by = c("BGC" = "BGC_name"))


# ------------------------- Step 5: Count number of BGCs per category per sample 统计每个样本每类 BGC 的数量 -------------------------
df_summary <- df_detected %>%
  group_by(SampleID = Run, 
           !!sym(class_col),   
           Category,           
           BGC_cluster) %>%    
  summarise(Count = n(), .groups = "drop")

# ------------------------- Step 6: Plot boxplot 按环境绘制箱线图 -------------------------
df_summary <- df_summary %>%
  mutate(BGC_cluster = ifelse(BGC_cluster == "Ribosomally synthesized and post-translationally modified peptide",
                              "RiPPs",
                              ifelse(BGC_cluster == "Nonribosomal peptide",
                                     "NRPs",
                                     BGC_cluster)))

# This ensures the x-axis order follows the given sequence 确定绘图时的 x 轴顺序为指定顺序

df_summary$BGC_cluster <- factor(df_summary$BGC_cluster,
                                  levels = c("Alkaloid","NRPs","Polyketide","RiPPs","Saccharide","Terpene","Other"))


p <- ggplot(df_summary, aes(x = BGC_cluster, 
                            y = Count, 
                            fill = BGC_cluster)) +  # 使用 BGC_cluster 作为填充颜色
  geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = BGC_cluster),  # 使用 BGC_cluster 作为散点颜色
              size = 1, alpha = 0.6,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  scale_fill_npg() +  # 应用 'npg' 颜色方案
  scale_color_npg() +  # 应用 'npg' 颜色方案
  theme_bw(base_size = 14) +
  labs(x = "BGC Cluster",
       y = "Number of BGCs detected") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Category, scales = "free_x")  # 按环境分面

# -------------------------  Step 7: Save the plot 保存文件 -------------------------
topptx(p,"BGC_cluster_count.pptx", width = 10, height = 6)
ggsave("BGC_cluster_count.pdf", plot = p, width = 10, height = 6)


# -------------------------  Step 8: Plot stacked bar 绘制堆叠柱状图-------------------------


# Aggregate per-sample counts to environment x cluster statistics -------------------------
# 按样本汇总并计算每个环境不同BGC类别的 mean 和 median -------------------------
env_col <- "Category"
env_cluster_stats <- df_summary %>%
  # ensure columns exist and are of correct types
  mutate(!!env_col := as.character(!!sym(env_col)),
         BGC_cluster = as.character(BGC_cluster),
         Count = as.numeric(Count)) %>%
  # group by environment and cluster, then compute mean and median across samples
  group_by(!!sym(env_col), BGC_cluster) %>%
  summarise(
    mean_count = mean(Count, na.rm = TRUE),
    median_count = median(Count, na.rm = TRUE),
    .groups = "drop"
  )

# Convert to proportions (per environment)  按环境把 mean/median 转为比例（占比）
env_cluster_props <- env_cluster_stats %>%
  group_by(!!sym(env_col)) %>%
  mutate(
    mean_prop   = mean_count / sum(mean_count, na.rm = TRUE),
    median_prop = median_count / sum(median_count, na.rm = TRUE)
  ) %>%
  ungroup()

# This ensures the  order follows the given sequence 确定绘图时的顺序为指定顺序

env_cluster_props$BGC_cluster <- factor(env_cluster_props$BGC_cluster,
                                 levels = c("Alkaloid","NRPs","Polyketide","RiPPs","Saccharide","Terpene","Other"))


env_cluster_props$Category <- factor(env_cluster_props$Category,
                                  levels = c("Sediment", "human", "water","Soil","Rhizosphere"))



p_mean <- ggplot(env_cluster_props, aes_string(x ="mean_prop" , y = env_col, fill = "BGC_cluster")) +
  geom_bar(stat = "identity", position = "fill") +  # position="fill" makes 100% stacked (proportion)
  # note: since we've precomputed proportions, you could use stat="identity" and y=mean_prop with position="stack"
  # using position="fill" plus mean_prop will still work and scale to 1 per env
  labs(x = "Environment",
       y = "Proportion (mean)",
       title = "Proportion of BGC clusters per Environment (mean)") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_npg() +  # 应用 'npg' 颜色方案
  scale_color_npg()   # 应用 'npg' 颜色方案

p_mean 

p_median <- ggplot(env_cluster_props, aes_string(x = "median_prop", y = env_col, fill = "BGC_cluster")) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Environment",
       y = "Proportion (median)",
       title = "Proportion of BGC clusters per Environment (median)") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_npg() +  # 应用 'npg' 颜色方案
  scale_color_npg()   # 应用 'npg' 颜色方案
p_median


# -------------------------  Step 9: Save the plot 保存文件 -------------------------
topptx(p_mean,"BGC_cluster_stacked_bar_mean.pptx", width = 8, height = 6)
ggsave("BGC_cluster_stacked_bar_mean.pdf", plot = p_mean, width = 8, height = 6)

topptx(p_median,"BGC_cluster_stacked_bar_median.pptx", width = 8, height = 6)
ggsave("BGC_cluster_stacked_bar_median.pdf", plot = p_median, width = 8, height = 6)

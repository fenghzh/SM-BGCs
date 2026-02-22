# ------------------------- Load packages 加载需要的包 -------------------------
library(dplyr)
library(Rtsne) 
library(ggplot2) 


#Set working directory & read data  设置路径，读取数据
setwd("D:/metagenome/result/2025_data_analysis/Plot")
df <-  read.csv("../all_BGC_result.csv")
rownames(df) <- df$Run
# ------------------------- Step 1: Separate abundance and metadata 分离丰度矩阵和注释信息 -------------------------
# Define the range of BGC abundance columns 定义 BGC 丰度所在的列范围
bgc_cols <- 2:(ncol(df)-7)
bgc_matrix <- df[, bgc_cols]
sample_info <- df[,c(1,(ncol(df)-6):ncol(df))]

#Make sure it is a numeric matrix  确保是数值矩阵
bgc_matrix <- as.matrix(bgc_matrix)

# Remove zero-only columns 删除全 0 的列
bgc_matrix <- bgc_matrix[, colSums(bgc_matrix) > 0]

# Remove duplicate rows 删除重复行
bgc_matrix <- unique(bgc_matrix)


# ------------------------- Step 2: Run t-SNE 运行 t-SNE -------------------------
set.seed(123) # 设置随机种子，保证结果可复现
tsne_out <- Rtsne(bgc_matrix, dims = 2, pca = FALSE,
                  max_iter = 1000, theta = 0.0,
                  perplexity = 20, verbose = FALSE)

# Prepare t-SNE result dataframe 整理 t-SNE 结果
tsne_result <- as.data.frame(tsne_out$Y)
colnames(tsne_result) <- c("tSNE1", "tSNE2")

# Merge with metadata 合并样本注释信息
tsne_result$Run <- rownames(bgc_matrix)
tsne_result <-  tsne_result  %>%   left_join(sample_info, by = c("Run" = "Run"))

# ------------------------- Step 3: Plot t-SNE 绘制 t-SNE 图 -------------------------
# Define colors 设置分类颜色 
sample_color <- data.frame(
  Category = unique(tsne_result$Category),
  colors   = c("#ccb994", "#8a0089", "#0000e2", "#8a460f", "#008b01")
)

p <- ggplot(tsne_result, aes(x = tSNE1, y = tSNE2, color = Category)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(name = "Category",
                     values = setNames(sample_color$colors, sample_color$Category)) +
  scale_shape_manual(values = c(15, 19), limits = c("saline", "non-saline")) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "t-SNE of BGC Abundance",
       x = "t-SNE 1",
       y = "t-SNE 2")
p
# -------------------------  Step 4: Save the plot 保存文件 -------------------------
topptx(p,"tSNE_BGC.pptx", width = 8, height = 6)
ggsave("tSNE_BGC.pdf", plot = p, width = 8, height = 6)


#For microbial composition data, you only need to import the corresponding table for analysis
#对于微生物组成数据，只需要导入对应的表格进行分析
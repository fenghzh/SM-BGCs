#load package  加载需要的包
library(dplyr)
library(ggplot2)
library(readr)
library(eoffice)

#Set working directory & read data  设置路径，读取数据
setwd("D:/metagenome/result/2025_data_analysis/Plot")
df <-  read.csv("../all_BGC_result.csv")



# ------------------------- Step 1: Define columns 定义列范围 -------------------------
# Define the range of BGC abundance columns 定义 BGC 丰度所在的列范围
bgc_cols <- 2:(ncol(df)-7)
# Define the classification column (environment type) 定义分类列（环境类型）
class_col <- "Category"   

# ------------------------- Step 2: Count detected BGCs per sample 统计每个样本的 BGC 个数 -------------------------
df <- df %>% mutate(across(bgc_cols, as.numeric))
df_summary <- df %>%
  rowwise() %>%
  mutate(BGC_count = sum(c_across(all_of(bgc_cols)) > 0, na.rm = TRUE)) %>%
  ungroup()

# ------------------------- Step 3: Set factor levels for classification 设置分类顺序 -------------------------


#Specify colors for each category  指定每个分类的颜色 
sample_colors <- c(
  "Rhizosphere" = "#008b01",
  "Soil"       = "#8a460f",
  "water"      = "#0000e2",
  "human"      = "#8a0089",
  "Sediment"   = "#ccb994"
)

# This ensures the x-axis order follows the given sequence 确定绘图时的 x 轴顺序为指定顺序

df_summary[[class_col]] <- factor(df_summary[[class_col]],
                                  levels = c("Rhizosphere", "Soil", "water", "human", "Sediment"))


# ------------------------- Step 4: Plot boxplot 绘制箱线图 -------------------------
p <-  ggplot(df_summary, aes(x = Category, y = BGC_count, color=Category)) +
      stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(0.8)) +
      geom_boxplot(position = position_dodge(0.8), alpha = 0.3, outlier.colour = NA) +
      geom_point(position = position_jitter(width = 0.2, height = 0), alpha =1,size=1) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(face = "plain", color = "black", 
                                   size = 10, angle = 90, hjust = 1),
          axis.text.y = element_text(face = "plain", color = "black", 
                                   size = 10),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank())  +
          scale_color_manual(values = sample_colors) 
p
# -------------------------  Step 5: Save the plot 保存文件 -------------------------
topptx(p,"BGC_count_boxplot.pptx", width = 7, height = 5)


#For microbial composition data, you only need to import the corresponding table for analysis
#对于微生物组成数据，只需要导入对应的表格进行分析
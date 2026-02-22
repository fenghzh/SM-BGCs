library(tidymodels)  # 加载 tidymodels 框架 (load tidymodels framework)
library(discrim)
library(klaR)
library(data.table)
library(doParallel)
library(foreach)
library(vip)
library(eoffice)
setwd("D:/metagenome/result/2025_data_analysis/Plot")


current_threads <- getDTthreads()
print(paste("当前data.table线程数:", current_threads))

# 设置data.table使用4个线程
setDTthreads(threads = 4)
print(paste("设置后data.table线程数:", getDTthreads()))

# 2. 使用doParallel设置并行计算
# 检测系统可用核心数
num_cores <- detectCores()
print(paste("可用核心数:", num_cores))

# 设置使用4个核心进行并行计算
cl <- makeCluster(4)
registerDoParallel(cl)




# 假设 bgc_data 数据框已读取并包含样本的BGC丰度特征和 Category 列
# 构建预处理配方 (create a recipe for preprocessing)
bgc_data <- read.csv("../all_BGC_result.csv",check.names = FALSE)
rownames(bgc_data) <- bgc_data[,1]
bgc_data <- bgc_data[,2:(ncol(bgc_data)-6)]
bgc_data$Category <- as.factor(bgc_data$Category)


#bgc_data <- bgc_data %>%
#mutate(across(2:(ncol(bgc_data)-6), ~ ifelse(. == 0, 0.0001, .)))
bgc_data <- bgc_data %>%
  mutate(across(2:(ncol(bgc_data) - 6), ~ log1p(.)))

rec <- recipe(Category ~ ., data = bgc_data) 

#  %>%
#  step_zv(all_predictors()) %>%                # 去除零方差特征 (remove zero-variance predictors)
#  step_nzv(all_predictors()) %>%               # 去除近零方差特征 (remove near-zero variance predictors)
#  step_normalize(all_numeric_predictors())     # 标准化数值特征 (zero mean, unit variance)


# 定义随机森林模型 (Random Forest)
rf_spec <- rand_forest(mode = "classification", trees = 500) %>%
  set_engine("ranger", importance = "impurity")  # importance="impurity" 以计算变量重要性:contentReference[oaicite:7]{index=7}

wf_rf <- workflow() %>% 
  add_model(rf_spec) %>% 
  add_recipe(rec)

# 定义 XGBoost 模型 (XGBoost)
xgb_spec <- boost_tree(mode = "classification", trees = 500, tree_depth = 6, learn_rate = 0.1) %>%
  set_engine("xgboost")

wf_xgb <- workflow() %>%
  add_model(xgb_spec) %>%
  add_recipe(rec)

# 定义支持向量机 SVM 模型(RBF kernel)
svm_spec <- svm_rbf(mode = "classification") %>%
  set_engine("kernlab")

wf_svm <- workflow() %>%
  add_model(svm_spec) %>% 
  add_recipe(rec)

# 生成10折交叉验证分组 (generate 10-fold CV splits)
set.seed(123)
cv_splits <- vfold_cv(bgc_data, v = 10, strata = Category)


# 定义评价指标 (define metrics)
metrics <- metric_set(accuracy, sensitivity, specificity)  # 同时计算准确率、敏感度、特异性:contentReference[oaicite:12]{index=12}

# 交叉验证训练 (CV resampling)
rf_res   <- fit_resamples(wf_rf,   cv_splits, metrics = metrics)
xgb_res  <- fit_resamples(wf_xgb,  cv_splits, metrics = metrics)
svm_res <- fit_resamples(wf_svm, cv_splits, metrics = metrics)

# 收集并汇总指标 (collect metrics)
rf_metrics  <- collect_metrics(rf_res)   %>% mutate(Model = "Random Forest")
xgb_metrics <- collect_metrics(xgb_res)  %>% mutate(Model = "XGBoost")
svm_metrics <- collect_metrics(svm_res) %>% mutate(Model = "SVM")

all_metrics <- bind_rows(rf_metrics, xgb_metrics, svm_metrics)





# 8. 绘制模型准确率比较图
# ==============================================
acc_data <- all_metrics %>% filter(.metric == "accuracy")

p <- ggplot(acc_data, aes(x = Model, y = mean, fill = Model)) +
  geom_col(width = 0.5) +
  labs(title = "模型准确率比较 (Accuracy Comparison)",
       y = "平均准确率 (Mean Accuracy)",
       x = "模型 (Model)") +
  theme_bw()
p
topptx(p,"BGC_Accuracy Comparison.pptx", width =6, height = 5)

# ==============================================
# 9. 拟合最终模型用于特征重要性分析 (RF 和 XGB)
# ==============================================
rf_fit  <- fit(wf_rf,  data = bgc_data)
xgb_fit <- fit(wf_xgb, data = bgc_data)
svm_fit <- fit(wf_svm, data = bgc_data)


# 9.1 随机森林变量重要性
p <- vip(rf_fit %>% extract_fit_parsnip(), num_features = 20) +
  ggtitle("随机森林前10个重要BGC特征 (Top 20 RF Features)")
p
topptx(p,"RF_TOP_20.pptx", width = 6, height = 5)

# 9.2 XGBoost变量重要性
p <- vip(xgb_fit %>% extract_fit_parsnip(), num_features = 20) +
  ggtitle("XGBoost前10个重要BGC特征 (Top 20 XGB Features)")
p
topptx(p,"xgb_TOP_20.pptx", width = 6, height = 5)


# ==============================================
# 10. 灵敏度和特异性比较
# ==============================================
# 选择 sensitivity 和 specificity
other_metrics <- all_metrics %>% filter(.metric %in% c("sensitivity", "specificity"))

p <- ggplot(other_metrics, aes(x = Model, y = mean, fill = .metric)) +
  geom_col(position = "dodge",width=0.5) +
  labs(title = "模型敏感性和特异性比较 (Sensitivity & Specificity Comparison)",
       y = "平均值 (Mean Value)",
       x = "模型 (Model)") +
  theme_bw()
p
topptx(p,"BGC_Sensitivity_Specificity.pptx", width = 6, height = 5)


# ==============================================
# 11. 混淆矩阵可视化
# ==============================================
library(yardstick)

# RF 混淆矩阵
rf_aug <- augment(rf_fit, new_data = bgc_data)
rf_aug$.pred_class <- as.factor(rf_aug$.pred_class)

rf_cm <- conf_mat(rf_aug, truth = Category, estimate = .pred_class)

# 绘制热力图
p <- autoplot(rf_cm, type = "heatmap") +
  ggtitle("随机森林混淆矩阵 (RF Confusion Matrix)")
p
topptx(p,"BGC_RF_Confusion_Matrix.pptx", width = 6, height = 5)

# XGB 混淆矩阵
xgb_aug <- augment(xgb_fit, new_data = bgc_data)
xgb_aug$.pred_class <- as.factor(xgb_aug$.pred_class)

xgb_cm <- conf_mat(xgb_aug, truth = Category, estimate = .pred_class)

# 绘制热力图
p <- autoplot(xgb_cm, type = "heatmap") +
  ggtitle("XGBoost混淆矩阵 (XGB Confusion Matrix)")
p
topptx(p,"BGC_XGB_Confusion_Matrix.pptx", width = 6, height = 5)

# SVM 混淆矩阵
svm_aug <- augment(svm_fit, new_data = bgc_data)
svm_aug$.pred_class <- as.factor(svm_aug$.pred_class)

svm_cm <- conf_mat(svm_aug, truth = Category, estimate = .pred_class)

# 绘制热力图
p <- autoplot(svm_cm, type = "heatmap") +
  ggtitle("支持向量机混淆矩阵 (SVM Confusion Matrix)")
p
topptx(p,"BGC_SVM_Confusion_Matrix.pptx", width = 6, height = 5)


# ==============================================http://127.0.0.1:27122/graphics/plot_zoom_png?width=1168&height=822
# 12. 特征重要性图
# ==============================================
library(vip)

# RF VIP
vip(rf_fit %>% extract_fit_parsnip(), num_features = 20) +
  ggtitle("随机森林前20个重要BGC特征 (Top 10 RF Features)")

# XGB VIP
vip(xgb_fit %>% extract_fit_parsnip(), num_features = 20) +
  ggtitle("XGBoost前10个重要BGC特征 (Top 10 XGB Features)")



#
##
##
##重要特征与环境的相关性
rf_vi <- vip::vi(xgb_fit %>% extract_fit_parsnip())
xgb_vi <- vip::vi(xgb_fit %>% extract_fit_parsnip())
rf_vi <- rf_vi$Variable[1:20]
xgb_vi <- xgb_vi$Variable[1:20]
common_elements <- intersect(rf_vi, xgb_vi )

X <- as.data.frame(unique(c(rf_vi,xgb_vi)))
top_features <- rf_vi$Variable[1:20]


# 假设 df_top 是你的前10重要特征 + Category
df_top <- bgc_data %>%
  dplyr::select(all_of(top_features), Category) %>%
  mutate(Category = as.factor(Category))

#对 BGC 丰度做 log1p 转换



categories <- levels(df_top$Category)

res_list <- list()

# 计算每个特征与每个类别的点二列相关系数
for(cat in categories){
  df_top$Category_bin <- ifelse(df_top$Category == cat, 1, 0)
  
  corr <- df_top %>%
    summarise(across(all_of(top_features),
                     ~ cor(.x, Category_bin, method = "spearman"),
                     .names = "{.col}")) %>%
    pivot_longer(cols = everything(), names_to = "Feature", values_to = "Corr") %>%
    mutate(Category = cat)
  
  res_list[[cat]] <- corr
}


# 合并结果
corr_df <- bind_rows(res_list)

# 找出每个特征正相关最强和负相关最强的类别
corr_summary <- corr_df %>%
  group_by(Feature) %>%
  summarise(
    Positive_Category = Category[which.max(Corr)],
    Positive_Corr = max(Corr),
    Negative_Category = Category[which.min(Corr)],
    Negative_Corr = min(Corr)
  ) %>%
  ungroup()

# 查看结果
corr_summary


# 4. 准备绘图数据
# ----------------------------
plot_corr <- corr_summary %>%
  pivot_longer(cols = c(Positive_Corr, Negative_Corr),
               names_to = "Type", values_to = "Corr") %>%
  mutate(Category = ifelse(Type == "Positive_Corr", Positive_Category, Negative_Category),
         Type = ifelse(Type == "Positive_Corr", "Positive", "Negative")) %>%
  left_join(rf_vi %>% dplyr::select(Variable, Importance),
            by = c("Feature" = "Variable"))


plot_corr$Type <- factor(plot_corr$Type, levels = c("Positive", "Negative"))

# ----------------------------
# 5. 绘图
# ----------------------------

p <- ggplot(plot_corr, aes(x = reorder(Feature, Importance), 
                           y = Corr, 
                           fill = Category)) +   # 用 Type 区分正负相关
  geom_bar(stat = "identity", position = "identity", alpha = 0.7) +  
  coord_flip() +
  labs(title = paste0("前", 10, "重要BGC特征与环境相关性 (Top ", 10, " BGCs Correlation with Environment)"),
       x = "BGC Feature",
       y = "Spearman相关系数 (Spearman Correlation)") +
  theme_minimal(base_size = 14) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")
p
topptx(p,"XGB_BGCs_Correlation_with_Environment.pptx", width = 8, height = 5)


##
##
##重要的teature和环境的相关性
library(patchwork)
bgc_class_table <- read.csv("../new_compound_BGC_cluster.csv")
rf_vi <- vip::vi(xgb_fit %>% extract_fit_parsnip())
rf_vi_20 <- rf_vi %>%
  arrange(desc(Importance)) %>%
  slice(1:20) %>%
  left_join(bgc_class_table, by = c("Variable" = "BGC_name"))

# 2. 绘制特征重要性点图（BGC 分类着色）
library(ggplot2)
library(patchwork)  # 用于拼图
library(ggsci)

# 左图：RF重要性（点图）
rf_vi_20 <- rf_vi %>%
  arrange(desc(Importance)) %>%
  slice(1:20) %>%
  left_join(bgc_class_table, by = c("Variable" = "BGC_name"))
rf_vi_20 <- rf_vi_20 %>%
  mutate(BGC_cluster = ifelse(BGC_cluster == "Ribosomally synthesized and post-translationally modified peptide",
                              "RiPPs",
                              ifelse(BGC_cluster == "Nonribosomal peptide",
                                     "NRPs",
                                     BGC_cluster)))

rf_vi_20$BGC_cluster <- factor(rf_vi_20$BGC_cluster,
                               levels = c("Alkaloid","NRPs","Polyketide","RiPPs","Saccharide","Terpene","Other"))



p1 <- ggplot(rf_vi_20, aes(x = reorder(Variable, Importance),
                           y = Importance,
                           color = BGC_cluster)) +
  geom_point(size = 4) +
  coord_flip() +
  labs(x = NULL, y = "Feature Importance") +
  # theme_minimal(base_size = 14) +
  #theme(axis.text.y = element_text(size = 10),
  #     legend.position = "right")  + 
  theme_bw()+
  scale_fill_npg() +  # 应用 'npg' 颜色方案
  scale_color_npg()   # 应用 'npg' 颜色方案
p1
# 右图：相关性（柱状图）
plot_corr$Category <- factor(plot_corr$Category,
                             levels = c("Rhizosphere", "Soil", "water", "human", "Sediment"))

sample_colors <- c(
  "Rhizosphere" = "#008b01",
  "Soil"       = "#8a460f",
  "water"      = "#0000e2",
  "human"      = "#8a0089",
  "Sediment"   = "#ccb994"
)

p2 <- ggplot(plot_corr, aes(x = reorder(Feature, Importance),
                            y = Corr,
                            fill = Category)) +
  geom_bar(stat = "identity", alpha = 0.7,width = 0.8) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1.2) +
  labs(x = NULL, y = "Spearman Correlation") +
  theme_minimal(base_size = 14) +
  
  theme_bw() +
  theme(axis.text.y = element_blank(),   # 去掉右边的Y轴文字
        axis.ticks.y = element_blank(),
        legend.position = "right")  +
  scale_fill_manual(values = sample_colors) 
p2

# 拼接：左右紧挨
final_plot <- p1 + p2 + plot_layout(guides = "collect") & 
  theme(legend.position = "right")


final_plot <- p1 + p2 + plot_layout(guides = "collect") & 
  theme(legend.position = "none",  # Remove the legend
        axis.text.y = element_blank(),  # Remove the y-axis labels
        axis.ticks.y = element_blank(),  # Remove the y-axis ticks
        axis.title.y = element_blank(),  # Remove the y-axis title
        axis.title.x = element_blank(),  # Remove the x-axis title
) 

final_plot 
topptx(final_plot,"BGC_XGB_PLOT.pptx", width = 4, height = 5)

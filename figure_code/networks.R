library(data.table)
library(dplyr)
library(Hmisc)
library(igraph)
library(ggraph)
library(ggplot2)
library(visNetwork)
library(htmlwidgets)
library(tidyr)

# --------------------------
# 1. 用户参数（请根据你的数据修改这些） / User parameters (edit as needed)
# --------------------------
# 1. 输入文件路径与列名（按需修改）
# 1. Input files & column names (modify as needed)
setwd("D:/metagenome/result/2025_data_analysis/Plot")

sample_col  <- "Run"            # 样本列名
category_col<- "Category"            # 分类列名（五类）
cat_to_test <- "Rhizosphere"                # <- 手动修改成你要分析的 Category 名称
outdir_base <- file.path("networks_by_category", cat_to_test)
dir.create(outdir_base, recursive = TRUE, showWarnings = FALSE)

#"Sediment"    "human"       "water"       "Soil"        "Rhizosphere"

# 分析阈值（可调）
prevalence_cut <- 0.05   # 特征至少在 5% 的样本中出现
pseudocount     <- 1e-6  # CLR 伪计数
cor_threshold   <- 0.50  # |rho| 阈值
q_threshold     <- 0.05  # FDR 阈值
min_degree_show <- 2     # 子图保留最小 degree（更严格时可设为2或3）
top_hubs_n      <- 10    # 在图上标注前 N 个 hub

# --------------------------
# 2. 读取数据并对齐样本 / Read data & align samples
# --------------------------
bgc_df <- read.csv("../all_BGC_result.csv",check.names = FALSE)
sp_df <- read.csv("../new_all_OTU_Genus_result.csv",check.names = FALSE)

sp_names <- gsub(".*\\.g__", "", colnames(sp_df)[2:(ncol(sp_df)-7)])
colnames(sp_df)[2:(ncol(sp_df)-7)] <- sp_names

# 选择该类样本
samps <- bgc_df %>% filter(.data[[category_col]] == cat_to_test) %>% pull(!!sym(sample_col))


bgc_num <- bgc_df %>% filter(.data[[sample_col]] %in% samps) 
sp_num  <- sp_df  %>% filter(.data[[sample_col]] %in% samps) 

# 确保行名是 SampleID
rownames(bgc_num) <- bgc_df %>% filter(.data[[sample_col]] %in% samps) %>% pull(!!sym(sample_col))
rownames(sp_num)  <- sp_df  %>% filter(.data[[sample_col]] %in% samps)  %>% pull(!!sym(sample_col))

bgc_num <- bgc_num[,2:(ncol(bgc_num)-7)]
sp_num <- sp_num[,2:(ncol(sp_num)-7)]


# --------------------------
# 4. 过滤低出现率特征（prevalence） / Prevalence filtering
# --------------------------
keep_bgc <- apply(bgc_num, 2, function(x) mean(x > 0, na.rm = TRUE)) >= prevalence_cut
keep_sp  <- apply(sp_num,  2, function(x) mean(x > 0, na.rm = TRUE)) >= prevalence_cut
bgc_filt <- bgc_num[, keep_bgc, drop = FALSE]
sp_filt  <- sp_num[, keep_sp,  drop = FALSE]
message("过滤后特征 — BGC: ", ncol(bgc_filt), "  species: ", ncol(sp_filt))
if(ncol(bgc_filt)==0 || ncol(sp_filt)==0) stop("过滤后没有特征，请放宽 prevalence_cut")

# --------------------------
# 5. CLR 变换（按样本） / CLR transform
# --------------------------
bgc_log <- log(as.matrix(bgc_filt) + pseudocount)
sp_log  <- log(as.matrix(sp_filt)  + pseudocount)
bgc_clr <- sweep(bgc_log, 1, rowMeans(bgc_log), "-")
sp_clr  <- sweep(sp_log,  1, rowMeans(sp_log),  "-")
# 列名保存
bgc_names <- colnames(bgc_clr); sp_names <- colnames(sp_clr)

# --------------------------
# 6. Spearman 相关 + P 值 (Hmisc::rcorr)
# --------------------------
combined <- cbind(bgc_clr, sp_clr)
# 注意：如果特征太多，下面可能较慢；先在一个 category 上调试
rc <- rcorr(as.matrix(combined), type = "spearman")
rmat <- rc$r; pmat <- rc$P

# 只取 BGC x Species 子矩阵
r_bx_s <- rmat[bgc_names, sp_names, drop = FALSE]
p_bx_s <- pmat[bgc_names, sp_names, drop = FALSE]

# --------------------------
# 7. FDR 校正并筛边 / BH correction & edge selection
# --------------------------
p_vec <- as.vector(p_bx_s)
q_vec <- p.adjust(p_vec, method = "BH")
q_mat <- matrix(q_vec, nrow = nrow(p_bx_s), ncol = ncol(p_bx_s),
                dimnames = list(rownames(p_bx_s), colnames(p_bx_s)))

sel_idx <- which((abs(r_bx_s) >= cor_threshold) & (q_mat <= q_threshold), arr.ind = TRUE)
if(nrow(sel_idx) == 0) stop("未筛到边：请放宽阈值")
edges_df <- data.frame(from = rownames(r_bx_s)[sel_idx[,1]],
                       to   = colnames(r_bx_s)[sel_idx[,2]],
                       cor  = r_bx_s[sel_idx],
                       pval = p_bx_s[sel_idx],
                       qval = q_mat[sel_idx],
                       stringsAsFactors = FALSE)
write.csv(edges_df, file = file.path(outdir_base, paste0("edges_cat_", cat_to_test, ".csv")), row.names = FALSE)
message("筛得边数: ", nrow(edges_df))


# 筛选正相关和负相关
edges_pos <- edges_df %>%
  filter(cor >= cor_threshold, qval <= q_threshold)

edges_neg <- edges_df %>%
  filter(cor <= -cor_threshold, qval <= q_threshold)


# --------------------------
# 8. 建图 igraph & 计算中心性 / Build igraph & compute centralities
# --------------------------
nodes <- unique(c(edges_df$from, edges_df$to))
node_type <- ifelse(nodes %in% bgc_names, "BGC", "Species")
node_df <- data.frame(id = nodes, label = nodes, type = node_type, stringsAsFactors = FALSE)

g_pos <- graph_from_data_frame(d = edges_pos, vertices = node_df, directed = FALSE)
g_neg <- graph_from_data_frame(d = edges_neg, vertices = node_df, directed = FALSE)


for(g in list(g_pos, g_neg)){
  V(g)$degree <- degree(g)
  V(g)$betweenness <- betweenness(g)
  V(g)$eigen <- tryCatch(evcent(g)$vector, error = function(e) rep(NA, vcount(g)))
}

# 标注 hub：按 degree 排序，top N 也会被标为 is_hub

hub_pos <- names(sort(degree(g_pos), decreasing = TRUE))[1:10]
hub_neg <- names(sort(degree(g_neg), decreasing = TRUE))[1:10]

V(g_pos)$label[ !(V(g_pos)$name %in% hub_pos) ] <- NA
V(g_neg)$label[ !(V(g_neg)$name %in% hub_neg) ] <- NA

# GraphML (保留属性)

write_graph(g_pos, file.path(outdir_base, paste0("network_gephi_pos_", cat_to_test, ".graphml")), format = "graphml")
message("已导出 Gephi 文件：", outdir_base)
write_graph(g_neg, file.path(outdir_base, paste0("network_gephi_neg_", cat_to_test, ".graphml")), format = "graphml")

library(igraph)

# 定义一个函数：计算网络指标
calc_network_stats <- function(g, name = "network"){
  # 节点数
  n <- vcount(g)
  # 边数
  L <- ecount(g)
  # 平均度
  avgK <- mean(degree(g))
  # giant component 占比（作为网络连通性指标）
  comp <- components(g)
  Con <- max(comp$csize) / vcount(g)
  # 平均聚类系数
  avgCC <- transitivity(g, type = "average")
  # 度相关系数（assortativity）
  rM <- assortativity_degree(g, directed = FALSE)
  
  # 计算中心性
  deg <- degree(g)
  btw <- betweenness(g)
  
  # 定义 keystone：度数前 5% 且介数中心性前 5%
  deg_cut <- quantile(deg, 0.95)
  btw_cut <- quantile(btw, 0.95)
  keystone_nodes <- V(g)$name[deg >= deg_cut & btw >= btw_cut]
  
  res <- list(
    network = name,
    n = n,
    L = L,
    avgK = avgK,
    Con = Con,
    avgCC = avgCC,
    rM = rM,
    n_keystone = length(keystone_nodes),
    keystone_nodes = keystone_nodes
  )
  return(res)
}

# 分别计算正负网络指标
stats_pos <- calc_network_stats(g_pos, "Positive network")
stats_neg <- calc_network_stats(g_neg, "Negative network")

stats_pos <- as.data.frame(stats_pos)
stats_neg <- as.data.frame(stats_neg)
write.csv(stats_pos, file.path(outdir_base, paste0("network_pos_stats_", cat_to_test, ".csv")), row.names = F)
write.csv(stats_neg, file.path(outdir_base, paste0("network_neg_stats", cat_to_test, ".csv")), row.names = F)



library(patchwork)
library(eoffice)
#统计每个网络的各种参数
net_stats <- read.csv("./networks_by_category/network_stats.csv")


metrics_vec <- c("n", "L", "avgK", "Con", "rM","n_keystone")
plot_df <- net_stats %>%
  select(network, Category, all_of(metrics_vec)) %>%
  pivot_longer(cols = all_of(metrics_vec), names_to = "metric", values_to = "value")


plot_df$Category <- factor(plot_df$Category,
                                  levels = c("Rhizosphere", "Soil", "water", "human", "Sediment"))

sample_colors <- c(
  "Rhizosphere" = "#008b01",
  "Soil"       = "#8a460f",
  "water"      = "#0000e2",
  "human"      = "#8a0089",
  "Sediment"   = "#ccb994"
)

plot_fun <- function(data, metric_name, net_type) {
  ggplot(filter(data, metric == metric_name, network == net_type),
         aes(x = Category, y = value, fill = Category)) +
    geom_col(width = 0.6, color = "black", show.legend = FALSE) +
    scale_fill_brewer(palette = "Set2") +
    theme_minimal(base_size = 12) +

    ggtitle(paste0(net_type, " - ", metric_name)) + 
    theme_bw() +
    scale_fill_manual(values = sample_colors) 
}

# ---------- 生成12个子图 ----------
plots <- lapply(metrics_vec, function(m) {
  list(
    plot_fun(plot_df, m, "Positive network"),
    plot_fun(plot_df, m, "Negative network")
  )
})

plots <- do.call(c, plots)

final_plot <- (plots[[1]] | plots[[3]] | plots[[5]] | plots[[7]] | plots[[9]] | plots[[11]]) /
  (plots[[2]] | plots[[4]] | plots[[6]] | plots[[8]] | plots[[10]] | plots[[12]])

final_plot

topptx(final_plot,"network_stats.pptx", width = 16, height = 5)

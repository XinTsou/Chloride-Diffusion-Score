#' =============================================================================
#' CDS v2.1 — Chloride Dynamic Score（改进算法）
#' =============================================================================
#'
#' 概述
#' --------
#' 基于空间转录组数据推断局部细胞外氯离子浓度动态场。
#' 核心理念：细胞膜上氯离子通道/转运体的净表达方向（外流 - 内流）决定了
#' 细胞向胞外排出 Cl⁻ 的倾向，经空间高斯扩散模型平滑后得到局部 Cl⁻ 累积估计。
#'
#' 基因来源
#' --------
#' 基因从此版本的 CDS_gene_catalog.csv（文献 curated）+ CDS_gene_catalog.r
#' （scRNA-seq 表达过滤）获取，通过 custom_genes 参数传入。
#' 这确保：(1) 基因方向基于文献 DOI；(2) 仅保留在 CRC scRNA-seq 中实际
#' 表达的基因；(3) GABRD 等靶点基因不参与 CDS 构建，避免循环论证。
#'
#' 相比 v1 的核心改进
#' --------------------
#' 1. 基因选择：文献目录 + scRNA-seq 表达过滤，方向有 DOI 支撑
#' 2. 平衡性：内流与外流取均值再求差，两侧权重平衡
#' 3. 靶点隔离：GABA/Glycine 受体亚基在 CSV 排除区，不可混入 CDS 输入
#' 4. 鲁棒性：处理缺失基因、常数值、零方差等边界条件
#' 5. 统计检验：置换检验 + Moran's I 空间自相关检测
#' 6. 标准化：同时输出 Z-score 和 rank-based 百分位
#'
#' 核心假设与局限性
#' ------------------
#' 1. mRNA 丰度 ≈ 蛋白丰度 ≈ 通道活性（忽略翻译调控和门控机制）
#' 2. Cl⁻ 在组织中以各向同性方式扩散（忽略组织微结构的方向性）
#' 3. 部分转运体的方向取决于电化学梯度
#' 4. 此分数为描述性推断指标，非验证的生物标志物
#' =============================================================================

library(Seurat)
library(Matrix)

`%||%` <- function(x, y) if (is.null(x)) y else x


# =============================================================================
# 第二部分：核心计算函数
# =============================================================================

#' 计算 Chloride Dynamic Score
#'
#' @param seurat_obj         Seurat 对象（需包含空间坐标，如 Visium 数据）
#' @param custom_genes       基因列表 list(influx = c(...), efflux = c(...))，
#'                           来自 filter_genes_by_scrna() 的输出。必需参数。
#' @param assay              使用哪个 Assay。默认自动检测：SCT > Spatial > DefaultAssay
#' @param layer              表达量数据层，默认 "data"（归一化值）
#' @param sigma              高斯扩散核带宽，单位是**中位最近邻点距的倍数**（默认 1.0，
#'                           Visium 上约合 100 µm）。σ 小于 1/3 点距时核在截断后只剩
#'                           自身权重，CDS 退化为 S_raw，此时发出 warning。
#' @param use_adaptive_sigma 是否启用基于 ECM 的自适应扩散半径（实验性功能，默认 FALSE）
#' @param ecm_genes          ECM 参考基因向量（仅在 use_adaptive_sigma = TRUE 时使用）
#' @param min_valid_genes    每侧最少需要的有效基因数，低于此值发出 warning
#'
#' @return 添加了 CDS_raw, CDS_zscore, CDS_percentile 三列的 Seurat 对象。
#'         参数记录在 seurat_obj@misc$CDS_params 中。
#' @export
calculate_cds <- function(
    seurat_obj,
    custom_genes,                        # REQUIRED — from filter_genes_by_scrna()
    assay              = NULL,
    layer              = "data",
    sigma              = 1.0,
    use_adaptive_sigma = FALSE,
    ecm_genes          = c("COL1A1", "COL1A2", "COL3A1", "FN1"),
    min_valid_genes    = 1L
) {

  # ---- 0. 参数验证 ----
  if (missing(custom_genes) || is.null(custom_genes)) {
    stop("[CDS] 必须提供 custom_genes 参数。请先运行 filter_genes_by_scrna() 获取。")
  }

  influx_genes <- custom_genes$influx %||% character(0)
  efflux_genes <- custom_genes$efflux %||% character(0)

  if (length(influx_genes) == 0 && length(efflux_genes) == 0) {
    stop("[CDS] custom_genes 中 influx 和 efflux 均为空。请检查过滤结果。")
  }

  message(sprintf("[CDS] 使用基因: %d 内流 / %d 外流",
                  length(influx_genes), length(efflux_genes)))

  if (is.null(assay)) {
    available_assays <- names(seurat_obj@assays)
    assay <- if ("SCT" %in% available_assays) {
      "SCT"
    } else if ("Spatial" %in% available_assays) {
      "Spatial"
    } else {
      DefaultAssay(seurat_obj)
    }
    message(sprintf("[CDS] 自动选择 assay: %s", assay))
  }

  if (!assay %in% names(seurat_obj@assays)) {
    stop(sprintf("[CDS] Assay '%s' 不存在。可用 assay: %s",
                 assay, paste(names(seurat_obj@assays), collapse = ", ")))
  }

  # 兼容 Seurat v4 (SCTAssay 无 layers) 和 v5
  if ("layers" %in% slotNames(seurat_obj@assays[[assay]])) {
    seurat_layers <- names(seurat_obj@assays[[assay]]@layers)
    if (!layer %in% seurat_layers) {
      stop(sprintf("[CDS] Layer '%s' 不存在于 assay '%s'。可用 layer: %s",
                   layer, assay, paste(seurat_layers, collapse = ", ")))
    }
  }

  if (sigma <= 0) stop("[CDS] sigma 必须 > 0")

  # ---- 1. 获取基因列表 ----
  # (gene validation performed above in parameter validation)

  # ---- 2. 匹配表达矩阵 ----
  exp_data      <- GetAssayData(seurat_obj, assay = assay, layer = layer)
  valid_influx  <- intersect(influx_genes, rownames(exp_data))
  valid_efflux  <- intersect(efflux_genes, rownames(exp_data))

  message(sprintf("[CDS] 数据中匹配到: %d/%d 内流基因, %d/%d 外流基因",
                  length(valid_influx), length(influx_genes),
                  length(valid_efflux), length(efflux_genes)))

  # 列出未匹配的基因
  missing_influx <- setdiff(influx_genes, rownames(exp_data))
  missing_efflux <- setdiff(efflux_genes, rownames(exp_data))
  if (length(missing_influx) > 0) {
    message(sprintf("[CDS] 未匹配的内流基因: %s", paste(missing_influx, collapse = ", ")))
  }
  if (length(missing_efflux) > 0) {
    message(sprintf("[CDS] 未匹配的外流基因: %s", paste(missing_efflux, collapse = ", ")))
  }

  if (length(valid_influx) == 0 && length(valid_efflux) == 0) {
    stop("[CDS] 数据中未找到任何目标基因。请检查基因符号或数据来源。")
  }
  if (length(valid_influx) < min_valid_genes) {
    warning(sprintf("[CDS] 有效内流基因仅 %d 个，内流端估计不可靠。",
                    length(valid_influx)))
  }
  if (length(valid_efflux) < min_valid_genes) {
    warning(sprintf("[CDS] 有效外流基因仅 %d 个，外流端估计不可靠。",
                    length(valid_efflux)))
  }

  # ---- 3. 计算净通量势 S_raw = mean(efflux) - mean(influx) ----
  # 两侧各取均值再求差，与基因数量无关

  s_influx <- if (length(valid_influx) > 0) {
    colMeans(exp_data[valid_influx, , drop = FALSE])
  } else {
    rep(0, ncol(exp_data))
  }

  s_efflux <- if (length(valid_efflux) > 0) {
    colMeans(exp_data[valid_efflux, , drop = FALSE])
  } else {
    rep(0, ncol(exp_data))
  }

  s_raw <- s_efflux - s_influx

  # 处理常数值（所有 spot 的 S_raw 相同 → CDS 无意义）
  if (sd(s_raw) < .Machine$double.eps) {
    warning("[CDS] S_raw 方差为零，所有细胞净通量势相同。CDS 无空间变异。")
    seurat_obj$CDS_raw       <- s_raw
    seurat_obj$CDS_zscore    <- rep(0, length(s_raw))
    seurat_obj$CDS_percentile <- rep(0.5, length(s_raw))
    seurat_obj@misc$CDS_params <- list(
      assay = assay, layer = layer, sigma = sigma,
      valid_influx = valid_influx, valid_efflux = valid_efflux,
      n_influx = length(valid_influx), n_efflux = length(valid_efflux),
      constant_s_raw = TRUE, timestamp = Sys.time()
    )
    return(seurat_obj)
  }

  # ---- 4. 空间高斯卷积 ----
  coords <- GetTissueCoordinates(seurat_obj)
  if (ncol(coords) < 2) {
    stop("[CDS] 空间坐标不足（需至少 2 列）。请确认对象为空间转录组数据。")
  }

  dist_mat <- as.matrix(dist(coords[, 1:2]))
  n_spots  <- nrow(dist_mat)

  # sigma 以点距为单位，换算成坐标单位；不同平台/样本的点距不同，这一步保证
  # 同一个 sigma 在所有数据上对应同一个物理半径
  d_nn    <- apply(dist_mat, 1, function(v) min(v[v > 0]))
  pitch_u <- median(d_nn)
  sigma_u <- sigma * pitch_u

  # 核在 3 sigma 处截断，sigma < 1/3 点距时每个 spot 只剩自身权重，CDS ≡ S_raw
  if (sigma_u < pitch_u / 3) {
    warning(sprintf(paste0("[CDS] sigma = %.3f 点距 = %.2f 坐标单位，小于 1/3 点距；",
                           "截断后核退化为单位阵，CDS 将等同于 S_raw。"),
                    sigma, sigma_u))
  }

  # 自适应扩散（实验性）
  if (use_adaptive_sigma) {
    message("[CDS] 启用实验性自适应扩散半径")
    ecm_valid <- intersect(ecm_genes, rownames(exp_data))
    if (length(ecm_valid) == 0) {
      warning("[CDS] 未找到 ECM 参考基因，回退到恒定 sigma。")
      adaptive_sigma_vec <- rep(sigma_u, n_spots)
    } else {
      ecm_signal <- colMeans(exp_data[ecm_valid, , drop = FALSE])
      ecm_range  <- max(ecm_signal) - min(ecm_signal)
      if (ecm_range < .Machine$double.eps) {
        adaptive_sigma_vec <- rep(sigma_u, n_spots)
      } else {
        ecm_norm <- (ecm_signal - min(ecm_signal)) / ecm_range
        adaptive_sigma_vec <- sigma_u * (1 - 0.5 * ecm_norm)
      }
    }
  } else {
    adaptive_sigma_vec <- rep(sigma_u, n_spots)
  }

  # 执行高斯卷积
  message(sprintf("[CDS] 执行空间卷积 (n = %d, sigma = %.2f 点距 = %.1f 坐标单位, 点距 = %.1f)...",
                  n_spots, sigma, sigma_u, pitch_u))
  cds_vector      <- numeric(n_spots)
  trunc_threshold <- exp(-3^2 / 2)   # > 3 sigma 截断，误差 < 1%

  for (i in seq_len(n_spots)) {
    si     <- adaptive_sigma_vec[i]
    w      <- exp(-(dist_mat[i, ]^2) / (2 * si^2))
    w[w < trunc_threshold] <- 0
    cds_vector[i] <- sum(s_raw * w)
  }

  # ---- 5. 标准化 ----
  cds_zscore    <- (cds_vector - mean(cds_vector)) / sd(cds_vector)
  cds_percentile <- rank(cds_vector, ties.method = "average") / length(cds_vector)

  # ---- 6. 写入 Seurat 对象 ----
  seurat_obj$CDS_raw       <- cds_vector
  seurat_obj$CDS_zscore    <- cds_zscore
  seurat_obj$CDS_percentile <- cds_percentile

  seurat_obj@misc$CDS_params <- list(
    assay            = assay,
    layer            = layer,
    sigma            = sigma,          # 点距倍数
    sigma_units      = sigma_u,        # 坐标单位
    pitch_units      = pitch_u,        # 中位最近邻点距（坐标单位）
    adaptive_sigma   = use_adaptive_sigma,
    valid_influx     = valid_influx,
    valid_efflux     = valid_efflux,
    n_influx         = length(valid_influx),
    n_efflux         = length(valid_efflux),
    constant_s_raw   = FALSE,
    timestamp        = Sys.time()
  )

  message("[CDS] 计算完成。CDS_zscore / CDS_percentile / CDS_raw 已写入 meta.data。")
  return(seurat_obj)
}


# =============================================================================
# 第三部分：诊断函数
# =============================================================================

#' CDS 诊断报告
#'
#' 报告基因匹配情况、基因平衡性、CDS 分布、空间自相关性等。
#'
#' @param seurat_obj 已运行 calculate_cds() 的 Seurat 对象
#' @export
cds_diagnose <- function(seurat_obj) {
  params <- seurat_obj@misc$CDS_params
  if (is.null(params)) {
    stop("未找到 CDS 参数记录。请先运行 calculate_cds()。")
  }

  cat("\n")
  cat("========== CDS 诊断报告 ==========\n\n")

  # -- 基因匹配 --
  cat("【基因匹配】\n")
  cat(sprintf("  内流: %d 个 — %s\n",
              params$n_influx,
              paste(params$valid_influx, collapse = ", ")))
  cat(sprintf("  外流: %d 个 — %s\n",
              params$n_efflux,
              paste(params$valid_efflux, collapse = ", ")))

  # -- 基因平衡性 --
  cat("\n【基因平衡性】\n")
  max_n  <- max(params$n_influx, params$n_efflux, 1)
  min_n  <- min(params$n_influx, params$n_efflux)
  ratio  <- min_n / max_n
  cat(sprintf("  平衡比 (min/max): %.2f (%d vs %d)\n",
              ratio, params$n_influx, params$n_efflux))
  if (ratio >= 0.5) {
    cat("  判定: 基本平衡\n")
  } else if (ratio >= 0.25) {
    cat("  判定: 轻度失衡，一侧基因可能主导 CDS，结果解读需谨慎\n")
  } else {
    cat("  判定: 严重失衡，一侧基因几乎完全主导 CDS 信号\n")
  }

  # -- CDS 分布 --
  if ("CDS_zscore" %in% colnames(seurat_obj@meta.data)) {
    cds_val <- seurat_obj$CDS_zscore
    cat("\n【CDS_zscore 分布】\n")
    cat(sprintf("  Mean = %+.3f, SD = %.3f\n", mean(cds_val), sd(cds_val)))
    cat(sprintf("  Range: [%.3f, %.3f]\n", min(cds_val), max(cds_val)))
    qs <- quantile(cds_val, probs = c(0.01, 0.25, 0.50, 0.75, 0.99))
    cat(sprintf("  P1=%.3f, Q1=%.3f, Median=%.3f, Q3=%.3f, P99=%.3f\n",
                qs[1], qs[2], qs[3], qs[4], qs[5]))

    # -- 空间自相关 --
    cat("\n【空间自相关】\n")
    cat("  计算 Moran's I (permutation test)... ")
    mi <- morans_i_test(seurat_obj, variable = "CDS_zscore", n_perm = 199)
    cat(sprintf("done.\n"))
    cat(sprintf("  Moran's I = %+.4f, p = %.4f (n_perm = 199)\n",
                mi$observed_I, mi$p_value))
    if (mi$p_value < 0.05) {
      cat("  判定: CDS 具有显著空间自相关性（高斯平滑的预期效果）\n")
    } else {
      cat("  注意: CDS 空间自相关不显著，sigma 可能过小或信号无空间结构\n")
    }
  } else {
    cat("\n【CDS_zscore】未找到，请确认 calculate_cds() 已成功运行。\n")
  }

  # -- 参数摘要 --
  cat("\n【运行参数】\n")
  cat(sprintf("  assay = %s, layer = %s\n", params$assay, params$layer))
  cat(sprintf("  sigma = %.1f, adaptive_sigma = %s\n",
              params$sigma, params$adaptive_sigma))
  if (isTRUE(params$constant_s_raw)) {
    cat("  警告: S_raw 为常数，CDS 值无意义。\n")
  }

  cat("\n====================================\n\n")
}

# =============================================================================
# 第四部分：空间统计检验
# =============================================================================

#' 计算 Moran's I 并做置换检验
#'
#' 使用逆距离权重矩阵 + 置换法评估空间自相关显著性。
#'
#' @param seurat_obj 包含 CDS 结果的 Seurat 对象
#' @param variable   meta.data 中的变量名，默认 "CDS_zscore"
#' @param n_perm     置换次数（建议 ≥ 199）
#' @param seed       随机种子
#' @return list(observed_I, p_value, null_distribution)
#' @export
morans_i_test <- function(seurat_obj, variable = "CDS_zscore",
                          n_perm = 999, seed = 42) {
  set.seed(seed)

  values <- seurat_obj@meta.data[[variable]]
  if (is.null(values)) {
    stop(sprintf("变量 '%s' 不存在于 meta.data 中。", variable))
  }

  coords   <- GetTissueCoordinates(seurat_obj)
  dist_mat <- as.matrix(dist(coords[, 1:2]))
  weights  <- 1 / (1 + dist_mat)
  diag(weights) <- 0

  obs_I <- .compute_morans_i(values, weights)

  null_I <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    null_I[i] <- .compute_morans_i(sample(values), weights)
  }

  p_value <- (sum(abs(null_I) >= abs(obs_I)) + 1) / (n_perm + 1)

  invisible(list(
    observed_I       = obs_I,
    p_value          = p_value,
    null_distribution = null_I
  ))
}

# Moran's I 核心计算
# I = (n / Σw) · [Σᵢⱼ wᵢⱼ zᵢ zⱼ / Σᵢ zᵢ²]
.compute_morans_i <- function(x, w) {
  n    <- length(x)
  z    <- x - mean(x)
  sum_w <- sum(w)
  if (sum_w < .Machine$double.eps) return(0)

  # 用矩阵乘法替代 outer(z, z)，效率更高
  numerator   <- n * sum(w * (z %*% t(z)))
  denominator <- sum_w * sum(z^2)
  if (denominator < .Machine$double.eps) return(0)

  numerator / denominator
}

#' CDS-靶基因关联检验（置换法）
#'
#' 通过随机打乱 CDS 值的空间对应关系来构建零分布。
#' 注意：此检验未完全控制空间自相关混淆，结果应结合生物学知识解读。
#'
#' @param seurat_obj   包含 CDS 结果的 Seurat 对象
#' @param target_genes 目标基因符号（字符向量）
#' @param assay        用于提取基因表达的 assay，默认自动检测
#' @param layer        数据层，默认 "data"
#' @param method       相关性方法，默认 "spearman"（对离群值更稳健）
#' @param n_perm       置换次数
#' @param seed         随机种子
#' @return data.frame(gene, cor_observed, p_value)
#' @export
test_cds_association <- function(
    seurat_obj,
    target_genes,
    assay  = NULL,
    layer  = "data",
    method = "spearman",
    n_perm = 999,
    seed   = 42
) {
  set.seed(seed)

  if (is.null(assay)) {
    assay <- seurat_obj@misc$CDS_params$assay %||% DefaultAssay(seurat_obj)
  }

  cds_vals <- seurat_obj$CDS_zscore
  if (is.null(cds_vals)) stop("请先运行 calculate_cds()。")

  exp_data <- GetAssayData(seurat_obj, assay = assay, layer = layer)
  n_cells  <- ncol(seurat_obj)

  results <- data.frame(
    gene         = character(0),
    cor_observed = numeric(0),
    p_value      = numeric(0),
    stringsAsFactors = FALSE
  )

  for (gene in target_genes) {
    if (!gene %in% rownames(exp_data)) {
      warning(sprintf("基因 '%s' 不在表达矩阵中，跳过。", gene))
      next
    }

    gene_expr <- exp_data[gene, ]

    # 跳过常数值基因
    if (sd(gene_expr) < .Machine$double.eps) {
      warning(sprintf("基因 '%s' 在所有细胞中表达值相同，跳过。", gene))
      next
    }

    obs_cor <- cor(cds_vals, gene_expr, method = method)

    null_cors <- numeric(n_perm)
    for (i in seq_len(n_perm)) {
      null_cors[i] <- cor(cds_vals[sample(n_cells)], gene_expr, method = method)
    }

    p_val <- (sum(abs(null_cors) >= abs(obs_cor)) + 1) / (n_perm + 1)

    results <- rbind(results, data.frame(
      gene         = gene,
      cor_observed = round(obs_cor, 4),
      p_value      = round(p_val, 4),
      stringsAsFactors = FALSE
    ))
  }

  return(results)
}


# =============================================================================
# 第五部分：可视化
# =============================================================================

#' 绘制 CDS 空间分布图
#'
#' @param seurat_obj 包含 CDS 结果的 Seurat 对象
#' @param mode       "zscore" (默认), "percentile", "raw"
#' @param ...        传递给 SpatialFeaturePlot 的额外参数
#' @return ggplot 对象
#' @export
plot_cds <- function(seurat_obj, mode = c("zscore", "percentile", "raw"), ...) {
  mode <- match.arg(mode)

  feature_map <- c(
    zscore     = "CDS_zscore",
    percentile = "CDS_percentile",
    raw        = "CDS_raw"
  )
  feature <- feature_map[mode]

  if (!feature %in% colnames(seurat_obj@meta.data)) {
    stop("请先运行 calculate_cds()。")
  }

  title_map <- c(
    zscore     = "CDS (Z-score)",
    percentile = "CDS (百分位)",
    raw        = "CDS (原始值)"
  )

  SpatialFeaturePlot(seurat_obj, features = feature, ...) +
    scale_fill_gradientn(
      colors = c("#2166AC", "#F7F7F7", "#B2182B"),
      name   = title_map[mode]
    ) +
    ggtitle(title_map[mode])
}

#' 绘制 CDS 与靶基因的散点图
#'
#' @param seurat_obj 包含 CDS 结果的 Seurat 对象
#' @param target_gene 目标基因
#' @param assay       assay，默认自动检测
#' @param layer       数据层，默认 "data"
#' @param method      相关性方法，默认 "spearman"
#' @return ggplot 对象
#' @export
plot_cds_scatter <- function(seurat_obj, target_gene,
                             assay = NULL, layer = "data",
                             method = "spearman") {
  if (is.null(assay)) {
    assay <- seurat_obj@misc$CDS_params$assay %||% DefaultAssay(seurat_obj)
  }

  df <- FetchData(seurat_obj, vars = c("CDS_zscore", target_gene))

  cor_val  <- cor(df[[1]], df[[2]], method = method)
  cor_text <- sprintf("%s = %.3f", method, cor_val)

  ggplot(df, aes(x = .data[[names(df)[1]]], y = .data[[names(df)[2]]])) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "#B2182B", linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf, label = cor_text,
             hjust = 1.1, vjust = 1.5, size = 4) +
    labs(x = "CDS (Z-score)", y = target_gene) +
    theme_minimal()
}


# =============================================================================
# 第六部分：使用示例（取消注释后运行）
# =============================================================================

# ---- 完整工作流 ----
# source("CDS_gene_catalog.r")
# source("CDS_v2_improved.r")
#
# # Step 1: 从 scRNA-seq 过滤基因
# filtered <- filter_genes_by_scrna(scrna_obj, cell_type_column = "cell_type")
# my_genes <- list(influx = filtered$influx_genes, efflux = filtered$efflux_genes)
#
# # Step 2: 在空间数据上计算 CDS
# spatial_obj <- calculate_cds(spatial_obj, custom_genes = my_genes, sigma = 1.0)
# cds_diagnose(spatial_obj)
# plot_cds(spatial_obj, mode = "zscore")

# ---- 靶基因关联检验（独立于 CDS 构建，无循环论证） ----
# results <- test_cds_association(
#   spatial_obj,
#   target_genes = c("GABRD", "SPARC", "STK39", "OSR1"),
#   n_perm = 999
# )
# print(results)

# ---- sigma 敏感性分析 ----
# sigma_values <- c(1, 5, 10, 20, 50)
# for (s in sigma_values) {
#   obj <- calculate_cds(spatial_obj, custom_genes = my_genes, sigma = s)
#   mi  <- morans_i_test(obj, n_perm = 99)
#   cat(sprintf("Sigma=%.0f: Moran's I=%.4f, p=%.4f\n", s, mi$observed_I, mi$p_value))
# }

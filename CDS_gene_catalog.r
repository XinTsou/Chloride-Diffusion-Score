# =============================================================================
# CDS 基因目录模块 — 文献基因目录 + scRNA-seq 表达过滤
# =============================================================================
#
# 职责
# --------
# 1. 读取 CDS_gene_catalog.csv，提供人类可审查的基因目录
# 2. 在 scRNA-seq 数据中按表达丰度过滤，筛选 CRC 组织中实际活跃的
#    氯离子通道/转运体基因
# 3. 输出可直接传入 calculate_cds(custom_genes = ...) 的基因列表
#
# 过滤逻辑
# --------
# 第 1 层：全局检出率 ≥ detection_rate（默认 10%）
# 第 2 层：在 ≥ min_cell_types 种细胞类型中被检测到（每类检出率 > 5%）
#
# 不做跨细胞类型表达变异度过滤——Visium spot 含多个细胞，
# 均匀表达的氯通道基因是 CDS 空间场的稳定基线。
# =============================================================================

library(Seurat)
library(Matrix)


# =============================================================================
# 第一部分：基因目录展示
# =============================================================================

#' 展示基因目录（从 CSV 读取并格式化打印）
#'
#' @param csv_path CDS_gene_catalog.csv 的路径
#' @export
show_gene_catalog <- function(csv_path = "CDS_gene_catalog.csv") {
  if (!file.exists(csv_path)) {
    stop(sprintf("基因目录文件不存在: %s", csv_path))
  }

  catalog <- read.csv(csv_path, stringsAsFactors = FALSE)

  cat("\n========== CDS 基因目录 ==========\n\n")

  for (cat in c("include", "controversial", "excluded")) {
    sub <- catalog[catalog$category == cat, , drop = FALSE]
    if (nrow(sub) == 0) next

    cat_label <- c(
      include      = "纳入 (include)",
      controversial = "争议区 (controversial)",
      excluded     = "排除 (excluded)"
    )[cat]

    cat(sprintf("--- %s (%d 基因) ---\n", cat_label, nrow(sub)))
    for (i in seq_len(nrow(sub))) {
      r <- sub[i, ]
      cat(sprintf(
        "  %-10s %-20s %-14s %-12s [%s]\n",
        r$gene, r$protein, r$direction, r$confidence, r$family
      ))
      if (!is.na(r$note) && nchar(r$note) > 0) {
        cat(sprintf("    %s\n", r$note))
      }
      if (!is.na(r$doi_direction) && nchar(r$doi_direction) > 0) {
        cat(sprintf("    DOI: https://doi.org/%s\n", r$doi_direction))
      }
    }
    cat("\n")
  }

  cat("====================================\n\n")
}


# =============================================================================
# 第二部分：scRNA-seq 表达过滤
# =============================================================================

#' 在 scRNA-seq 数据中过滤氯离子通道基因
#'
#' 读取 CDS_gene_catalog.csv 中 category="include" 的基因（可选争议区），
#' 在 scRNA-seq 数据中执行两层过滤：
#'   第 1 层：全局检出率 ≥ detection_rate
#'   第 2 层：覆盖 ≥ min_cell_types 种细胞类型（每类检出率 > 5%）
#'
#' @param scrna_obj             整合后的 CRC scRNA-seq Seurat 对象
#' @param catalog_file          基因目录 CSV 路径
#' @param assay                 assay 名称
#' @param layer                 表达量数据层
#' @param detection_rate        第 1 层：最低全局检出率（默认 0.10）
#' @param min_cell_types        第 2 层：最少覆盖细胞类型数（默认 3）
#' @param cell_type_column      meta.data 中的细胞类型注释列名
#' @param include_controversial 是否纳入争议区基因（默认 FALSE）
#'
#' @return 列表，包含：
#'   - influx_genes: 通过过滤的内流基因
#'   - efflux_genes: 通过过滤的外流基因
#'   - influx_report: 内流候选全量表
#'   - efflux_report: 外流候选全量表
#'   - parameters: 运行参数记录
#' @export
filter_genes_by_scrna <- function(
    scrna_obj,
    catalog_file          = "CDS_gene_catalog.csv",
    assay                 = "RNA",
    layer                 = "data",
    detection_rate        = 0.10,
    min_cell_types        = 3,
    cell_type_column      = "cell_type",
    include_controversial = FALSE
) {

  # ---- 0. 输入验证 ----
  if (!file.exists(catalog_file)) {
    stop(sprintf("[filter] 基因目录文件不存在: %s", catalog_file))
  }
  if (!cell_type_column %in% colnames(scrna_obj@meta.data)) {
    stop(sprintf("[filter] 列 '%s' 不在 meta.data 中。可用列: %s",
                 cell_type_column,
                 paste(colnames(scrna_obj@meta.data), collapse = ", ")))
  }
  if (!assay %in% names(scrna_obj@assays)) {
    stop(sprintf("[filter] Assay '%s' 不存在。可用: %s",
                 assay, paste(names(scrna_obj@assays), collapse = ", ")))
  }

  message("===== CDS 基因过滤 =====")

  # ---- 1. 读取基因目录 ----
  catalog <- read.csv(catalog_file, stringsAsFactors = FALSE)

  if (include_controversial) {
    catalog <- catalog[catalog$category %in% c("include", "controversial"), , drop = FALSE]
    message("[filter] 纳入争议区基因")
  } else {
    catalog <- catalog[catalog$category == "include", , drop = FALSE]
  }

  # 按方向分类（LRRC8 bidirectional_efflux 和 SLC26A9 efflux_leaning 计入外流）
  influx_candidates <- catalog[catalog$direction == "influx", , drop = FALSE]
  efflux_candidates <- catalog[
    catalog$direction %in% c("efflux", "bidirectional_efflux", "efflux_leaning"),
    , drop = FALSE
  ]

  message(sprintf("[filter] 候选基因: %d 内流 / %d 外流",
                  nrow(influx_candidates), nrow(efflux_candidates)))

  # ---- 2. 获取表达矩阵和细胞类型 ----
  exp_data   <- LayerData(scrna_obj, assay = assay, layer = layer)
  cell_types <- scrna_obj@meta.data[[cell_type_column]]
  unique_ct  <- unique(cell_types)
  n_cells    <- ncol(scrna_obj)

  message(sprintf("[filter] 表达矩阵: %d 基因 x %d 细胞, %d 种细胞类型",
                  nrow(exp_data), n_cells, length(unique_ct)))

  # ---- 3. 计算每个基因的检出率和细胞类型覆盖度 ----
  compute_gene_stats <- function(gene) {
    if (!gene %in% rownames(exp_data)) {
      return(c(detection_rate = 0, n_cell_types = 0))
    }
    expr_vec <- exp_data[gene, ]

    # 全局检出率
    dr <- mean(expr_vec > 0)

    # 每种细胞类型的检出率
    ct_rates <- vapply(unique_ct, function(ct) {
      idx <- which(cell_types == ct)
      if (length(idx) == 0) return(0)
      mean(expr_vec[idx] > 0)
    }, numeric(1))

    c(detection_rate = dr, n_cell_types = sum(ct_rates > 0.05))
  }

  # ---- 4. 对每个候选基因应用过滤 ----
  apply_filters <- function(candidates_df, direction_label) {
    if (nrow(candidates_df) == 0) {
      return(list(
        report = data.frame(
          gene = character(0), detection_rate = numeric(0),
          n_cell_types = integer(0), passed_L1 = logical(0),
          passed_L2 = logical(0), passed = logical(0),
          stringsAsFactors = FALSE
        ),
        passed_genes = character(0)
      ))
    }

    genes  <- candidates_df$gene
    stats  <- t(vapply(genes, compute_gene_stats, numeric(2)))
    report <- data.frame(
      gene           = genes,
      detection_rate = round(stats[, "detection_rate"], 4),
      n_cell_types   = as.integer(stats[, "n_cell_types"]),
      stringsAsFactors = FALSE
    )
    report$passed_L1 <- report$detection_rate >= detection_rate
    report$passed_L2 <- report$n_cell_types >= min_cell_types
    report$passed    <- report$passed_L1 & report$passed_L2

    # 打印排除信息
    excluded_L1 <- report$gene[!report$passed_L1]
    excluded_L2 <- report$gene[report$passed_L1 & !report$passed_L2]

    if (length(excluded_L1) > 0) {
      details <- sprintf("%s(%.1f%%)", excluded_L1,
                         report$detection_rate[!report$passed_L1] * 100)
      message(sprintf("[filter] %s 第1层排除 (检出率<%.0f%%): %s",
                      direction_label, detection_rate * 100,
                      paste(details, collapse = ", ")))
    }
    if (length(excluded_L2) > 0) {
      idx <- which(report$passed_L1 & !report$passed_L2)
      details <- sprintf("%s(%dCT)", report$gene[idx], report$n_cell_types[idx])
      message(sprintf("[filter] %s 第2层排除 (覆盖<%d种细胞类型): %s",
                      direction_label, min_cell_types,
                      paste(details, collapse = ", ")))
    }

    list(
      report       = report,
      passed_genes = report$gene[report$passed]
    )
  }

  influx_result <- apply_filters(influx_candidates, "内流")
  efflux_result <- apply_filters(efflux_candidates, "外流")

  influx_genes <- influx_result$passed_genes
  efflux_genes <- efflux_result$passed_genes

  # ---- 5. 打印摘要 ----
  message(sprintf("\n===== CDS 基因过滤报告 ====="))
  message(sprintf("候选: %d 内流 / %d 外流",
                  nrow(influx_candidates), nrow(efflux_candidates)))
  message(sprintf("通过: %d 内流 / %d 外流",
                  length(influx_genes), length(efflux_genes)))
  if (length(influx_genes) > 0) {
    message(sprintf("内流基因: %s", paste(influx_genes, collapse = ", ")))
  } else {
    message("内流基因: (无)")
  }
  if (length(efflux_genes) > 0) {
    message(sprintf("外流基因: %s", paste(efflux_genes, collapse = ", ")))
  } else {
    message("外流基因: (无)")
  }
  message(sprintf("=============================\n"))

  # ---- 6. 返回 ----
  list(
    influx_genes   = influx_genes,
    efflux_genes   = efflux_genes,
    influx_report  = influx_result$report,
    efflux_report  = efflux_result$report,
    parameters     = list(
      catalog_file          = catalog_file,
      assay                 = assay,
      layer                 = layer,
      detection_rate        = detection_rate,
      min_cell_types        = min_cell_types,
      cell_type_column      = cell_type_column,
      include_controversial = include_controversial,
      n_total_cells         = ncol(scrna_obj),
      n_cell_types          = length(unique_ct),
      timestamp             = Sys.time()
    )
  )
}


# =============================================================================
# 第三部分：使用示例
# =============================================================================

# ---- 示例 1：预览基因目录 ----
# show_gene_catalog("CDS_gene_catalog.csv")

# ---- 示例 2：在 scRNA-seq 数据中过滤基因 ----
# filtered <- filter_genes_by_scrna(
#   scrna_obj,
#   cell_type_column = "cell_type",
#   detection_rate   = 0.10,
#   min_cell_types   = 3
# )
#
# print(filtered$influx_genes)
# print(filtered$efflux_genes)
# View(filtered$influx_report)  # 查看完整过滤表

# ---- 示例 3：纳入争议区基因 ----
# filtered <- filter_genes_by_scrna(
#   scrna_obj,
#   cell_type_column      = "cell_type",
#   include_controversial = TRUE
# )

# ---- 示例 4：传入 CDS 计算 ----
# source("CDS_v2_improved.r")
# filtered <- filter_genes_by_scrna(scrna_obj, cell_type_column = "cell_type")
# my_genes <- list(influx = filtered$influx_genes, efflux = filtered$efflux_genes)
# spatial_obj <- calculate_cds(spatial_obj, custom_genes = my_genes, sigma = 1.0)

# CDS 基因过滤 — 适配多层数据（data.1 + data.2）
# 当 scRNA-seq 对象包含多个 data.* 层时使用
library(Seurat)

filter_genes_multilayer <- function(
    scrna_obj,
    catalog_file          = "CDS_gene_catalog.csv",
    detection_rate        = 0.10,
    min_cell_types        = 3,
    cell_type_column      = "cell_clusters",
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

  message("===== CDS 基因过滤 (多层数据适配) =====")

  # ---- 1. 读取基因目录 ----
  catalog <- read.csv(catalog_file, stringsAsFactors = FALSE)

  if (include_controversial) {
    catalog <- catalog[catalog$category %in% c("include", "controversial"), , drop = FALSE]
    message("[filter] 纳入争议区基因")
  } else {
    catalog <- catalog[catalog$category == "include", , drop = FALSE]
  }

  # 分类
  influx_candidates <- catalog[catalog$direction == "influx", , drop = FALSE]
  efflux_candidates <- catalog[
    catalog$direction %in% c("efflux", "bidirectional_efflux", "efflux_leaning"),
    , drop = FALSE
  ]

  message(sprintf("[filter] 候选基因: %d 内流 / %d 外流",
                  nrow(influx_candidates), nrow(efflux_candidates)))

  # ---- 2. 获取细胞类型 ----
  cell_types <- scrna_obj@meta.data[[cell_type_column]]
  unique_ct  <- unique(cell_types)
  n_cells    <- ncol(scrna_obj)

  message(sprintf("[filter] %d 细胞, %d 种细胞类型", n_cells, length(unique_ct)))

  # ---- 3. 用 FetchData 批量提取候选基因（自动处理多层数据） ----
  all_candidates <- unique(c(catalog$gene))
  message(sprintf("[filter] 用 FetchData 提取 %d 个候选基因的表达...", length(all_candidates)))

  expr_all <- FetchData(scrna_obj, vars = all_candidates, layer = "data")
  message(sprintf("[filter] 成功提取: %d / %d 个基因",
                  sum(colnames(expr_all) %in% all_candidates),
                  length(all_candidates)))

  # ---- 4. 按基因计算统计 ----
  compute_stats <- function(gene) {
    if (!gene %in% colnames(expr_all)) {
      return(c(detection_rate = 0, n_cell_types = 0))
    }
    expr_vec <- expr_all[[gene]]

    # 处理全 NA（基因在数据层中完全缺失）
    if (all(is.na(expr_vec))) {
      return(c(detection_rate = 0, n_cell_types = 0))
    }

    # 全局检出率（排除 NA）
    dr <- mean(expr_vec > 0, na.rm = TRUE)

    # 每种细胞类型的检出率
    ct_rates <- vapply(unique_ct, function(ct) {
      idx <- which(cell_types == ct)
      if (length(idx) == 0) return(0)
      mean(expr_vec[idx] > 0, na.rm = TRUE)
    }, numeric(1))

    c(detection_rate = dr, n_cell_types = sum(ct_rates > 0.05, na.rm = TRUE))
  }

  # ---- 5. 应用过滤 ----
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
    stats  <- t(vapply(genes, compute_stats, numeric(2)))
    report <- data.frame(
      gene           = genes,
      detection_rate = round(stats[, "detection_rate"], 4),
      n_cell_types   = as.integer(stats[, "n_cell_types"]),
      stringsAsFactors = FALSE
    )
    report$passed_L1 <- report$detection_rate >= detection_rate
    report$passed_L2 <- report$n_cell_types >= min_cell_types
    report$passed    <- report$passed_L1 & report$passed_L2
    # 排除 NA (无法评估的基因)
    report$passed[is.na(report$passed)] <- FALSE

    # 打印排除信息
    excluded_L1 <- report$gene[!report$passed_L1 %in% TRUE]
    excluded_L2 <- report$gene[report$passed_L1 %in% TRUE & !report$passed_L2 %in% TRUE]

    if (length(excluded_L1) > 0) {
      details <- sprintf("%s(%.1f%%)", excluded_L1,
                         report$detection_rate[!report$passed_L1 %in% TRUE] * 100)
      message(sprintf("[filter] %s 第1层排除 (检出率<%.0f%%): %s",
                      direction_label, detection_rate * 100,
                      paste(details, collapse = ", ")))
    }
    if (length(excluded_L2) > 0) {
      idx <- which(report$passed_L1 %in% TRUE & !report$passed_L2 %in% TRUE)
      details <- sprintf("%s(%dCT)", report$gene[idx], report$n_cell_types[idx])
      message(sprintf("[filter] %s 第2层排除 (覆盖<%d种细胞类型): %s",
                      direction_label, min_cell_types,
                      paste(details, collapse = ", ")))
    }

    list(
      report       = report,
      passed_genes = report$gene[report$passed %in% TRUE]
    )
  }

  influx_result <- apply_filters(influx_candidates, "内流")
  efflux_result <- apply_filters(efflux_candidates, "外流")

  influx_genes <- influx_result$passed_genes
  efflux_genes <- efflux_result$passed_genes

  # ---- 6. 打印摘要 ----
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

  # ---- 7. 返回 ----
  list(
    influx_genes   = influx_genes,
    efflux_genes   = efflux_genes,
    influx_report  = influx_result$report,
    efflux_report  = efflux_result$report,
    parameters     = list(
      catalog_file          = catalog_file,
      detection_rate        = detection_rate,
      min_cell_types        = min_cell_types,
      cell_type_column      = cell_type_column,
      include_controversial = include_controversial,
      n_total_cells         = n_cells,
      n_cell_types          = length(unique_ct),
      timestamp             = Sys.time()
    )
  )
}

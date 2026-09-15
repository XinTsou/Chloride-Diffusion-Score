# =============================================================================
# CDS demo — full pipeline on a simulated spatial transcriptomics dataset
# =============================================================================
# Self-contained: needs no external data. Run from the repository root:
#
#     Rscript demo/demo_cds.R
#
# The demo builds a Visium-like hexagonal spot grid on which a chloride efflux
# signature and a stromal gene (COL1A1) share a spatial gradient, then runs
# calculate_cds() and the diagnostic / statistical helpers on it.
#
# Output: a console report plus figures in demo/output/
# =============================================================================

suppressPackageStartupMessages({library(Seurat); library(Matrix); library(ggplot2)})

source("CDS_v2_improved.r")

OUT <- "demo/output"
dir.create(OUT, showWarnings = FALSE)

set.seed(2026)

# =============================================================================
# 1. Simulated hexagonal spot grid
# =============================================================================
N_SIDE <- 30          # N_SIDE^2 spots
PITCH  <- 100         # coordinate units between neighbouring spots

grid   <- expand.grid(i = 0:(N_SIDE - 1), j = 0:(N_SIDE - 1))
coords <- data.frame(
  x = PITCH * grid$i,
  y = PITCH * grid$j + ifelse(grid$i %% 2 == 1, PITCH / 2, 0)
)
rownames(coords) <- sprintf("spot_%04d", seq_len(nrow(coords)))
n <- nrow(coords)

cat(sprintf("\nSimulated grid: %d spots, nearest-neighbour pitch = %d coordinate units\n",
            n, PITCH))

# =============================================================================
# 2. Simulated expression
# =============================================================================
# A stromal band occupies the right-hand side of the tissue. Chloride efflux
# genes follow the stroma; the influx gene is expressed in the complementary
# epithelium-like region. Sequencing depth varies across the grid.

u       <- (coords$x - min(coords$x)) / diff(range(coords$x))   # 0 (left) .. 1 (right)
stroma  <- plogis(6 * (u - 0.62))
libsize <- rgamma(n, shape = 4, scale = 2500)
depth   <- libsize / mean(libsize)

draw <- function(pattern, base) rpois(n, depth * base * pattern)

counts <- cbind(
  # signature genes
  SLC12A2 = draw(1.05 - stroma,        12),   # influx
  SLC12A7 = draw(stroma + 0.05,         8),   # efflux
  LRRC8A  = draw(0.5 * stroma + 0.15,   8),   # efflux
  BEST1   = draw(stroma + 0.02,         5),   # efflux
  # stromal marker, co-localised with the efflux arm
  COL1A1  = draw(stroma + 0.05,        12),
  # epithelium-like marker, anti-localised
  EPCAM   = draw(1.05 - stroma,        20),
  # decoy gene with no spatial structure
  ACTB    = draw(1,                    40)
)

# background genes so the object resembles a real assay
bg <- matrix(rpois(n * 200, depth * 6), nrow = n,
             dimnames = list(rownames(coords), sprintf("BG%03d", 1:200)))
counts <- cbind(counts, bg)

obj <- CreateSeuratObject(counts = t(counts), assay = "RNA", project = "CDS_demo")
obj <- NormalizeData(obj, verbose = FALSE)

# attach the coordinates so GetTissueCoordinates() works
obj[["slice1"]] <- CreateFOV(coords[, c("x", "y")], type = "centroids",
                             assay = "RNA", key = "slice1")

cat("\nCells:", ncol(obj), "| layers:", paste(Layers(obj[["RNA"]]), collapse = ", "), "\n")

# =============================================================================
# 3. Run CDS
# =============================================================================
# sigma is expressed in units of the median nearest-neighbour spot pitch, so
# sigma = 1.0 means "smooth over one spot spacing" on any platform.
genes <- list(influx = "SLC12A2", efflux = c("SLC12A7", "BEST1", "LRRC8A"))

obj <- calculate_cds(obj, custom_genes = genes, assay = "RNA", sigma = 1.0)

# =============================================================================
# 4. Diagnostics
# =============================================================================
cds_diagnose(obj)

cat("\n-- Moran's I of CDS (n_perm = 199) --\n")
mi <- morans_i_test(obj, variable = "CDS_zscore", n_perm = 199)
cat(sprintf("   I = %+.4f, p = %.4f\n", mi$observed_I, mi$p_value))

cat("\n-- CDS vs marker genes (spot-shuffle permutation, n_perm = 199) --\n")
print(test_cds_association(obj, target_genes = c("COL1A1", "EPCAM", "ACTB"),
                           assay = "RNA", n_perm = 199))

# =============================================================================
# 5. Effect of sigma
# =============================================================================
# As sigma grows the field is smoothed further and CDS departs from the
# unsmoothed S_raw. Below ~1/3 pitch the 3-sigma truncation leaves only the
# self-weight and CDS collapses onto S_raw.
ed    <- LayerData(obj, assay = "RNA", layer = "data")
s_raw <- colMeans(ed[genes$efflux, , drop = FALSE]) -
         colMeans(ed[genes$influx, , drop = FALSE])

cat("\n-- cor(CDS_raw, S_raw) as a function of sigma --\n")
sweep <- do.call(rbind, lapply(c(0.1, 1.0, 2.0, 3.0), function(s) {
  o <- calculate_cds(obj, custom_genes = genes, assay = "RNA", sigma = s)
  data.frame(sigma_pitch = s,
             sigma_units = o@misc$CDS_params$sigma_units,
             cor_vs_Sraw = cor(o$CDS_raw, s_raw, method = "spearman"))
}))
print(sweep, row.names = FALSE)
write.csv(sweep, file.path(OUT, "sigma_sweep.csv"), row.names = FALSE)

# =============================================================================
# 6. Figures
# =============================================================================
obj <- calculate_cds(obj, custom_genes = genes, assay = "RNA", sigma = 1.0)
df  <- data.frame(x = coords[colnames(obj), "x"], y = coords[colnames(obj), "y"],
                  CDS = obj$CDS_zscore,
                  COL1A1 = as.numeric(ed["COL1A1", colnames(obj)]),
                  SLC12A2 = as.numeric(ed["SLC12A2", colnames(obj)]))

map_plot <- function(dd, value, title) {
  ggplot(dd, aes(x = x, y = y, colour = .data[[value]])) +
    geom_point(size = 1.6, stroke = 0) +
    scale_colour_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"),
                           name = NULL) +
    coord_fixed() +
    labs(title = title) +
    theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

p <- patchwork::wrap_plots(
  map_plot(df, "SLC12A2", "SLC12A2 (influx)"),
  map_plot(df, "COL1A1",  "COL1A1 (stroma)"),
  map_plot(df, "CDS",     "CDS (sigma = 1.0)"),
  nrow = 1)

ggsave(file.path(OUT, "cds_demo_map.png"), p, width = 9, height = 3.4, dpi = 200)
cat(sprintf("\nFigure written to %s\n", file.path(OUT, "cds_demo_map.png")))

# CDS: Chloride Diffusion Score

A spatially adaptive algorithm that computationally infers local extracellular chloride (Cl⁻) potential from spatial transcriptomics data. CDS integrates unbiased Cl⁻ influx/efflux gene signatures with a stromal-constrained Gaussian diffusion kernel to map otherwise invisible ionic gradients in tumour tissue.

## Reference

Zou X, Wang Q, et al. *GABA/GABRD-mediated chloride efflux drives endothelial-to-mesenchymal transition and immune exclusion in colorectal cancer*. (2025)

## Dependencies

- R ≥ 4.0
- [Seurat](https://satijalab.org/seurat/) ≥ 5.0
- Matrix

## Quick start

```r
# 1. Source the algorithm
source("CDS_v2_improved.r")

# 2. Filter genes using scRNA-seq data (optional, or provide your own gene list)
source("CDS_gene_catalog.r")
filtered <- filter_genes_by_scrna(
  scrna_obj,
  catalog_file   = "CDS_gene_catalog.csv",
  cell_type_column = "cell_type",
  detection_rate   = 0.10,
  min_cell_types   = 3
)

# 3. Compute CDS on a 10x Visium Seurat object
my_genes <- list(
  influx = filtered$influx_genes,
  efflux = filtered$efflux_genes
)
spatial_obj <- calculate_cds(spatial_obj, custom_genes = my_genes, sigma = 5)

# 4. Visualise
SpatialFeaturePlot(spatial_obj, features = "CDS")
```

## Files

| File | Description |
|------|-------------|
| `CDS_v2_improved.r` | Core CDS algorithm with permutation testing and Moran's I spatial autocorrelation |
| `CDS_gene_catalog.r` | Gene catalog utilities: display, scRNA-seq expression filtering |
| `CDS_gene_catalog.csv` | Curated catalogue of Cl⁻ channels and transporters with literature DOI support |
| `CDS_filter_genes_multi_layer.r` | Alternative gene filtering adapted for multi-layer scRNA-seq data |

## Algorithm overview

For each Visium spot, a raw chloride potential is defined as:

> *S_raw = mean(efflux genes) − mean(influx genes)*

This pointwise potential is then spatially convolved with an adaptive Gaussian kernel whose bandwidth is modulated by local stromal density (e.g., COL1A1 expression), simulating passive Cl⁻ diffusion across the heterogeneous tissue. The resulting field is Z-score standardised and transformed via a sigmoid activation to yield the final CDS.

## License

MIT

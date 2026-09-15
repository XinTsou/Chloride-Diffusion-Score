# CDS: Chloride Diffusion Score

A spatially adaptive algorithm that infers a local extracellular chloride (Cl⁻) potential
from spatial transcriptomics data. CDS combines a curated Cl⁻ influx/efflux gene signature
with a Gaussian diffusion kernel to map ionic gradients that are not directly observable in
tissue sections.

## Reference

Zou X, Wang Q, et al. *GABA/GABRD-mediated chloride efflux drives endothelial-to-mesenchymal
transition and immune exclusion in colorectal cancer*. (2025)

## Dependencies

- R ≥ 4.0
- [Seurat](https://satijalab.org/seurat/) ≥ 5.0
- Matrix

## Quick start

```r
# 1. Source the algorithm
source("CDS_v2_improved.r")

# 2. Derive the influx/efflux gene list from scRNA-seq (optional — you may
#    supply your own gene list instead)
source("CDS_gene_catalog.r")
filtered <- filter_genes_by_scrna(
  scrna_obj,
  catalog_file     = "CDS_gene_catalog.csv",
  cell_type_column = "cell_type",
  detection_rate   = 0.10,
  min_cell_types   = 3
)
my_genes <- list(influx = filtered$influx_genes, efflux = filtered$efflux_genes)

# 3. Compute CDS on a spatial (e.g. 10x Visium) Seurat object.
#    sigma is in units of the median nearest-neighbour spot pitch (1.0 ≈ 100 µm on Visium).
spatial_obj <- calculate_cds(spatial_obj, custom_genes = my_genes, sigma = 1.0)

# 4. Inspect
cds_diagnose(spatial_obj)
plot_cds(spatial_obj, mode = "zscore")
```

## Demo

`demo/demo_cds.R` runs the whole pipeline — simulation, CDS, diagnostics, statistics and
plots — on a simulated Visium-like dataset. It needs no external data:

```sh
Rscript demo/demo_cds.R
```

It writes a console report, a σ sweep table and a figure to `demo/output/`.
(The demo uses `CreateFOV()` and therefore needs Seurat ≥ 5.1.)

## Algorithm

For each spot, a raw chloride potential is defined as

> *S_raw = mean(efflux genes) − mean(influx genes)*

The two arms are averaged before subtraction, so the score does not depend on how many
genes were retained on either side.

This pointwise potential is then smoothed with a Gaussian kernel of bandwidth σ, so that
each spot accumulates the net flux of its neighbourhood — a discrete approximation to
passive Cl⁻ diffusion over a length scale σ. The kernel is truncated at 3σ.

Three values are written to `meta.data`:

| Column | Contents |
|---|---|
| `CDS_raw` | convolved field, before standardisation |
| `CDS_zscore` | `CDS_raw` standardised to zero mean and unit variance |
| `CDS_percentile` | rank-based percentile of `CDS_raw` within the sample |

### The `sigma` parameter

`sigma` is measured in **units of the median nearest-neighbour spot pitch**, not in
platform coordinates. The default `sigma = 1.0` therefore means "smooth over one spot
spacing" — about 100 µm on 10x Visium — whatever coordinate system the platform reports,
so the same value denotes the same physical radius across samples and platforms.

Because the distance matrix has a zero diagonal, the kernel always keeps a self-weight of
1 at each spot. If `sigma` falls below one third of the spot pitch, the 3σ truncation
removes every neighbour and the kernel collapses to the identity matrix: CDS then equals
`S_raw` and no smoothing takes place. `calculate_cds()` warns when this is the case.

Setting `use_adaptive_sigma = TRUE` (off by default) additionally modulates the bandwidth
by local ECM gene expression, using `ecm_genes`.

## API

| Function | Purpose |
|---|---|
| `calculate_cds()` | compute CDS; writes `CDS_raw` / `CDS_zscore` / `CDS_percentile` to `meta.data` and records the run parameters in `seurat_obj@misc$CDS_params` |
| `cds_diagnose()` | console report: gene matching, arm balance, distribution, Moran's I |
| `morans_i_test()` | Moran's I with an inverse-distance weight matrix and a permutation p-value |
| `test_cds_association()` | correlation of CDS with target genes, with a spot-shuffle permutation p-value |
| `plot_cds()` | spatial map of CDS (`mode = "zscore"`, `"percentile"` or `"raw"`) |
| `plot_cds_scatter()` | CDS against a target gene, with a fitted trend line |
| `filter_genes_by_scrna()` | derive the influx/efflux gene list from an scRNA-seq object |
| `filter_genes_multilayer()` | the same, for multi-layer scRNA-seq objects |
| `show_gene_catalog()` | print the curated gene catalogue |

## Gene catalogue

`CDS_gene_catalog.csv` lists curated chloride channels and transporters together with their
direction of transport (influx / efflux), subcellular localisation and supporting literature
(DOI). `filter_genes_by_scrna()` keeps the entries that are detected above a detection-rate
threshold in at least `min_cell_types` cell types, and returns the surviving influx and
efflux lists.

## Interpretation

CDS is a **descriptive** score. It summarises where a tissue's transcriptional configuration
favours net Cl⁻ efflux, under the assumptions that mRNA abundance tracks transporter
activity, that diffusion is isotropic, and that the direction of transport is fixed for each
transporter. It is not a measurement of chloride concentration and not a validated
biomarker.

`test_cds_association()` builds its null by shuffling CDS values across spots. Since CDS and
any spatially structured target gene are both autocorrelated fields, that null is
anti-conservative — it destroys the autocorrelation of CDS while leaving that of the target
intact. Read the resulting p-value as a test for spatial structure, not as evidence of
specific co-localisation; the latter requires a spatially constrained null such as a block
or shift permutation.

## Files

| File | Description |
|---|---|
| `CDS_v2_improved.r` | core algorithm, diagnostics, spatial statistics and plotting |
| `CDS_gene_catalog.r` | catalogue display and scRNA-seq expression filtering |
| `CDS_gene_catalog.csv` | curated Cl⁻ channel/transporter catalogue with literature DOIs |
| `CDS_filter_genes_multi_layer.r` | gene filtering adapted to multi-layer scRNA-seq objects |
| `demo/demo_cds.R` | self-contained demo on simulated data |

## License

MIT

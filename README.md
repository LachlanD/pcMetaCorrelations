# pcMetaCorrelations

`pcMetaCorrelations` provides tools to identify metadata columns that are driving principal component embeddings.

## What it does

- Identifies metadata variables that are associated with PCA axes.
- Supports numeric and categorical metadata.
- Uses linear-model regression by default for more interpretable results and easier extension.

## Installation

```r
# install.packages("devtools")
devtools::install_local("../pcMetaCorrelations")
```

## Quick start

```r
library(pcMetaCorrelations)
res <- pc_meta_correlations(pca_matrix, metadata, mode = "lm")
plot_pc_meta_heatmap(res)
plot_top_pc_meta(res, n = 10, value = "effect_size")
```

## How it works

- `mode = "lm"` fits a linear model for each metadata/PC combination.
  - numeric metadata: `lm(PC ~ metadata)`
  - categorical metadata: `lm(PC ~ factor(metadata))`
- `mode = "correlation"` preserves the original correlation-based approach.

## Result columns

- `metadata`: metadata field name
- `pc`: PCA axis name
- `type`: `numeric` or `categorical`
- `statistic`: model coefficient or test statistic
- `p.value`: raw p-value for the association
- `adj.p.value`: adjusted p-value
- `effect_size`: magnitude of the association
- `direction`: direction or top factor level

## Seurat example

```r
library(pcMetaCorrelations)
if (requireNamespace("Seurat", quietly = TRUE)) {
  seurat_obj <- Seurat::pbmc_small
  res <- pc_meta_correlations(seurat_obj, reduction = "pca", mode = "lm")
  plot_pc_meta_heatmap(res)
}
```

## Rich metadata example

A richer metadata example is available from `SeuratData::celegans.embryo`, which contains more than 20 metadata columns.

```r
if (requireNamespace("SeuratData", quietly = TRUE)) {
  SeuratData::InstallData("celegans.embryo")
  data("celegans.embryo", package = "celegans.embryo.SeuratData")
  celegans.embryo <- Seurat::UpdateSeuratObject(celegans.embryo)
  res2 <- pc_meta_correlations(celegans.embryo, reduction = "pca", mode = "lm")
  plot_pc_meta_heatmap(res2, top_n = 15, top_pcs = 8)
  plot_top_pc_meta(res2, n = 10, value = "effect_size")
}
```

## Learn more

Read the full package website at [https://lachland.github.io/pcMetaCorrelations/index.html](https://lachland.github.io/pcMetaCorrelations/index.html).

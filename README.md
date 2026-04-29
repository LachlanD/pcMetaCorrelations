# pcMetaCorrelations

`pcMetaCorrelations` provides tools to identify metadata columns that are correlated with principal component axes.

## Installation

```r
# install.packages("devtools")
devtools::install_local("../pcMetaCorrelations")
```

## Example

```r
library(pcMetaCorrelations)
res <- pc_meta_correlations(pca_matrix, metadata)
plot_pc_meta_heatmap(res)
plot_top_pc_meta(res, n = 10, value = "effect_size")
```

Read the full package website at [https://lachland.github.io/pcMetaCorrelations/index.html](https://lachland.github.io/pcMetaCorrelations/index.html).

## Seurat example

```r
library(pcMetaCorrelations)
if (requireNamespace("Seurat", quietly = TRUE)) {
  seurat_obj <- Seurat::pbmc_small
  res <- pc_meta_correlations(seurat_obj, reduction = "pca")
  plot_pc_meta_heatmap(res)
}
```

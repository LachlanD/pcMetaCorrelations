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
```

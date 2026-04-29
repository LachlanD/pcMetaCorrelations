library(testthat)

test_that("pc_meta_correlations works for numeric and categorical metadata", {
  set.seed(42)
  pc <- matrix(rnorm(50), nrow = 10, ncol = 5)
  colnames(pc) <- paste0("PC", 1:5)
  meta <- data.frame(
    age = 1:10,
    batch = factor(rep(c("A", "B"), each = 5)),
    stringsAsFactors = FALSE
  )

  res <- pc_meta_correlations(pc, metadata = meta, method = "pearson", adjust = "none", min.cells = 5)
  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) >= 10)
  expect_true(all(c("metadata", "pc", "type", "statistic", "p.value", "adj.p.value", "effect_size", "direction") %in% colnames(res)))
  expect_true(all(res$type %in% c("numeric", "categorical")))
  expect_equal(sum(res$metadata == "age"), 5)
  expect_equal(sum(res$metadata == "batch"), 5)
})

test_that("plot helpers run without error", {
  set.seed(42)
  pc <- matrix(rnorm(50), nrow = 10, ncol = 5)
  colnames(pc) <- paste0("PC", 1:5)
  meta <- data.frame(
    age = 1:10,
    batch = factor(rep(c("A", "B"), each = 5)),
    stringsAsFactors = FALSE
  )

  res <- pc_meta_correlations(pc, metadata = meta, method = "pearson", adjust = "none", min.cells = 5)
  p1 <- plot_pc_meta_heatmap(res, value = "statistic", top_n = 2, top_pcs = 3)
  p2 <- plot_top_pc_meta(res, n = 2, value = "effect_size")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

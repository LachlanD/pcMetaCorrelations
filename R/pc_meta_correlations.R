utils::globalVariables(c("pc", "metadata", "type"))

#' Identify associations between principal components and metadata variables
#'
#' `pc_meta_correlations()` computes associations between PCA axes and metadata
#' columns. Numeric metadata are tested by correlation, while categorical
#' metadata are tested by ANOVA and reported with eta-squared effect sizes.
#'
#' @param object A Seurat object, a numeric PCA matrix, or a data frame containing
#'   principal component embeddings.
#' @param metadata A data frame of metadata variables. If `object` is a Seurat
#'   object, this defaults to the object's metadata. Otherwise this argument is
#'   required.
#' @param reduction Name of the PCA reduction to use when `object` is a Seurat
#'   object. Default: `"pca"`.
#' @param meta.cols Character or integer vector of metadata columns to analyse.
#'   Default: all columns.
#' @param method Correlation method for numeric metadata or ranking method for
#'   linear-model mode: `"pearson"` or `"spearman"`.
#' @param mode Analysis mode: `"lm"` for linear regression or `"correlation"`
#'   for simple correlation/ANOVA.
#' @param adjust p-value adjustment method passed to `p.adjust()`.
#'   Default: `"BH"`.
#' @param min.cells Minimum non-missing values required for a metadata column.
#'   Default: `10`.
#' @param verbose If `TRUE`, print progress messages.
#'
#' @return A data frame with one row per metadata variable and PC axis pair.
#'   Columns include `metadata`, `pc`, `type`, `statistic`, `p.value`,
#'   `adj.p.value`, `effect_size`, and `direction`.
#'
#' @import ggplot2
#' @export
pc_meta_correlations <- function(
  object,
  metadata = NULL,
  reduction = "pca",
  meta.cols = NULL,
  method = c("pearson", "spearman"),
  mode = c("lm", "correlation"),
  adjust = c("BH", "bonferroni", "none"),
  min.cells = 10,
  verbose = TRUE
) {
  method <- match.arg(method)
  mode <- match.arg(mode)
  adjust <- match.arg(adjust)

  is_seurat <- inherits(object, "Seurat") || inherits(object, "SeuratObject")

  if (is_seurat) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat must be installed to use a Seurat object as input.",
           call. = FALSE)
    }

    if (is.null(metadata)) {
      metadata <- object@meta.data
      if (!is.data.frame(metadata)) {
        stop("Unable to extract metadata from Seurat object.", call. = FALSE)
      }
    }

    if (!reduction %in% names(object@reductions)) {
      stop(sprintf("Seurat reduction '%s' not found.", reduction), call. = FALSE)
    }

    pc_mat <- Seurat::Embeddings(object[[reduction]])
  } else if (is.matrix(object) || is.data.frame(object)) {
    pc_mat <- as.matrix(object)
    if (is.null(metadata)) {
      stop("'metadata' is required when 'object' is not a Seurat object.",
           call. = FALSE)
    }
  } else {
    stop("'object' must be a Seurat object or a PCA matrix/data.frame.",
         call. = FALSE)
  }

  metadata <- as.data.frame(metadata)

  if (!is.null(meta.cols)) {
    metadata <- metadata[, meta.cols, drop = FALSE]
  }

  if (nrow(metadata) != nrow(pc_mat)) {
    if (!is.null(rownames(metadata)) && !is.null(rownames(pc_mat))) {
      if (!all(rownames(pc_mat) %in% rownames(metadata))) {
        stop("Row names of metadata and PCA matrix do not match.", call. = FALSE)
      }
      metadata <- metadata[rownames(pc_mat), , drop = FALSE]
    } else {
      stop("Number of rows in metadata must match the PCA matrix.", call. = FALSE)
    }
  }

  keep <- vapply(metadata, function(x) {
    x <- x[!is.na(x)]
    length(x) >= min.cells && length(unique(x)) > 1
  }, logical(1))

  if (!any(keep)) {
    stop("No metadata columns passed the minimum-value and variance filters.",
         call. = FALSE)
  }

  metadata <- metadata[, keep, drop = FALSE]
  pcs <- colnames(pc_mat)
  if (is.null(pcs)) {
    pcs <- paste0("PC", seq_len(ncol(pc_mat)))
    colnames(pc_mat) <- pcs
  }

  results <- list()
  row_id <- 1L

  for (meta_name in colnames(metadata)) {
    meta_values <- metadata[[meta_name]]
    var_type <- if (is.numeric(meta_values) && length(unique(meta_values[!is.na(meta_values)])) > 2) {
      "numeric"
    } else {
      "categorical"
    }

    for (pc_name in pcs) {
      pc_values <- pc_mat[, pc_name]
      non_na <- !is.na(pc_values) & !is.na(meta_values)
      pc_sub <- pc_values[non_na]
      meta_sub <- meta_values[non_na]

      if (length(pc_sub) < min.cells) {
        next
      }

      if (var_type == "numeric") {
        if (length(unique(meta_sub)) < 2L) next
        if (mode == "correlation") {
          test <- suppressWarnings(stats::cor.test(pc_sub, meta_sub, method = method))
          statistic <- as.numeric(test$estimate)
          pval <- test$p.value
          effect <- abs(statistic)
          direction <- if (statistic > 0) "positive" else if (statistic < 0) "negative" else "zero"
        } else {
          predictor <- if (method == "spearman") base::rank(meta_sub) else meta_sub
          fit <- stats::lm(pc_sub ~ predictor)
          coef_tab <- summary(fit)$coefficients
          statistic <- as.numeric(coef_tab[2, "Estimate"])
          pval <- as.numeric(coef_tab[2, "Pr(>|t|)"])
          effect <- abs(statistic)
          direction <- if (statistic > 0) "positive" else if (statistic < 0) "negative" else "zero"
        }
      } else {
        meta_fac <- factor(meta_sub)
        if (nlevels(meta_fac) < 2L) next
        fit <- stats::lm(pc_sub ~ meta_fac)
        an <- summary(stats::aov(fit))[[1]]
        ss <- an["meta_fac", "Sum Sq"]
        pval <- an["meta_fac", "Pr(>F)"]
        fstat <- an["meta_fac", "F value"]
        total <- sum(an[, "Sum Sq"])
        effect <- if (total > 0) ss / total else NA_real_
        direction <- levels(meta_fac)[which.max(tapply(pc_sub, meta_fac, mean))][1]
        statistic <- as.numeric(fstat)
      }

      results[[row_id]] <- data.frame(
        metadata = meta_name,
        pc = pc_name,
        type = var_type,
        statistic = statistic,
        p.value = pval,
        effect_size = effect,
        direction = direction,
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }
  }

  results <- do.call(rbind, results)
  results$adj.p.value <- stats::p.adjust(results$p.value, method = adjust)
  results$metadata <- as.character(results$metadata)
  results$pc <- as.character(results$pc)
  results$type <- factor(results$type, levels = c("numeric", "categorical"))
  results
}

#' Plot PC–metadata correlation heatmap
#'
#' @param results Data frame returned by `pc_meta_correlations()`.
#' @param value Value to plot: `"statistic"` or `"effect_size"`.
#' @param top_n Number of top metadata variables to show by median absolute value.
#' @param top_pcs Number of top PCs to show by index.
#' @return A ggplot object.
#' @export
plot_pc_meta_heatmap <- function(results,
                                 value = c("statistic", "effect_size"),
                                 top_n = 20,
                                 top_pcs = 10) {
  value <- match.arg(value)
  if (!value %in% names(results)) {
    stop("Value must be one of 'statistic' or 'effect_size'.", call. = FALSE)
  }

  top_meta <- unique(results$metadata[order(-abs(results[[value]]))])[1:min(top_n, nrow(results))]
  top_pc <- unique(results$pc)[1:min(top_pcs, length(unique(results$pc)))]
  plot_data <- results[results$metadata %in% top_meta & results$pc %in% top_pc, , drop = FALSE]

  plot_data$pc <- factor(plot_data$pc, levels = unique(top_pc))
  plot_data$metadata <- factor(plot_data$metadata, levels = rev(unique(top_meta)))

  ggplot2::ggplot(plot_data, ggplot2::aes(x = pc, y = metadata, fill = .data[[value]])) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    ggplot2::labs(x = "PC", y = "Metadata", fill = value) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Plot top PC–metadata associations
#'
#'
#' @param results Data frame returned by `pc_meta_correlations()`.
#' @param metadata Character vector of metadata columns to plot. If `NULL`, the
#'   top variables by effect size are used.
#' @param n Number of metadata variables to display.
#' @param value Value to rank by: `"effect_size"` or `"statistic"`.
#' @return A ggplot object.
#' @export
plot_top_pc_meta <- function(results,
                             metadata = NULL,
                             n = 10,
                             value = c("effect_size", "statistic")) {
  value <- match.arg(value)
  if (!value %in% names(results)) {
    stop("Value must be one of 'effect_size' or 'statistic'.", call. = FALSE)
  }

  if (is.null(metadata)) {
    top_meta <- unique(results$metadata[order(-abs(results[[value]]))])[1:min(n, length(unique(results$metadata)))]
  } else {
    top_meta <- metadata
  }

  plot_data <- results[results$metadata %in% top_meta, , drop = FALSE]
  plot_data$metadata <- factor(plot_data$metadata, levels = unique(top_meta))

  ggplot2::ggplot(plot_data, ggplot2::aes(x = metadata, y = .data[[value]], fill = type)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~pc, scales = "free_y") +
    ggplot2::labs(x = "Metadata", y = value, fill = "Type") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

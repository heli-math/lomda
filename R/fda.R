#' PCA plus Functional Data Analysis over age
#'
#' Fits a PCA plus FDA workflow for longitudinal omics data. Stage 1 applies
#' PCA to the omics block, as in \code{\link{lomda}}. Stage 2 smooths each
#' selected PC score as a function of age, producing age-specific mean PC
#' trajectories that can be plotted, predicted, and mapped back to important
#' metabolites through the PCA loadings.
#'
#' When \pkg{face} is installed and \code{method = "auto"} or
#' \code{method = "face"}, \code{\link[face]{face.sparse}} is used. Otherwise
#' \code{method = "auto"} falls back to \code{\link[stats]{smooth.spline}} on
#' the observed PC scores.
#'
#' @param x A \code{lomda} object or a data frame formatted as described in
#'   \code{\link{lomda}}.
#' @param n_pc Integer. Number of PCs to smooth when \code{x} is a data frame.
#'   Ignored when \code{x} is a \code{lomda} object.
#' @param age_col Character. Name of the age/time column used as the functional
#'   argument. Defaults to \code{"age"}.
#' @param id_col Character. Subject identifier column. Defaults to \code{"ID"}.
#' @param method Character. One of \code{"auto"}, \code{"face"}, or
#'   \code{"spline"}. Defaults to \code{"auto"}.
#' @param grid_length Integer. Number of age-grid points for fitted curves.
#' @param scale Logical. Scale omics features before PCA when \code{x} is a
#'   data frame. Defaults to \code{FALSE}, matching \code{\link{lomda}}.
#' @param center Logical. Center omics features before PCA when \code{x} is a
#'   data frame. Defaults to \code{TRUE}.
#' @param ... Additional arguments passed to \code{face::face.sparse()} or
#'   \code{stats::smooth.spline()}.
#'
#' @return An object of class \code{"lomda_fda"} containing PCA results,
#'   PC-score data, FDA fits for each PC, the age grid, and the matched call.
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 40, n_visits = 3,
#'                            n_features = 20, seed = 42)
#' fda <- lomda_fda(dat, n_pc = 3, method = "spline")
#' predict(fda, ages = c(40, 50, 60))
#' lomda_fda_important(fda, pc = 1, n_top = 5)
#' plot_fda_trajectory(fda, pc = 1)
#'
#' @export
lomda_fda <- function(x,
                      n_pc = 3,
                      age_col = "age",
                      id_col = "ID",
                      method = c("auto", "face", "spline"),
                      grid_length = 100,
                      scale = FALSE,
                      center = TRUE,
                      ...) {
  cl <- match.call()
  method <- match.arg(method)

  if (inherits(x, "lomda")) {
    pca_res <- list(
      pca = x$fit_pca,
      scores = x$scores,
      var_explained = x$var_explained,
      loadings = x$loadings
    )
    n_pc <- x$n_pc
    omics_cols <- x$omics_cols
  } else if (is.data.frame(x)) {
    pca_res <- lomda_pca(x, n_pc = n_pc, scale = scale, center = center)
    meta_cols <- c(id_col, "visit", age_col)
    omics_cols <- setdiff(names(x), meta_cols)
  } else {
    stop("x must be a 'lomda' object or a data frame.")
  }

  scores_df <- pca_res$scores
  if (!age_col %in% names(scores_df))
    stop(age_col, " not found in score data frame.")
  if (!id_col %in% names(scores_df))
    stop(id_col, " not found in score data frame.")

  age <- scores_df[[age_col]]
  if (!is.numeric(age))
    stop(age_col, " must be numeric for FDA smoothing.")
  if (anyNA(age))
    stop(age_col, " contains missing values.")

  age_grid <- seq(min(age), max(age), length.out = grid_length)
  use_method <- .lomda_choose_fda_method(method)

  pc_names <- paste0("PC", seq_len(n_pc))
  missing <- setdiff(pc_names, names(scores_df))
  if (length(missing) > 0)
    stop("score data are missing columns: ", paste(missing, collapse = ", "))

  fda_fits <- lapply(pc_names, function(pc) {
    .lomda_fit_pc_curve(
      scores_df = scores_df,
      pc = pc,
      age_col = age_col,
      id_col = id_col,
      age_grid = age_grid,
      method = use_method,
      ...
    )
  })
  names(fda_fits) <- pc_names

  structure(
    list(
      fit_pca = pca_res$pca,
      scores = scores_df,
      var_explained = pca_res$var_explained,
      loadings = pca_res$loadings,
      fda_fits = fda_fits,
      method = use_method,
      age_col = age_col,
      id_col = id_col,
      age_grid = age_grid,
      n_pc = n_pc,
      omics_cols = omics_cols,
      call = cl
    ),
    class = "lomda_fda"
  )
}

#' Print method for PCA plus FDA fits
#'
#' @param x A \code{lomda_fda} object.
#' @param ... Ignored.
#' @export
print.lomda_fda <- function(x, ...) {
  cat("====================================================\n")
  cat("  LOMDA-FDA: PCA + age-indexed FDA\n")
  cat("====================================================\n\n")
  cat("Call:\n  "); print(x$call); cat("\n")
  cat("Stage 1 - PCA\n")
  cat("  Features      :", length(x$omics_cols), "\n")
  cat("  PCs smoothed  :", x$n_pc, "\n")
  ve <- round(x$var_explained * 100, 2)
  cat("  Var. explained:", paste0(names(ve), " ", ve, "%", collapse = ", "), "\n\n")
  cat("Stage 2 - FDA over ", x$age_col, "\n", sep = "")
  cat("  Method        :", x$method, "\n")
  cat("  Age range     :", paste(round(range(x$age_grid), 2), collapse = " to "), "\n")
  cat("  Grid points   :", length(x$age_grid), "\n\n")
  cat("Use predict(), lomda_fda_important(), plot_fda_trajectory(), or plot().\n")
  invisible(x)
}

#' Predict age-specific PC trajectories from a PCA plus FDA fit
#'
#' @param object A \code{lomda_fda} object.
#' @param ages Numeric vector of ages where the fitted mean PC curves should
#'   be evaluated. If \code{NULL}, the fitted age grid is used.
#' @param newdata Optional new omics data. When supplied, PC scores are
#'   predicted with the stored PCA model and returned alongside age-matched
#'   FDA mean scores and deviations.
#' @param pcs Integer vector of PCs to return. Defaults to all fitted PCs.
#' @param type Character. \code{"mean"} returns the fitted mean curves;
#'   \code{"deviation"} requires \code{newdata} and returns observed PC scores
#'   minus the age-specific fitted mean.
#' @param ... Ignored.
#'
#' @return A data frame with one row per requested age or new observation.
#'
#' @export
predict.lomda_fda <- function(object,
                              ages = NULL,
                              newdata = NULL,
                              pcs = NULL,
                              type = c("mean", "deviation"),
                              ...) {
  if (!inherits(object, "lomda_fda"))
    stop("object must be a 'lomda_fda' object.")

  type <- match.arg(type)
  pc_names <- .lomda_pc_names(object, pcs)

  if (is.null(newdata)) {
    ages <- ages %||% object$age_grid
    out <- data.frame(age = ages)
    names(out)[1] <- object$age_col
    for (pc in pc_names) {
      out[[pc]] <- .lomda_predict_pc_curve(object$fda_fits[[pc]], ages)
    }
    return(out)
  }

  if (!object$age_col %in% names(newdata))
    stop(object$age_col, " must be present in newdata.")

  omics <- newdata[, object$omics_cols, drop = FALSE]
  pc_scores <- as.data.frame(stats::predict(object$fit_pca, newdata = omics))
  pc_scores <- pc_scores[, pc_names, drop = FALSE]

  ages <- newdata[[object$age_col]]
  out <- newdata[, intersect(c(object$id_col, "visit", object$age_col), names(newdata)),
                 drop = FALSE]
  for (pc in pc_names) {
    mean_name <- paste0(pc, "_mean")
    score_name <- paste0(pc, "_score")
    out[[mean_name]] <- .lomda_predict_pc_curve(object$fda_fits[[pc]], ages)
    out[[score_name]] <- pc_scores[[pc]]
    if (type == "deviation") {
      out[[paste0(pc, "_deviation")]] <- out[[score_name]] - out[[mean_name]]
    }
  }
  out
}

#' Rank metabolites by contribution to an FDA-smoothed PC
#'
#' Maps age-smoothed PC trajectories back to the original omics features using
#' PCA loadings. For a single PC, features are ranked by absolute loading.
#' For multiple PCs, loadings are combined using each PC's variance explained.
#'
#' @param x A \code{lomda_fda} or \code{lomda} object.
#' @param pc Integer vector. PCs used for ranking. Defaults to \code{1}.
#' @param n_top Integer. Number of features to return.
#' @param weighted Logical. If \code{TRUE}, weight each PC by its proportion of
#'   variance explained before combining rankings.
#'
#' @return A data frame containing feature names, loadings, absolute loadings,
#'   optional weighted importance scores, and rank.
#'
#' @examples
#' dat <- simulate_lomda_data(seed = 1)
#' fda <- lomda_fda(dat, n_pc = 3, method = "spline")
#' lomda_fda_important(fda, pc = 1, n_top = 10)
#'
#' @export
lomda_fda_important <- function(x, pc = 1, n_top = 10, weighted = TRUE) {
  if (!inherits(x, "lomda_fda") && !inherits(x, "lomda"))
    stop("x must be a 'lomda_fda' or 'lomda' object.")

  pc_names <- paste0("PC", pc)
  missing <- setdiff(pc_names, colnames(x$loadings))
  if (length(missing) > 0)
    stop("Requested PCs not found in loadings: ", paste(missing, collapse = ", "))

  loadings <- x$loadings[, pc_names, drop = FALSE]
  if (weighted) {
    weights <- x$var_explained[pc_names]
    importance <- as.vector(abs(loadings) %*% weights)
  } else {
    importance <- rowSums(abs(loadings))
  }

  out <- data.frame(
    feature = rownames(loadings),
    importance = importance,
    stringsAsFactors = FALSE
  )
  if (length(pc_names) == 1) {
    out$loading <- loadings[, 1]
    out$abs_loading <- abs(loadings[, 1])
  }

  out <- out[order(-out$importance), , drop = FALSE]
  out$rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  utils::head(out, n_top)
}

#' Plot FDA-smoothed PC trajectory over age
#'
#' @param x A \code{lomda_fda} object.
#' @param pc Integer. PC to plot.
#' @param show_subjects Logical. Show individual subject trajectories.
#' @param color_by Optional score-data column used to colour observed points.
#'   Defaults to \code{NULL}, which draws points in a single colour.
#'
#' @return A \code{ggplot2} object.
#'
#' @export
plot_fda_trajectory <- function(x, pc = 1, show_subjects = TRUE,
                                color_by = NULL) {
  if (!inherits(x, "lomda_fda"))
    stop("x must be a 'lomda_fda' object.")

  pc_name <- paste0("PC", pc)
  if (!pc_name %in% names(x$fda_fits))
    stop(pc_name, " not found in FDA fit.")

  scores_df <- x$scores
  curve <- x$fda_fits[[pc_name]]$mean_curve
  if (!is.null(color_by) && !color_by %in% names(scores_df))
    stop(color_by, " not found in score data frame.")

  p <- ggplot2::ggplot(scores_df, ggplot2::aes(x = .data[[x$age_col]],
                                               y = .data[[pc_name]]))
  if (show_subjects) {
    p <- p + ggplot2::geom_line(
      ggplot2::aes(group = .data[[x$id_col]]),
      color = "grey78", linewidth = 0.35, alpha = 0.65
    )
  }

  point_layer <- if (is.null(color_by)) {
    ggplot2::geom_point(color = "#2166AC", size = 2.2, alpha = 0.72)
  } else {
    ggplot2::geom_point(ggplot2::aes(color = .data[[color_by]]),
                        size = 2.2, alpha = 0.72)
  }

  p +
    point_layer +
    ggplot2::geom_line(data = curve,
                       ggplot2::aes(x = .data$age, y = .data$mean),
                       inherit.aes = FALSE, color = "#B23A48",
                       linewidth = 1.35) +
    ggplot2::labs(
      x = x$age_col,
      y = paste(pc_name, "score"),
      color = color_by,
      title = paste("Age-smoothed trajectory for", pc_name),
      subtitle = paste("FDA method:", x$method)
    ) +
    .lomda_theme()
}

#' Plot important metabolites for an FDA-smoothed PC
#'
#' @param x A \code{lomda_fda} or \code{lomda} object.
#' @param pc Integer. PC to use for ranking.
#' @param n_top Integer. Number of metabolites to show.
#'
#' @return A \code{ggplot2} object.
#'
#' @export
plot_fda_importance <- function(x, pc = 1, n_top = 15) {
  imp <- lomda_fda_important(x, pc = pc, n_top = n_top, weighted = FALSE)
  imp$feature <- factor(imp$feature, levels = rev(imp$feature))

  ggplot2::ggplot(imp, ggplot2::aes(x = .data$feature,
                                    y = .data$importance,
                                    fill = .data$loading)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
                                  high = "#B23A48", midpoint = 0) +
    ggplot2::labs(
      x = "Feature",
      y = "Absolute loading",
      fill = "Loading",
      title = paste("Top metabolites contributing to PC", pc)
    ) +
    .lomda_theme()
}

#' Plot method for PCA plus FDA objects
#'
#' @param x A \code{lomda_fda} object.
#' @param type Character. One of \code{"trajectory"}, \code{"importance"}, or
#'   \code{"all"}.
#' @param pc Integer. PC to plot.
#' @param pause Logical. If \code{TRUE}, wait for Return between plots when
#'   multiple plots are requested.
#' @param ... Additional arguments passed to the specialized plotting function.
#'
#' @return A \code{ggplot2} object or a list of objects, invisibly.
#'
#' @export
plot.lomda_fda <- function(x, type = "all", pc = 1,
                           pause = interactive(), ...) {
  type <- match.arg(type, c("trajectory", "importance", "all"))
  types <- if (type == "all") c("trajectory", "importance") else type

  plots <- lapply(types, function(one_type) {
    p <- switch(one_type,
      trajectory = plot_fda_trajectory(x, pc = pc, ...),
      importance = plot_fda_importance(x, pc = pc, ...)
    )
    print(p)
    if (pause && one_type != tail(types, 1)) {
      readline("Press Return to show the next plot...")
    }
    p
  })

  invisible(if (length(plots) == 1) plots[[1]] else plots)
}

.lomda_choose_fda_method <- function(method) {
  if (method == "auto") {
    if (requireNamespace("face", quietly = TRUE)) "face" else "spline"
  } else {
    method
  }
}

.lomda_fit_pc_curve <- function(scores_df, pc, age_col, id_col, age_grid,
                                method, ...) {
  pc_df <- data.frame(
    y = scores_df[[pc]],
    argvals = scores_df[[age_col]],
    subj = scores_df[[id_col]]
  )
  pc_df <- pc_df[stats::complete.cases(pc_df), , drop = FALSE]

  if (method == "face") {
    if (!requireNamespace("face", quietly = TRUE)) {
      stop("Package 'face' is required for method = 'face'. ",
           "Install it or use method = 'spline'.")
    }
    fit <- face::face.sparse(pc_df, argvals.new = age_grid, ...)
    mean_y <- as.numeric(fit$mu.new)
    if (length(mean_y) != length(age_grid)) {
      stop("face::face.sparse() did not return a mean curve on age_grid.")
    }
    curve <- data.frame(age = age_grid, mean = mean_y)
    return(list(method = "face", fit = fit, mean_curve = curve, pc = pc))
  }

  fit <- stats::smooth.spline(x = pc_df$argvals, y = pc_df$y, ...)
  mean_y <- as.numeric(stats::predict(fit, x = age_grid)$y)
  curve <- data.frame(age = age_grid, mean = mean_y)
  list(method = "spline", fit = fit, mean_curve = curve, pc = pc)
}

.lomda_predict_pc_curve <- function(pc_fit, ages) {
  if (pc_fit$method == "spline") {
    return(as.numeric(stats::predict(pc_fit$fit, x = ages)$y))
  }
  stats::approx(
    x = pc_fit$mean_curve$age,
    y = pc_fit$mean_curve$mean,
    xout = ages,
    rule = 2
  )$y
}

.lomda_pc_names <- function(object, pcs = NULL) {
  pcs <- pcs %||% seq_len(object$n_pc)
  pc_names <- paste0("PC", pcs)
  missing <- setdiff(pc_names, names(object$fda_fits))
  if (length(missing) > 0)
    stop("Requested PCs not found in FDA fit: ", paste(missing, collapse = ", "))
  pc_names
}

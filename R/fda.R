#' PCA plus Functional Data Analysis over age
#'
#' Lower-level PCA plus FDA workflow for longitudinal omics data. Most users
#' can call \code{lomda(..., time = "age")} instead. Stage 1 applies PCA to
#' the omics block, as in \code{\link{lomda}}. Stage 2 smooths each selected
#' PC score as a function of age, producing age-specific mean PC trajectories
#' that can be plotted, predicted, and mapped back to important metabolites
#' through the PCA loadings.
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
#' @param visit_col Character. Visit column retained in score and prediction
#'   output when present. Defaults to \code{"visit"}.
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
#' dat <- simulate_lomda_data(n_subjects = 10, n_visits = 3,
#'                            n_features = 8, seed = 42)
#' fda <- lomda_fda(dat, n_pc = 2, method = "spline")
#' predict(fda, ages = c(40, 50, 60))
#' lomda_fda_important(fda, pc = 1, n_top = 5)
#' plot_fda_trajectory(fda, pc = 1)
#'
#' @export
lomda_fda <- function(x,
                      n_pc = 3,
                      age_col = "age",
                      id_col = "ID",
                      visit_col = "visit",
                      method = c("auto", "face", "spline"),
                      grid_length = 100,
                      scale = FALSE,
                      center = TRUE,
                      ...) {
  cl <- match.call()
  method <- match.arg(method)

  if (inherits(x, "lomda")) {
    age_col <- x$age_col %||% age_col
    id_col <- x$id_col %||% id_col
    visit_col <- x$visit_col %||% visit_col
    pca_res <- list(
      pca = x$fit_pca,
      scores = x$scores,
      var_explained = x$var_explained,
      loadings = x$loadings
    )
    n_pc <- x$n_pc
    omics_cols <- x$omics_cols
  } else if (is.data.frame(x)) {
    pca_res <- lomda_pca(x, n_pc = n_pc, id_col = id_col,
                         visit_col = visit_col, age_col = age_col,
                         scale = scale, center = center)
    meta_cols <- c(id_col, visit_col, age_col)
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
      visit_col = visit_col,
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
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 10, n_features = 8, seed = 1)
#' fit <- lomda(dat, n_pc = 2, time = "age", method_fda = "spline")
#' print(fit)
#'
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

#' Summary method for PCA plus FDA fits
#'
#' @param object A \code{lomda_fda} object.
#' @param ... Ignored.
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 10, n_features = 8, seed = 2)
#' fit <- lomda(dat, n_pc = 2, time = "age", method_fda = "spline")
#' summary(fit)
#'
#' @export
summary.lomda_fda <- function(object, ...) {
  cat("====================================================\n")
  cat("  LOMDA-FDA Summary - age-smoothed PC trajectories\n")
  cat("====================================================\n\n")
  cat("Stage 1 - PCA\n")
  ve <- round(object$var_explained * 100, 2)
  cat("  Var. explained:", paste0(names(ve), " ", ve, "%", collapse = ", "), "\n\n")
  cat("Stage 2 - FDA over ", object$age_col, "\n", sep = "")
  cat("  Method      :", object$method, "\n")
  cat("  Age range   :", paste(round(range(object$age_grid), 2), collapse = " to "), "\n")
  cat("  Grid points :", length(object$age_grid), "\n\n")

  curve_summary <- do.call(rbind, lapply(names(object$fda_fits), function(pc) {
    curve <- object$fda_fits[[pc]]$mean_curve$mean
    data.frame(
      PC = pc,
      min_mean = round(min(curve), 4),
      max_mean = round(max(curve), 4),
      range = round(diff(range(curve)), 4),
      stringsAsFactors = FALSE
    )
  }))
  print(curve_summary)
  invisible(object)
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
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 10, n_features = 8, seed = 3)
#' fit <- lomda(dat, n_pc = 2, time = "age", method_fda = "spline")
#' predict(fit, ages = c(35, 45, 55))
#' predict(fit, newdata = dat[1:3, ], type = "deviation")
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
  out <- newdata[, intersect(c(object$id_col, object$visit_col %||% "visit",
                               object$age_col), names(newdata)),
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

#' Rank metabolites by contribution to FDA age-varying PCs
#'
#' Maps age-smoothed PC trajectories back to the original omics features.
#' Each PC is weighted by its age-trajectory variability
#' \deqn{S_k = G^{-1}\sum_g \{\hat\mu_k(a_g) - \bar\mu_k\}^2,}
#' then feature importance is computed by combining these weights with the
#' absolute PCA loadings across the selected PCs.
#'
#' @param x A \code{lomda_fda} or \code{lomda} object.
#' @param pc Integer vector. PCs used for ranking. Defaults to all fitted PCs.
#' @param n_top Integer. Number of features to return.
#' @param weighted Logical. If \code{TRUE}, weight each PC by its age-curve
#'   variability \eqn{S_k}. If \code{FALSE}, combine absolute loadings without
#'   FDA weights.
#'
#' @return A data frame containing feature names, loadings, absolute loadings,
#'   optional weighted importance scores, and rank.
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 10, n_features = 8, seed = 1)
#' fda <- lomda_fda(dat, n_pc = 2, method = "spline")
#' lomda_fda_important(fda, n_top = 10)
#'
#' @export
lomda_fda_important <- function(x, pc = NULL, n_top = 10, weighted = TRUE) {
  if (!inherits(x, "lomda_fda") && !inherits(x, "lomda"))
    stop("x must be a 'lomda_fda' or 'lomda' object.")
  if (!inherits(x, "lomda_fda"))
    stop("x must be a PCA+FDA fit. Use lomda_lmm_important() for PCA+LMM fits.")

  pc <- pc %||% seq_len(x$n_pc)
  pc_names <- paste0("PC", pc)
  missing <- setdiff(pc_names, colnames(x$loadings))
  if (length(missing) > 0)
    stop("Requested PCs not found in loadings: ", paste(missing, collapse = ", "))

  loadings <- x$loadings[, pc_names, drop = FALSE]
  if (weighted) {
    weights <- .lomda_fda_curve_scores(x)[pc_names]
    importance <- as.vector(abs(loadings) %*% weights)
    signed_importance <- as.vector(loadings %*% weights)
  } else {
    importance <- rowSums(abs(loadings))
    signed_importance <- rowSums(loadings)
  }

  out <- data.frame(
    feature = rownames(loadings),
    importance = importance,
    signed_importance = signed_importance,
    stringsAsFactors = FALSE
  )
  if (weighted) {
    out$fda_weighted <- TRUE
  }
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
#' @param pc Integer vector. PCs to plot. Defaults to \code{1}. The newer
#'   \code{pcs} argument can also be used.
#' @param pcs Integer vector. PCs to plot. If supplied, overrides \code{pc}.
#' @param show_subjects Logical. Show individual subject trajectories.
#' @param n_subjects Integer or \code{NULL}. Maximum number of subjects to
#'   draw. Defaults to \code{NULL}, meaning all subjects. Use a smaller number
#'   for large datasets.
#' @param subject_ids Optional vector of subject IDs to draw. Overrides
#'   \code{n_subjects}.
#' @param color_by Optional score-data column used to colour observed points.
#'   Defaults to \code{NULL}, which draws points in a single colour.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 10, n_features = 8, seed = 4)
#' fit <- lomda(dat, n_pc = 2, time = "age", method_fda = "spline")
#' plot_fda_trajectory(fit, pcs = 1:2, n_subjects = 5)
#'
#' @export
plot_fda_trajectory <- function(x, pc = 1, pcs = NULL, show_subjects = TRUE,
                                n_subjects = NULL, subject_ids = NULL,
                                color_by = NULL) {
  if (!inherits(x, "lomda_fda"))
    stop("x must be a 'lomda_fda' object.")

  pcs <- pcs %||% pc
  pc_names <- .lomda_pc_names(x, pcs)

  scores_df <- x$scores
  if (!is.null(color_by) && !color_by %in% names(scores_df))
    stop(color_by, " not found in score data frame.")

  scores_df <- .lomda_select_subjects(
    scores_df, id_col = x$id_col, n_subjects = n_subjects,
    subject_ids = subject_ids
  )

  keep_cols <- unique(c(x$id_col, x$age_col, color_by, pc_names))
  long <- tidyr::pivot_longer(
    scores_df[, keep_cols, drop = FALSE],
    cols = dplyr::all_of(pc_names),
    names_to = "PC",
    values_to = "score"
  )

  curve <- do.call(rbind, lapply(pc_names, function(pc_name) {
    out <- x$fda_fits[[pc_name]]$mean_curve
    out$PC <- pc_name
    out
  }))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data[[x$age_col]],
                                          y = .data$score))
  if (show_subjects) {
    p <- p + ggplot2::geom_line(
      ggplot2::aes(group = interaction(.data[[x$id_col]], .data$PC)),
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
    ggplot2::facet_wrap(~ .data$PC, scales = "free_y") +
    ggplot2::labs(
      x = x$age_col,
      y = "PC score",
      color = color_by,
      title = "Age-smoothed PC trajectories",
      subtitle = paste("FDA method:", x$method)
    ) +
    .lomda_theme()
}

#' Plot important metabolites for an FDA-smoothed PC
#'
#' @param x A \code{lomda_fda} or \code{lomda} object.
#' @param pc Integer vector. PCs to use for ranking. Defaults to all fitted PCs.
#' @param n_top Integer. Number of metabolites to show.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 10, n_features = 8, seed = 5)
#' fit <- lomda(dat, n_pc = 2, time = "age", method_fda = "spline")
#' plot_fda_importance(fit, n_top = 5)
#'
#' @export
plot_fda_importance <- function(x, pc = NULL, n_top = 15) {
  imp <- lomda_fda_important(x, pc = pc, n_top = n_top, weighted = TRUE)
  imp$feature <- factor(imp$feature, levels = rev(imp$feature))

  ggplot2::ggplot(imp, ggplot2::aes(x = .data$feature,
                                    y = .data$importance)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "Feature",
      y = "FDA-weighted loading importance",
      title = "Top metabolites contributing to age-varying PCs"
    ) +
    .lomda_theme()
}

#' Plot method for PCA plus FDA objects
#'
#' @param x A \code{lomda_fda} object.
#' @param type Character. One of \code{"variance"}, \code{"trajectory"},
#'   \code{"importance"}, or \code{"all"}. Defaults to \code{"all"}.
#' @param pc Integer vector. PCs to plot. Defaults to \code{1}.
#' @param pcs Integer vector. PCs to plot. If supplied, overrides \code{pc}.
#' @param n_subjects Integer or \code{NULL}. Maximum number of subjects to
#'   draw in trajectory plots. Defaults to all subjects.
#' @param pause Logical. If \code{TRUE}, wait for Return between plots when
#'   multiple plots are requested.
#' @param ... Additional arguments passed to the specialized plotting function.
#'
#' @return A \code{ggplot2} object or a list of objects, invisibly.
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 10, n_features = 8, seed = 6)
#' fit <- lomda(dat, n_pc = 2, time = "age", method_fda = "spline")
#' plot(fit, type = "trajectory", pcs = 1:2, n_subjects = 5)
#' plot(fit, type = "variance")
#'
#' @export
plot.lomda_fda <- function(x, type = "all", pc = 1, pcs = NULL,
                           n_subjects = NULL, pause = interactive(), ...) {
  type <- match.arg(type, c("variance", "trajectory", "importance", "all"))
  types <- if (type == "all") c("variance", "trajectory", "importance") else type
  trajectory_pcs <- pcs %||% pc
  importance_pcs <- pcs

  plots <- lapply(types, function(one_type) {
    p <- switch(one_type,
      variance = plot_variance_explained(x),
      trajectory = plot_fda_trajectory(x, pcs = trajectory_pcs,
                                       n_subjects = n_subjects, ...),
      importance = plot_fda_importance(x, pc = importance_pcs, ...)
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

.lomda_select_subjects <- function(scores_df, id_col, n_subjects = NULL,
                                   subject_ids = NULL) {
  if (!is.null(subject_ids)) {
    return(scores_df[scores_df[[id_col]] %in% subject_ids, , drop = FALSE])
  }
  if (is.null(n_subjects)) {
    return(scores_df)
  }
  n_subjects <- max(0L, as.integer(n_subjects))
  keep <- utils::head(unique(scores_df[[id_col]]), n_subjects)
  scores_df[scores_df[[id_col]] %in% keep, , drop = FALSE]
}

.lomda_fda_curve_scores <- function(x) {
  scores <- vapply(x$fda_fits, function(fit) {
    mu <- fit$mean_curve$mean
    mean((mu - mean(mu))^2)
  }, numeric(1))
  scores[names(x$fda_fits)]
}

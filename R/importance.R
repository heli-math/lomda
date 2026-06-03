#' Rank important metabolites from a fitted LOMDA model
#'
#' Computes feature-level importance scores using the Stage 2 model attached to
#' a fitted \code{\link{lomda}} object. For PCA+LMM fits, feature loadings are
#' weighted by the visit-effect slope for each PC. For PCA+FDA fits, feature
#' loadings are weighted by the age-trajectory variability of each PC.
#'
#' @param x A fitted \code{lomda} or \code{lomda_fda} object.
#' @param n_top Integer. Number of features to return.
#' @param pc Integer vector of PCs to use. Defaults to all fitted PCs.
#' @param use_abs Logical. If \code{TRUE}, rank by absolute weighted
#'   contributions. Defaults to \code{TRUE}.
#'
#' @return A data frame of ranked omics features.
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 10, n_features = 8, seed = 1)
#' fit_visit <- lomda(dat, n_pc = 2, time = "visit")
#' lomda_important(fit_visit, n_top = 5)
#' fit_age <- lomda(dat, n_pc = 2, time = "age", method_fda = "spline")
#' lomda_important(fit_age, n_top = 5)
#'
#' @export
lomda_important <- function(x, n_top = 10, pc = NULL, use_abs = TRUE) {
  if (inherits(x, "lomda_fda") || identical(x$stage2, "fda")) {
    return(lomda_fda_important(x, pc = pc, n_top = n_top, weighted = TRUE))
  }
  lomda_lmm_important(x, n_top = n_top, pc = pc, use_abs = use_abs)
}

#' Rank metabolites by LMM-weighted PCA loadings
#'
#' For PCA+LMM fits, metabolite importance is computed across selected PCs as
#' \deqn{I_r = \sum_k |w_{rk}\hat\beta_{1k}|,}
#' where \eqn{w_{rk}} is the PCA loading of feature \eqn{r} on PC \eqn{k}, and
#' \eqn{\hat\beta_{1k}} is the fitted visit-effect slope for PC \eqn{k}.
#'
#' @param x A \code{lomda} object fitted with \code{time = "visit"}.
#' @param n_top Integer. Number of features to return.
#' @param pc Integer vector of PCs to use. Defaults to all fitted PCs.
#' @param use_abs Logical. If \code{TRUE}, rank by absolute weighted
#'   contributions. Defaults to \code{TRUE}.
#'
#' @return A data frame containing feature names, importance scores, signed
#'   scores, ranks, and the PC slopes used as weights.
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 10, n_features = 8, seed = 2)
#' fit <- lomda(dat, n_pc = 2, time = "visit")
#' lomda_lmm_important(fit, n_top = 5)
#'
#' @export
lomda_lmm_important <- function(x, n_top = 10, pc = NULL, use_abs = TRUE) {
  if (!inherits(x, "lomda"))
    stop("x must be a 'lomda' object.")
  if (!is.null(x$stage2) && !identical(x$stage2, "lmm"))
    stop("lomda_lmm_important() is available for time = 'visit' fits only.")

  pc <- pc %||% seq_len(x$n_pc)
  pc_names <- paste0("PC", pc)
  missing <- setdiff(pc_names, colnames(x$loadings))
  if (length(missing) > 0)
    stop("Requested PCs not found in loadings: ", paste(missing, collapse = ", "))

  slopes <- .lomda_lmm_slopes(x)[pc_names]
  loadings <- x$loadings[, pc_names, drop = FALSE]
  signed <- as.vector(loadings %*% slopes)
  importance <- if (use_abs) {
    as.vector(abs(loadings) %*% abs(slopes))
  } else {
    signed
  }

  out <- data.frame(
    feature = rownames(loadings),
    importance = importance,
    signed_importance = signed,
    stringsAsFactors = FALSE
  )
  out <- out[order(-abs(out$importance)), , drop = FALSE]
  out$rank <- seq_len(nrow(out))
  attr(out, "pc_slopes") <- slopes
  rownames(out) <- NULL
  utils::head(out, n_top)
}

.lomda_lmm_slopes <- function(x) {
  time_col <- x$time_col %||% "visit"
  time_row_name <- .lomda_bt(time_col)
  slopes <- vapply(names(x$fit_lmm), function(pc) {
    coef_tbl <- as.data.frame(stats::coef(summary(x$fit_lmm[[pc]])))
    row_name <- if (time_col %in% rownames(coef_tbl)) time_col else time_row_name
    if (!row_name %in% rownames(coef_tbl))
      stop("Time coefficient ", time_col, " not found in LMM for ", pc, ".")
    coef_tbl[row_name, "Estimate"]
  }, numeric(1))
  slopes
}

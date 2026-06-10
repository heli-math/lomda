#' Stage 1: PCA on longitudinal omics data
#'
#' Applies PCA to the omics block of a longitudinal dataset. This is the
#' first stage of the LOMDA pipeline and can be called independently for
#' exploratory purposes.
#'
#' @param data A data frame formatted as described in \code{\link{lomda}}.
#' @param n_pc Integer. Number of principal components to extract.
#' @param id_col Character. Name of the subject identifier column.
#' @param visit_col Character. Name of the visit column.
#' @param age_col Character. Name of the age column.
#' @param scale Logical. Scale features to unit variance? Default \code{TRUE}.
#' @param center Logical. Center features? Default \code{TRUE}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{pca}}{The \code{\link[stats]{prcomp}} object.}
#'   \item{\code{scores}}{Data frame: metadata columns + PC score columns.}
#'   \item{\code{var_explained}}{Proportion of variance per PC.}
#'   \item{\code{loadings}}{Matrix of variable loadings (rotation).}
#' }
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 40, n_features = 30, seed = 1)
#' pca_res <- lomda_pca(dat, n_pc = 5)
#' head(pca_res$scores)
#'
#' @importFrom stats prcomp
#' @export
lomda_pca <- function(data,
                      n_pc       = 3,
                      id_col     = "ID",
                      visit_col  = "visit",
                      age_col    = "age",
                      scale      = TRUE,
                      center     = TRUE) {

  meta_cols  <- c(id_col, visit_col, age_col)
  missing_cols <- setdiff(meta_cols, names(data))
  if (length(missing_cols) > 0)
    stop("data must contain columns: ", paste(missing_cols, collapse = ", "))

  omics_cols <- setdiff(names(data), meta_cols)

  if (length(omics_cols) < 2)
    stop("Need at least 2 omics feature columns.")
  if (n_pc > length(omics_cols))
    stop("n_pc exceeds the number of features.")

  pca_obj <- prcomp(data[, omics_cols, drop = FALSE],
                    center = center, scale. = scale)

  pc_names  <- paste0("PC", seq_len(n_pc))
  pc_scores <- as.data.frame(pca_obj$x[, seq_len(n_pc), drop = FALSE])
  names(pc_scores) <- pc_names

  scores_df     <- cbind(data[, meta_cols, drop = FALSE], pc_scores)
  var_explained <- (pca_obj$sdev^2 / sum(pca_obj$sdev^2))[seq_len(n_pc)]
  names(var_explained) <- pc_names

  list(
    pca           = pca_obj,
    scores        = scores_df,
    var_explained = var_explained,
    loadings      = pca_obj$rotation[, seq_len(n_pc), drop = FALSE]
  )
}


#' Stage 2: LMM on PC scores
#'
#' Fits a Linear Mixed Model with random intercept to each PC score derived
#' from Stage 1. This function can be called independently if you already
#' have PC scores.
#'
#' The model for the \eqn{k}-th PC is:
#' \deqn{t_{ijk} = \beta_{0k} + \beta_{1k} \cdot j +
#'   b_{ik} + e_{ijk}.}
#'
#' @param scores_df A data frame containing at least \code{ID}, the time
#'   variable, and columns named \code{PC1}, \code{PC2}, etc.
#' @param n_pc Integer. Number of PC columns to model.
#' @param time_col Character. Time variable used as the fixed effect in the
#'   LMM. Defaults to \code{"visit"}.
#' @param id_col Character. Subject identifier column used for the random
#'   intercept. Defaults to \code{"ID"}.
#' @param REML Logical. Use REML? Default \code{FALSE} (required for LRT).
#'
#' @return A list with:
#' \describe{
#'   \item{\code{fits}}{Named list of \code{lmerMod} full-model fits.}
#'   \item{\code{null_fits}}{Named list of null-model fits with the time
#'     variable dropped.}
#' }
#'
#' @examples
#' dat <- simulate_lomda_data(seed = 7)
#' pca_res <- lomda_pca(dat, n_pc = 3)
#' lmm_res <- lomda_lmm(pca_res$scores, n_pc = 3)
#' names(lmm_res$fits)
#'
#' @importFrom lmerTest lmer
#' @importFrom stats as.formula
#' @export
lomda_lmm <- function(scores_df,
                      n_pc       = 3,
                      time_col   = "visit",
                      id_col     = "ID",
                      REML       = FALSE) {
  pc_names <- paste0("PC", seq_len(n_pc))
  missing  <- setdiff(pc_names, names(scores_df))
  if (length(missing) > 0)
    stop("scores_df is missing columns: ", paste(missing, collapse = ", "))

  if (!time_col %in% names(scores_df))
    stop(time_col, " not found in score data frame.")
  if (!id_col %in% names(scores_df))
    stop(id_col, " not found in score data frame.")

  time_term <- .lomda_bt(time_col)
  fixed_rhs <- time_term
  null_rhs  <- "1"
  random_rhs <- paste0("(1|", .lomda_bt(id_col), ")")

  fits      <- vector("list", n_pc); names(fits)      <- pc_names
  null_fits <- vector("list", n_pc); names(null_fits) <- pc_names

  for (pc in pc_names) {
    fits[[pc]] <- suppressMessages(
      lmerTest::lmer(as.formula(paste0(pc, " ~ ", fixed_rhs, " + ", random_rhs)),
                     data = scores_df, REML = REML)
    )
    null_fits[[pc]] <- suppressMessages(
      lmerTest::lmer(as.formula(paste0(pc, " ~ ", null_rhs, " + ", random_rhs)),
                     data = scores_df, REML = REML)
    )
  }

  list(fits = fits, null_fits = null_fits)
}

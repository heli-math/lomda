#' Longitudinal Omics Multivariate Dimension-reduction Analysis (LOMDA)
#'
#' @description
#' Fits a two-stage model to longitudinal omics data.
#'
#' \strong{Stage 1}: PCA is applied to the omics measurements (metabolites,
#' proteins, etc.) stacked across all visits. The resulting PC scores
#' will be further analysed in the second stage.
#'
#' \strong{Stage 2}: The PC scores are modelled using one of these two workflows.
#'
#' If \code{time = "visit"}, each PC\eqn{j} is modelled by a Linear Mixed Model
#' (LMM) with subject-specific random intercepts:
#'
#' \deqn{t_{ijk} = \beta_{0k} + \beta_{1k} \cdot j + b_{ik} + e_{ijk}.}
#'
#' If \code{time = "age"}, each PC is smoothed as a function of age using
#' Functional Data Analysis (FDA):
#'
#' \deqn{t_{ijk} = \mu_k(age_{ij}) + u_{ik}(age_{ij}) + e_{ijk}.}
#'
#' Here, \eqn{i=1,...,N} indexes subjects,
#' \eqn{j=1,...,M_i} indexes the observed visits for subject \eqn{i},
#' \eqn{k=1,...,K} indexes PCs,
#' \eqn{age_{ij}} is the age of subject \eqn{i} at visit \eqn{j}.
#'
#' @param data A data frame with one row per subject-visit observation,
#' including meta data such as ID, visit, age.
#' @param n_pc Integer. Number of PCs to extract and model in Stage 2.
#' Defaults to \code{3}.
#' @param time Character. Stage 2 time scale. Use \code{"visit"} for PCA+LMM
#'   or \code{"age"} for PCA+FDA.
#' @param method_fda Character. FDA smoother used only when \code{time = "age"}.
#'   One of \code{"auto"}, \code{"face"}, or \code{"spline"}.
#'   Defaults to \code{"auto"}.
#' @param grid_length Integer. Number of age-grid points for FDA curves when
#'   \code{time = "age"}. Defaults to \code{100}.
#' @param scale Logical. Whether to scale the omics features to unit variance
#'   before PCA. Defaults to \code{FALSE}.
#' @param center Logical. Whether to center the omics features before PCA.
#'   Defaults to \code{TRUE}.
#' @param REML Logical. Whether to use Restricted Maximum Likelihood (REML)
#'   in Stage 2. Defaults to \code{FALSE}.
#'
#' @return An object of class \code{"lomda"} with the following components:
#' \describe{
#'   \item{\code{fit_pca}}{The \code{\link[stats]{prcomp}} object from Stage 1.}
#'   \item{\code{scores}}{Data frame of PC scores merged with metadata
#'     (ID, visit, age).}
#'   \item{\code{var_explained}}{Named numeric vector of variance explained by each PC.}
#'   \item{\code{loadings}}{Matrix of PCA loadings (features x PCs).}
#'   \item{\code{stage2}}{Either \code{"lmm"} or \code{"fda"}.}
#'   \item{\code{fit_lmm}}{For \code{time = "visit"}, named list of
#'     \code{lmerMod} objects, one per PC.}
#'   \item{\code{fit_lmm_null}}{For \code{time = "visit"}, named list of null
#'     \code{lmerMod} objects used for LRT.}
#'   \item{\code{fda_fits}}{For \code{time = "age"}, named list of FDA fits.}
#'   \item{\code{n_pc}}{Number of PCs modelled.}
#'   \item{\code{time}}{Time scale used in Stage 2.}
#'   \item{\code{omics_cols}}{Names of the omics feature columns used in PCA.}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 50, n_visits = 3, n_features = 20,
#'                            seed = 42)
#' fit_visit <- lomda(dat, n_pc = 3, time = "visit")
#' print(fit_visit)
#' summary(fit_visit)
#'
#' fit_age <- lomda(dat, n_pc = 3, time = "age", method_fda = "spline")
#' print(fit_age)
#' summary(fit_age)
#'
#' @seealso \code{\link{lomda_lrt}}, \code{\link{lomda_wald}},
#'   \code{\link{plot.lomda}}, \code{\link{simulate_lomda_data}}
#'
#' @importFrom lmerTest lmer
#' @importFrom stats prcomp
#' @export
lomda <- function(data,
                  n_pc      = 3,
                  time      = c("visit", "age"),
                  method_fda = "auto",
                  grid_length = 100,
                  scale     = FALSE,
                  center    = TRUE,
                  REML      = FALSE) {

  cl <- match.call()
  time <- match.arg(time)
  id_col <- "ID"
  visit_col <- "visit"
  age_col <- "age"

  ## ---- Input checks --------------------------------------------------------
  required_cols <- c(id_col, visit_col, age_col)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop("data must contain columns: ", paste(missing_cols, collapse = ", "))

  meta_cols <- c(id_col, visit_col, age_col)
  omics_cols <- setdiff(names(data), meta_cols)
  pca_data <- data[, c(meta_cols, omics_cols), drop = FALSE]

  if (length(omics_cols) < 2)
    stop("Need at least 2 omics feature columns (columns 4+).")

  if (n_pc > length(omics_cols))
    stop("n_pc (", n_pc, ") cannot exceed the number of omics features (",
         length(omics_cols), ").")

  if (time == "age") {
    method_fda <- if (is.null(method_fda)) "auto" else method_fda
    method_fda <- match.arg(method_fda, c("auto", "face", "spline"))
    message("Stage 1: PCA on ", length(omics_cols), " features...")
    message("Stage 2: FDA over age for each of ", n_pc, " PCs...")
    out <- lomda_fda(
      pca_data,
      n_pc = n_pc,
      age_col = age_col,
      id_col = id_col,
      visit_col = visit_col,
      method = method_fda,
      grid_length = grid_length,
      scale = scale,
      center = center
    )
    out$call <- cl
    out$time <- time
    out$time_col <- age_col
    out$stage2 <- "fda"
    class(out) <- c("lomda_fda", "lomda")
    return(out)
  }

  if (!is.null(method_fda) && !identical(method_fda, "auto")) {
    stop("method_fda is only used when time = 'age'. For time = 'visit', ",
         "LOMDA fits the PCA+LMM workflow and no FDA method is used.")
  }

  ## ---- Stage 1: PCA --------------------------------------------------------
  message("Stage 1: PCA on ", length(omics_cols), " features...")
  fit_pca <- lomda_pca(pca_data, n_pc = n_pc, id_col = id_col,
                       visit_col = visit_col, age_col = age_col,
                       scale = scale, center = center)
  scores_df <- fit_pca$scores

  ## ---- Stage 2: LMM --------------------------------------------------------
  message("Stage 2: LMM for each of ", n_pc, " PCs...")
  fit_lmm <- lomda_lmm(scores_df, n_pc = n_pc, time_col = visit_col,
                       id_col = id_col, REML = REML)

  ## ---- Return --------------------------------------------------------------
  structure(
    list(
      fit_pca       = fit_pca$pca,
      scores        = scores_df,
      var_explained = fit_pca$var_explained,
      loadings      = fit_pca$loadings,
      fit_lmm       = fit_lmm$fits,
      fit_lmm_null  = fit_lmm$null_fits,
      stage2        = "lmm",
      time          = time,
      time_col      = visit_col,
      id_col        = id_col,
      visit_col     = visit_col,
      age_col       = age_col,
      n_pc          = n_pc,
      omics_cols    = omics_cols,
      call          = cl
    ),
    class = "lomda"
  )
}


#' Print method for lomda objects
#'
#' @param x A \code{lomda} object.
#' @param ... Ignored.
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 10, n_features = 8, seed = 1)
#' fit <- lomda(dat, n_pc = 2, time = "visit")
#' print(fit)
#'
#' @export
print.lomda <- function(x, ...) {
  cat("====================================================\n")
  cat("  LOMDA: Longitudinal Omics Multivariate DA\n")
  cat("====================================================\n\n")
  cat("Call:\n  "); print(x$call); cat("\n")
  cat("Stage 1 - PCA\n")
  cat("  Features      :", length(x$omics_cols), "\n")
  cat("  PCs extracted :", x$n_pc, "\n")
  ve <- round(x$var_explained * 100, 2)
  cat("  Var. explained:", paste0(names(ve), " ", ve, "%", collapse = ", "), "\n\n")
  cat("Stage 2 - LMM over visit (random intercept)\n")
  cat("  Formula : PC ~ ", x$time_col %||% "visit",
      " + (1 | ", x$id_col %||% "ID", ")\n", sep = "")
  cat("  N obs.  :", nrow(x$scores), "\n")
  cat("  N subj. :", length(unique(x$scores[[x$id_col %||% "ID"]])), "\n\n")
  cat("Use summary(), lomda_lrt(), lomda_wald(), or plot() for results.\n")
  invisible(x)
}


#' Summary method for lomda objects
#'
#' Prints a table of Stage 2 LMM coefficient estimates (fixed effects) for
#' each PC, including the time-effect slope \eqn{\hat\beta_1}.
#'
#' @param object A \code{lomda} object.
#' @param ... Ignored.
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 10, n_features = 8, seed = 2)
#' fit <- lomda(dat, n_pc = 2, time = "visit")
#' summary(fit)
#'
#' @export
summary.lomda <- function(object, ...) {
  if (!is.null(object$stage2) && !identical(object$stage2, "lmm"))
    stop("summary.lomda() is for PCA+LMM fits. Use summary() on the ",
         "lomda_fda object returned by lomda(..., time = 'age').")

  cat("====================================================\n")
  cat("  LOMDA Summary - Stage 2 Fixed Effects\n")
  cat("====================================================\n\n")

  for (pc in colnames(object$loadings)) {
    cat("---", pc,
        sprintf("(%.1f%% variance)", object$var_explained[pc] * 100), "---\n")
    coef_tbl <- as.data.frame(coef(summary(object$fit_lmm[[pc]])))
    ## round for display
    coef_tbl[, sapply(coef_tbl, is.numeric)] <-
      round(coef_tbl[, sapply(coef_tbl, is.numeric)], 5)
    print(coef_tbl)
    cat("\n")
  }
  invisible(object)
}

#' Longitudinal Omics Multivariate Dimension-reduction Analysis (LOMDA)
#'
#' @description
#' Fits the two-stage \eqn{lomda} model to longitudinal omics data.
#'
#' \strong{Stage 1}: PCA is applied to the omics measurements (metabolites,
#' proteins, etc.) pooled across all visits. The resulting PC scores capture
#' the major axes of variation in the high-dimensional data.
#'
#' \strong{Stage 2}: For each selected PC, a Linear Mixed Model (LMM) with
#' subject-specific random intercepts is fitted:
#'
#' \deqn{t_{ik} = \beta_0 + \beta_1 \cdot k + b_i + g_{ik}}
#'
#' where \eqn{b_i \sim N(0, \Sigma_b)} is the random intercept and
#' \eqn{g_{ik} \sim N(0, \Sigma_g)}.
#'
#' @param data A data frame with columns: \code{ID} (subject identifier),
#'   \code{visit} (visit; integer 1, 2, 3, ...), \code{age} (or other
#'   time-invariant covariates), and one column per omics feature. The omics
#'   columns must start from column 4 onward.
#' @param n_pc Integer. Number of PCs to extract and model in Stage 2.
#'   Defaults to \code{3}.
#' @param covariates Character vector of additional covariate column names to
#'   include in the Stage 2 LMM (e.g., \code{c("age", "visit")}). Defaults to
#'   \code{"visit"}.
#' @param scale Logical. Whether to scale the omics features to unit variance
#'   before PCA. Defaults to \code{TRUE}.
#' @param center Logical. Whether to center the omics features before PCA.
#'   Defaults to \code{TRUE}.
#' @param REML Logical. Whether to use Restricted Maximum Likelihood (REML)
#'   in Stage 2. Defaults to \code{FALSE} (ML, required for LRT).
#'
#' @return An object of class \code{"lomda"} with the following components:
#' \describe{
#'   \item{\code{pca}}{The \code{\link[stats]{prcomp}} object from Stage 1.}
#'   \item{\code{scores}}{Data frame of PC scores merged with metadata
#'     (ID, visit, covariates).}
#'   \item{\code{lmm_fits}}{Named list of \code{lmerMod} objects, one per PC.}
#'   \item{\code{lmm_null_fits}}{Named list of null \code{lmerMod} objects
#'     (time dropped) used for LRT.}
#'   \item{\code{n_pc}}{Number of PCs modelled.}
#'   \item{\code{covariates}}{Covariate names used in Stage 2.}
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{data}}{The input data (metadata + PC scores), stored for
#'     downstream plotting.}
#' }
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 50, n_visits = 3, n_features = 20,
#'                            seed = 42)
#' fit <- lomda(dat, n_pc = 3, covariates = "visit")
#' print(fit)
#' summary(fit)
#'
#' @seealso \code{\link{lomda_lrt}}, \code{\link{lomda_wald}},
#'   \code{\link{plot.lomda}}, \code{\link{simulate_lomda_data}}
#'
#' @importFrom lmerTest lmer
#' @importFrom stats prcomp
#' @export
lomda <- function(data,
                  n_pc      = 3,
                  covariates = NULL,
                  scale     = TRUE,
                  center    = TRUE,
                  REML      = FALSE) {

  cl <- match.call()

  ## ---- Input checks --------------------------------------------------------
  required_cols <- c("ID", "visit")
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop("data must contain columns: ", paste(missing_cols, collapse = ", "))

  if (!all(covariates %in% names(data)) && length(covariates) > 0)
    stop("Some covariates not found in data: ",
         paste(setdiff(covariates, names(data)), collapse = ", "))

  meta_cols   <- c("ID", "visit", "age")
  omics_cols  <- setdiff(names(data), meta_cols)

  if (length(omics_cols) < 2)
    stop("Need at least 2 omics feature columns (columns 4+).")

  if (n_pc > length(omics_cols))
    stop("n_pc (", n_pc, ") cannot exceed the number of omics features (",
         length(omics_cols), ").")

  ## ---- Stage 1: PCA --------------------------------------------------------
  message("Stage 1: PCA on ", length(omics_cols), " features...")
  pca_obj <- prcomp(data[, omics_cols, drop = FALSE],
                    center = center, scale. = scale)

  pc_scores <- as.data.frame(pca_obj$x[, seq_len(n_pc), drop = FALSE])
  pc_names  <- paste0("PC", seq_len(n_pc))
  names(pc_scores) <- pc_names

  scores_df <- cbind(data[, meta_cols, drop = FALSE], pc_scores)

  ## ---- Stage 2: LMM --------------------------------------------------------
  message("Stage 2: LMM for each of ", n_pc, " PCs...")

  cov_str   <- if (length(covariates) > 0) paste(covariates, collapse = " + ")
               # else NULL
  fixed_rhs <- if (!is.null(cov_str)) paste("visit +", cov_str) else "visit"
  null_rhs  <- if (!is.null(cov_str)) cov_str else "1"

  lmm_fits      <- vector("list", n_pc)
  lmm_null_fits <- vector("list", n_pc)
  names(lmm_fits)      <- pc_names
  names(lmm_null_fits) <- pc_names

  for (pc in pc_names) {
    full_formula <- as.formula(
      paste0(pc, " ~ ", fixed_rhs, " + (1 | ID)"))
    null_formula <- as.formula(
      paste0(pc, " ~ ", null_rhs, " + (1 | ID)"))

    lmm_fits[[pc]] <- suppressMessages(
      lmerTest::lmer(full_formula, data = scores_df, REML = REML)
    )
    lmm_null_fits[[pc]] <- suppressMessages(
      lmerTest::lmer(null_formula, data = scores_df, REML = REML)
    )
  }

  ## ---- Variance explained --------------------------------------------------
  var_explained <- (pca_obj$sdev^2 / sum(pca_obj$sdev^2))[seq_len(n_pc)]
  names(var_explained) <- pc_names

  ## ---- Return --------------------------------------------------------------
  structure(
    list(
      pca           = pca_obj,
      scores        = scores_df,
      lmm_fits      = lmm_fits,
      lmm_null_fits = lmm_null_fits,
      n_pc          = n_pc,
      covariates    = covariates,
      var_explained = var_explained,
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
  cat("Stage 2 - LMM (random intercept)\n")
  cat("  Formula : PC ~ ", paste(x$covariates, collapse = " + "),
      "+ (1 | ID)\n")
  cat("  N obs.  :", nrow(x$scores), "\n")
  cat("  N subj. :", length(unique(x$scores$ID)), "\n\n")
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
#' @export
summary.lomda <- function(object, ...) {
  cat("====================================================\n")
  cat("  LOMDA Summary — Stage 2 Fixed Effects\n")
  cat("====================================================\n\n")

  for (pc in names(object$lmm_fits)) {
    cat("---", pc,
        sprintf("(%.1f%% variance)", object$var_explained[pc] * 100), "---\n")
    coef_tbl <- as.data.frame(coef(summary(object$lmm_fits[[pc]])))
    # Round for display
    coef_tbl[, sapply(coef_tbl, is.numeric)] <-
      round(coef_tbl[, sapply(coef_tbl, is.numeric)], 5)
    print(coef_tbl)
    cat("\n")
  }
  invisible(object)
}

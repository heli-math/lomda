#' Likelihood Ratio Test for the time effect in LOMDA
#'
#' Tests whether the time-effect slope (\eqn{\bm{\beta}_1}) is significantly
#' different from zero for each PC, by comparing the full LMM (with time)
#' against a nested null model (without time). Both models are fitted by
#' maximum likelihood (\code{REML = FALSE}).
#'
#' The LRT statistic is:
#' \deqn{LRT = -2(\ell_{\text{null}} - \ell_{\text{full}}) \sim \chi^2(1)}
#'
#' @param x A \code{lomda} object (fitted with \code{REML = FALSE}).
#'
#' @return A data frame with one row per PC and columns:
#' \describe{
#'   \item{\code{PC}}{PC name.}
#'   \item{\code{logLik_full}}{Log-likelihood of the full model.}
#'   \item{\code{logLik_null}}{Log-likelihood of the null model.}
#'   \item{\code{LRT_stat}}{Likelihood ratio test statistic.}
#'   \item{\code{df}}{Degrees of freedom (1).}
#'   \item{\code{p_value}}{P-value from the chi-squared distribution.}
#' }
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 60, seed = 99)
#' fit <- lomda(dat, n_pc = 3, REML = FALSE)
#' lomda_lrt(fit)
#'
#' @importFrom stats logLik pchisq
#' @export
lomda_lrt <- function(x) {
  if (!inherits(x, "lomda"))
    stop("x must be a 'lomda' object.")

  results <- lapply(names(x$lmm_fits), function(pc) {
    ll_full <- as.numeric(logLik(x$lmm_fits[[pc]]))
    ll_null <- as.numeric(logLik(x$lmm_null_fits[[pc]]))
    lrt     <- -2 * (ll_null - ll_full)
    pval    <- pchisq(lrt, df = 1, lower.tail = FALSE)
    data.frame(
      PC          = pc,
      logLik_full = round(ll_full, 4),
      logLik_null = round(ll_null, 4),
      LRT_stat    = round(lrt,     4),
      df          = 1L,
      p_value     = signif(pval, 4),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  class(out) <- c("lomda_lrt", "data.frame")
  out
}


#' Wald Test for the time-effect slope in LOMDA
#'
#' Performs a Wald test on the time-effect slope \eqn{\hat{\bm{\beta}}_1} for each
#' PC. The test statistic is:
#' \deqn{W = \frac{\hat{\bm{\beta}}_1}{\text{SE}(\hat{\bm{\beta}}_1)}}
#' which is compared against a standard normal distribution (z-test) or,
#' when \code{use_t = TRUE}, against a t-distribution with Satterthwaite
#' degrees of freedom (from \pkg{lmerTest}).
#'
#' @param x A \code{lomda} object.
#' @param use_t Logical. If \code{TRUE} (default), use the t-statistic and
#'   p-value from \pkg{lmerTest}'s Satterthwaite approximation. If
#'   \code{FALSE}, use a z-test (large-sample approximation).
#'
#' @return A data frame with one row per PC and columns:
#' \describe{
#'   \item{\code{PC}}{PC name.}
#'   \item{\code{beta1}}{Estimated time-effect slope.}
#'   \item{\code{se}}{Standard error of the estimate.}
#'   \item{\code{statistic}}{Test statistic (t or z).}
#'   \item{\code{df}}{Degrees of freedom (NA for z-test).}
#'   \item{\code{p_value}}{Two-sided p-value.}
#'   \item{\code{ci_lower}}{95\% confidence interval lower bound.}
#'   \item{\code{ci_upper}}{95\% confidence interval upper bound.}
#' }
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 60, seed = 99)
#' fit <- lomda(dat, n_pc = 3)
#' lomda_wald(fit)
#' lomda_wald(fit, use_t = FALSE)
#'
#' @importFrom stats coef vcov qnorm pchisq
#' @export
lomda_wald <- function(x, use_t = TRUE) {
  if (!inherits(x, "lomda"))
    stop("x must be a 'lomda' object.")

  results <- lapply(names(x$lmm_fits), function(pc) {
    fit <- x$lmm_fits[[pc]]

    if (use_t) {
      # lmerTest provides Satterthwaite df + t-stat + p-value in coef(summary())
      cs      <- as.data.frame(coef(summary(fit)))
      time_row <- cs["visit", , drop = FALSE]
      beta1   <- time_row[1, "Estimate"]
      se      <- time_row[1, "Std. Error"]
      t_stat  <- time_row[1, "t value"]
      df_sw   <- if ("df" %in% names(time_row)) time_row[1, "df"] else NA_real_
      pval    <- if ("Pr(>|t|)" %in% names(time_row)) time_row[1, "Pr(>|t|)"] else {
        2 * pt(-abs(t_stat), df = df_sw)
      }
      stat_label <- "t_stat"
      stat_val   <- t_stat
    } else {
      # z-test (large-sample)
      cs      <- as.data.frame(coef(summary(fit)))
      beta1   <- cs["visit", "Estimate"]
      se      <- cs["visit", "Std. Error"]
      z_stat  <- beta1 / se
      pval    <- 2 * pnorm(-abs(z_stat))
      df_sw   <- NA_real_
      stat_label <- "z_stat"
      stat_val   <- z_stat
    }

    ci_lower <- beta1 - qnorm(0.975) * se
    ci_upper <- beta1 + qnorm(0.975) * se

    data.frame(
      PC        = pc,
      beta1     = round(beta1,    5),
      se        = round(se,       5),
      statistic = round(stat_val, 4),
      df        = round(df_sw,    2),
      p_value   = signif(pval,    4),
      ci_lower  = round(ci_lower, 5),
      ci_upper  = round(ci_upper, 5),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  class(out) <- c("lomda_wald", "data.frame")
  out
}

# Needed for use_t = FALSE z-test path
pnorm <- stats::pnorm
pt    <- stats::pt

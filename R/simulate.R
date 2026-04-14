#' Simulate longitudinal omics data for lomda
#'
#' Generates a synthetic longitudinal omics dataset suitable for testing and
#' demonstrating the \eqn{lomda} two-stage pipeline. The data structure mirrors a
#' typical metabolomics study: multiple subjects measured at several visits,
#' with a latent time trend embedded in the first few features.
#'
#' @details
#' The simulation scheme:
#' \enumerate{
#'   \item Subject-level random intercepts \eqn{b_i \sim N(0, \sigma_b^2)}
#'     are drawn for each subject.
#'   \item For each subject \eqn{i} and visit \eqn{j}:
#'     \eqn{x_{ij,k} = \mu_k + \lambda_k b_i + \beta_k \cdot t_j + \epsilon_{ij,k}}
#'     where \eqn{\lambda_k} are feature-level loadings (large for the first
#'     \code{n_signal} features) and \eqn{\beta_k} is a time trend (non-zero
#'     for the first \code{n_signal} features).
#'   \item Age is generated as a cross-sectional covariate
#'     \eqn{\text{age}_i \sim N(45, 10^2)}.
#' }
#'
#' @param n_subjects Integer. Number of subjects. Default \code{50}.
#' @param n_visits Integer. Number of visits per subject. Default \code{3}.
#' @param n_features Integer. Total number of omics features. Default \code{50}.
#' @param n_signal Integer. Number of features with a true time trend.
#'   Default \code{5}.
#' @param time_effect Numeric scalar. Slope of the time trend on the signal
#'   features. Default \code{0.5}.
#' @param sigma_b Numeric. Standard deviation of subject random intercepts.
#'   Default \code{1.5}.
#' @param sigma_e Numeric. Residual standard deviation. Default \code{1}.
#' @param seed Integer. Random seed for reproducibility. Default \code{2024}.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{\code{ID}}{Subject identifier (character, e.g. \code{"S001"}).}
#'   \item{\code{visit}}{Visit number (integer 1, 2, ..., \code{n_visits}).}
#'   \item{\code{age}}{Subject age (numeric).}
#'   \item{\code{M1}, \code{M2}, ...}{Omics feature values.}
#' }
#' The data frame has \code{n_subjects * n_visits} rows.
#'
#' @examples
#' dat <- simulate_lomda_data(n_subjects = 30, n_visits = 3,
#'                            n_features = 20, seed = 1)
#' dim(dat)
#' head(dat[, 1:7])
#'
#' @export
simulate_lomda_data <- function(n_subjects  = 50,
                                n_visits    = 3,
                                n_features  = 50,
                                n_signal    = 5,
                                time_effect = 0.5,
                                sigma_b     = 1.5,
                                sigma_e     = 1.0,
                                seed        = 2026) {

  set.seed(seed)
  if (n_signal > n_features)
    stop("n_signal cannot exceed n_features.")

  # Subject-level attributes
  subject_ids <- sprintf("S%03d", seq_len(n_subjects))
  ages        <- rnorm(n_subjects, mean = 45, sd = 10)
  rand_int    <- rnorm(n_subjects, mean = 0,  sd = sigma_b)

  # Feature-level loadings and time effects
  lambda <- c(rep(1.5, n_signal), rep(0.2, n_features - n_signal))
  beta_k <- c(rep(time_effect, n_signal), rep(0, n_features - n_signal))
  mu_k   <- rnorm(n_features, mean = 0, sd = 2)

  rows <- vector("list", n_subjects * n_visits)
  # rows <- vector("character", n_subjects * n_visits)
  idx  <- 1L

  for (i in seq_len(n_subjects)) {
    for (j in seq_len(n_visits)) {
      noise  <- rnorm(n_features, mean = 0, sd = sigma_e)
      omics  <- mu_k + lambda * rand_int[i] + beta_k * j + noise
      # rows[[idx]] <- c(
      #   list(ID = subject_ids[i], time = j, age = ages[i] + 2 * (j - 1)),
      #   setNames(as.list(omics), paste0("M", seq_len(n_features)))
      # )
      rows[[idx]] <- c(
        list(
          ID   = as.character(subject_ids[i]),
          visit = j,
          age  = ages[i] + 2 * (j - 1)
        ),
        setNames(as.list(omics), paste0("M", seq_len(n_features)))
      )
      # rows[idx] <- c(
      #   list(ID = subject_ids[i], time = j, age = ages[i]),
      #   setNames(as.list(omics), paste0("M", seq_len(n_features)))
      # )
      idx <- idx + 1L
    }
  }
  out <- bind_rows(rows)
  out
}

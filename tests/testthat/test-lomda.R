test_that("simulate_lomda_data returns correct structure", {
  dat <- simulate_lomda_data(n_subjects = 10, n_visits = 3,
                             n_features = 15, seed = 1)
  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 10 * 3)
  expect_true("ID" %in% names(dat))
  expect_true("visit" %in% names(dat))
  expect_true("age" %in% names(dat))
  expect_true("M1" %in% names(dat))
  expect_equal(ncol(dat), 3 + 15)  # ID, visit, age + features
  expect_equal(sort(unique(dat$visit)), 1:3)
})

test_that("lomda() fits without error on simulated data", {
  dat <- simulate_lomda_data(n_subjects = 20, n_visits = 3,
                             n_features = 10, seed = 2)
  fit <- expect_no_error(lomda(dat, n_pc = 2, time = "visit"))
  expect_s3_class(fit, "lomda")
  expect_equal(fit$stage2, "lmm")
  expect_equal(fit$time, "visit")
  expect_equal(fit$n_pc, 2)
  expect_named(fit$fit_lmm, c("PC1", "PC2"))
  expect_named(fit$fit_lmm_null, c("PC1", "PC2"))
})

test_that("lomda() supports adjustment covariates for visit LMM", {
  dat <- simulate_lomda_data(n_subjects = 20, n_visits = 3,
                             n_features = 10, seed = 12)
  dat$batch <- rep(c("A", "B"), length.out = nrow(dat))
  fit <- expect_no_error(lomda(dat, n_pc = 2, time = "visit",
                               adjust = "batch"))
  expect_equal(fit$adjust, "batch")
  expect_true("batch" %in% names(fit$scores))
  expect_false("batch" %in% fit$omics_cols)
})

test_that("lomda() supports user-supplied metadata column names", {
  dat <- simulate_lomda_data(n_subjects = 20, n_visits = 3,
                             n_features = 10, seed = 17)
  names(dat)[1:3] <- c("subject_id", "wave", "age_years")

  fit_visit <- expect_no_error(
    lomda(dat, ID = "subject_id", visit = "wave", age = "age_years",
          n_pc = 2, time = "visit")
  )
  expect_equal(fit_visit$id_col, "subject_id")
  expect_equal(fit_visit$time_col, "wave")
  expect_true(all(c("subject_id", "wave", "age_years") %in%
                    names(fit_visit$scores)))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(plot_scores(fit_visit), "ggplot")
    expect_s3_class(plot_trajectory(fit_visit), "ggplot")
  }

  fit_age <- expect_no_error(
    lomda(dat, ID = subject_id, visit = wave, age = age_years,
          n_pc = 2, time = "age", method_fda = "spline",
          grid_length = 25)
  )
  expect_equal(fit_age$id_col, "subject_id")
  expect_equal(fit_age$age_col, "age_years")
  pred <- predict(fit_age, ages = c(35, 45, 55))
  expect_true(all(c("age_years", "PC1", "PC2") %in% names(pred)))
})

test_that("lomda() routes age analysis to FDA", {
  dat <- simulate_lomda_data(n_subjects = 20, n_visits = 3,
                             n_features = 10, seed = 13)
  fit <- expect_no_error(lomda(dat, n_pc = 2, time = "age",
                               method_fda = "spline", grid_length = 25))
  expect_s3_class(fit, "lomda_fda")
  expect_s3_class(fit, "lomda")
  expect_equal(fit$stage2, "fda")
  expect_equal(fit$time, "age")
  expect_named(fit$fda_fits, c("PC1", "PC2"))
  pred <- predict(fit, ages = c(35, 45, 55))
  expect_equal(nrow(pred), 3)
  expect_s3_class(plot(fit, type = "variance"), "ggplot")
})

test_that("lomda() rejects FDA method for visit LMM", {
  dat <- simulate_lomda_data(n_subjects = 20, n_visits = 3,
                             n_features = 10, seed = 14)
  expect_error(lomda(dat, n_pc = 2, time = "visit", method_fda = "face"),
               "method_fda is only used")
})

test_that("lomda() rejects unknown FDA method", {
  dat <- simulate_lomda_data(n_subjects = 20, n_visits = 3,
                             n_features = 10, seed = 15)
  expect_error(lomda(dat, n_pc = 2, time = "age", method_fda = "bad"),
               "should be one of")
})

test_that("lomda_pca returns correct components", {
  dat <- simulate_lomda_data(n_subjects = 15, n_features = 12, seed = 3)
  res <- lomda_pca(dat, n_pc = 4)
  expect_true(all(c("pca", "scores", "var_explained", "loadings") %in%
                    names(res)))
  expect_equal(ncol(res$loadings), 4)
  expect_equal(length(res$var_explained), 4)
  expect_true(all(res$var_explained >= 0))
  expect_true(sum(res$var_explained) <= 1 + 1e-9)
})

test_that("lomda_lmm returns fits and null_fits", {
  dat <- simulate_lomda_data(n_subjects = 15, n_features = 10, seed = 4)
  pca_res <- lomda_pca(dat, n_pc = 2)
  lmm_res <- lomda_lmm(pca_res$scores, n_pc = 2)
  expect_named(lmm_res, c("fits", "null_fits"))
  expect_named(lmm_res$fits,      c("PC1", "PC2"))
  expect_named(lmm_res$null_fits, c("PC1", "PC2"))
})

test_that("lomda_lrt returns a data frame with expected columns", {
  dat <- simulate_lomda_data(n_subjects = 20, n_features = 10, seed = 5)
  fit <- lomda(dat, n_pc = 2, REML = FALSE)
  lrt <- lomda_lrt(fit)
  expect_s3_class(lrt, "data.frame")
  expect_equal(nrow(lrt), 2)
  expect_true(all(c("PC", "LRT_stat", "p_value") %in% names(lrt)))
  expect_true(all(lrt$LRT_stat >= 0))
  expect_true(all(lrt$p_value >= 0) && all(lrt$p_value <= 1))
})

test_that("lomda_wald returns a data frame with expected columns", {
  dat <- simulate_lomda_data(n_subjects = 20, n_features = 10, seed = 6)
  fit <- lomda(dat, n_pc = 2)
  wald <- lomda_wald(fit)
  expect_s3_class(wald, "data.frame")
  expect_equal(nrow(wald), 2)
  expect_true(all(c("PC", "beta1", "se", "statistic", "p_value",
                    "ci_lower", "ci_upper") %in% names(wald)))
  # ci_lower < beta1 < ci_upper
  expect_true(all(wald$ci_lower < wald$beta1))
  expect_true(all(wald$beta1 < wald$ci_upper))
})

test_that("lomda_wald z-test option works", {
  dat <- simulate_lomda_data(n_subjects = 20, n_features = 10, seed = 7)
  fit <- lomda(dat, n_pc = 2)
  wald_z <- lomda_wald(fit, use_t = FALSE)
  expect_true(all(is.na(wald_z$df)))
})

test_that("plot functions return ggplot objects", {
  skip_if_not_installed("ggplot2")
  dat <- simulate_lomda_data(n_subjects = 20, n_features = 10, seed = 8)
  fit <- lomda(dat, n_pc = 3)
  expect_s3_class(plot_scores(fit),             "ggplot")
  expect_s3_class(plot_loadings(fit),           "ggplot")
  expect_s3_class(plot_loadings(fit, pc_y=NULL),"ggplot")
  expect_s3_class(plot_variance_explained(fit), "ggplot")
  expect_s3_class(plot_trajectory(fit),         "ggplot")
})

test_that("lomda_fda estimates trajectories and predicts on age", {
  dat <- simulate_lomda_data(n_subjects = 20, n_features = 10, seed = 10)
  fda <- lomda_fda(dat, n_pc = 2, method = "spline", grid_length = 25)
  expect_s3_class(fda, "lomda_fda")
  expect_named(fda$fda_fits, c("PC1", "PC2"))
  pred <- predict(fda, ages = c(35, 45, 55))
  expect_s3_class(pred, "data.frame")
  expect_equal(nrow(pred), 3)
  expect_true(all(c("age", "PC1", "PC2") %in% names(pred)))
  imp <- lomda_fda_important(fda, pc = 1, n_top = 5)
  expect_equal(nrow(imp), 5)
  expect_true(all(c("feature", "importance", "rank") %in% names(imp)))
})

test_that("FDA plot functions return ggplot objects", {
  skip_if_not_installed("ggplot2")
  dat <- simulate_lomda_data(n_subjects = 20, n_features = 10, seed = 11)
  fda <- lomda_fda(dat, n_pc = 2, method = "spline")
  expect_s3_class(plot_fda_trajectory(fda, pcs = 1:2, n_subjects = 5), "ggplot")
  expect_s3_class(plot_fda_importance(fda), "ggplot")
})

test_that("importance ranking works for LMM and FDA fits", {
  dat <- simulate_lomda_data(n_subjects = 20, n_features = 10, seed = 16)
  fit_visit <- lomda(dat, n_pc = 2, time = "visit")
  imp_visit <- lomda_important(fit_visit, n_top = 5)
  expect_equal(nrow(imp_visit), 5)
  expect_true(all(c("feature", "importance", "signed_importance", "rank") %in%
                    names(imp_visit)))

  fit_age <- lomda(dat, n_pc = 2, time = "age", method_fda = "spline")
  imp_age <- lomda_important(fit_age, n_top = 5)
  expect_equal(nrow(imp_age), 5)
  expect_true(all(c("feature", "importance", "signed_importance", "rank") %in%
                    names(imp_age)))
})

test_that("lomda errors on bad input", {
  dat <- simulate_lomda_data(n_subjects = 10, n_features = 5, seed = 9)
  # Missing ID column
  expect_error(lomda(dat[, -1], n_pc = 2), "ID")
  # n_pc too large
  expect_error(lomda(dat, n_pc = 100), "n_pc")
  # Bad covariate
  expect_error(lomda(dat, adjust = "nonexistent"), "nonexistent")
})

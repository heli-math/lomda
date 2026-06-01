#' @importFrom utils tail
NULL

#' Plot method for lomda objects
#'
#' Convenience wrapper that dispatches to specialized plot functions. The
#' \code{type} argument selects which plot to produce.
#'
#' @param x A \code{lomda} object.
#' @param type Character vector. One or more of \code{"scores"} (PC1 vs PC2 scatter),
#'   \code{"loadings"} (feature loadings bar chart), \code{"variance"}
#'   (scree / variance-explained plot), or \code{"trajectory"} (mean PC
#'   score over visit). Use \code{"all"} to show all plots. Default
#'   \code{"all"}.
#' @param pc_x Integer. PC for the x-axis in score/loading plots. Default 1.
#' @param pc_y Integer. PC for the y-axis in score/loading plots. Default 2.
#' @param color_by Character. Column name in the score data frame to use for
#'   colouring points (e.g., \code{"visit"} or \code{"age"}).
#'   Default \code{"visit"}.
#' @param n_top Integer. Number of top-loading features to label in the
#'   loadings plot. Default \code{10}.
#' @param pause Logical. If \code{TRUE}, wait for Return between plots when
#'   multiple plots are requested.
#' @param ... Additional arguments passed to underlying plot functions.
#'
#' @return A \code{ggplot} object or list of \code{ggplot} objects (invisibly).
#'
#' @examples
#' dat <- simulate_lomda_data(seed = 42)
#' fit <- lomda(dat, n_pc = 3)
#' plot(fit, type = "scores")
#' plot(fit, type = "variance")
#' plot(fit, type = "loadings")
#' plot(fit, type = "trajectory")
#'
#' @export
plot.lomda <- function(x, type = "all",
                       pc_x = 1L, pc_y = 2L,
                       color_by = "visit",
                       n_top = 10L,
                       pause = interactive(), ...) {
  valid_types <- c("scores", "loadings", "variance", "trajectory")
  if (identical(type, "all")) {
    type <- valid_types
  }
  type <- match.arg(type, valid_types, several.ok = TRUE)

  plots <- lapply(type, function(one_type) {
    p <- switch(one_type,
      scores     = plot_scores(x,     pc_x = pc_x, pc_y = pc_y,
                               color_by = color_by, ...),
      loadings   = plot_loadings(x,   pc_x = pc_x, pc_y = pc_y,
                                 n_top = n_top, ...),
      variance   = plot_variance_explained(x, ...),
      trajectory = plot_trajectory(x, ...)
    )
    print(p)
    if (pause && one_type != tail(type, 1)) {
      readline("Press Return to show the next plot...")
    }
    p
  })

  invisible(if (length(plots) == 1) plots[[1]] else plots)
}


#' PC Score Plot (PC1 vs PC2 or any two PCs)
#'
#' Produces a scatter plot of PC scores, coloured by a metadata variable.
#' Each point is one observation (subject x visit).
#'
#' @param x A \code{lomda} object, or a scores data frame (output of
#'   \code{\link{lomda_pca}}).
#' @param pc_x Integer. Index of the PC for the x-axis. Default \code{1}.
#' @param pc_y Integer. Index of the PC for the y-axis. Default \code{2}.
#' @param color_by Character. Column to colour by. Default \code{"visit"}.
#' @param label_ids Logical. Whether to add subject ID labels with
#'   \pkg{ggrepel}. Default \code{FALSE}.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' dat <- simulate_lomda_data(seed = 1)
#' fit <- lomda(dat)
#' plot_scores(fit)
#' plot_scores(fit, pc_x = 1, pc_y = 3, color_by = "age")
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw scale_color_brewer
#'   scale_color_manual theme element_text
#' @importFrom ggrepel geom_text_repel
#' @export
plot_scores <- function(x,
                        pc_x      = 1L,
                        pc_y      = 2L,
                        color_by  = "visit",
                        label_ids = FALSE) {

  scores_df <- if (inherits(x, "lomda")) x$scores else x
  var_exp   <- if (inherits(x, "lomda")) x$var_explained else NULL

  xvar <- paste0("PC", pc_x)
  yvar <- paste0("PC", pc_y)

  if (!xvar %in% names(scores_df))
    stop(xvar, " not found in score data frame.")
  if (!yvar %in% names(scores_df))
    stop(yvar, " not found in score data frame.")

  if (!color_by %in% names(scores_df))
    stop(color_by, " not found in score data frame.")

  scores_df[[color_by]] <- if (is.numeric(scores_df[[color_by]]) &&
                                length(unique(scores_df[[color_by]])) > 6)
    scores_df[[color_by]]
  else
    factor(scores_df[[color_by]])

  x_lab <- if (!is.null(var_exp) && xvar %in% names(var_exp))
    sprintf("%s (%.1f%%)", xvar, var_exp[xvar] * 100) else xvar
  y_lab <- if (!is.null(var_exp) && yvar %in% names(var_exp))
    sprintf("%s (%.1f%%)", yvar, var_exp[yvar] * 100) else yvar

  p <- ggplot2::ggplot(scores_df,
         ggplot2::aes(x = .data[[xvar]], y = .data[[yvar]],
                      color = .data[[color_by]])) +
    ggplot2::geom_point(size = 2.5, alpha = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "grey60", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        color = "grey60", linewidth = 0.4) +
    ggplot2::labs(x = x_lab, y = y_lab, color = color_by,
                  title = paste("PC Score Plot:", xvar, "vs", yvar)) +
    .lomda_theme() +
    ggplot2::theme(legend.position = "right")

  if (is.factor(scores_df[[color_by]])) {
    p <- p + ggplot2::scale_color_brewer(palette = "Set2")
  }

  if (label_ids && "ID" %in% names(scores_df)) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = .data[["ID"]]),
      size = 2.5, max.overlaps = 20
    )
  }
  p
}


#' Feature Loadings Plot
#'
#' Plots a biplot-style bar chart of feature loadings for one or two PCs,
#' highlighting the top-loading features.
#'
#' @param x A \code{lomda} object.
#' @param pc_x Integer. PC index for x-axis loadings. Default \code{1}.
#' @param pc_y Integer. PC index for y-axis loadings (biplot). If \code{NULL},
#'   a simple bar chart of \code{pc_x} loadings is produced.
#' @param n_top Integer. Number of features to label / highlight. Default 10.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' dat <- simulate_lomda_data(seed = 5)
#' fit <- lomda(dat, n_pc = 3)
#' plot_loadings(fit)
#' plot_loadings(fit, pc_x = 1, pc_y = NULL)
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_point geom_text labs
#'   theme_bw scale_fill_gradient2 coord_fixed
#' @importFrom ggrepel geom_text_repel
#' @export
plot_loadings <- function(x, pc_x = 1L, pc_y = 2L, n_top = 10L) {
  if (!inherits(x, "lomda"))
    stop("x must be a 'lomda' object.")

  rot <- as.data.frame(x$fit_pca$rotation)
  xvar <- paste0("PC", pc_x)
  yvar <- if (!is.null(pc_y)) paste0("PC", pc_y) else NULL

  if (!xvar %in% names(rot))
    stop(xvar, " not in rotation matrix.")

  if (!is.null(yvar) && !yvar %in% names(rot))
    stop(yvar, " not in rotation matrix.")

  rot$feature <- rownames(rot)
  rot$loading_abs <- abs(rot[[xvar]])

  if (is.null(yvar)) {
    # Bar chart for one PC
    rot <- rot[order(-rot$loading_abs), ]
    rot$feature <- factor(rot$feature, levels = rot$feature)
    p <- ggplot2::ggplot(rot,
           ggplot2::aes(x = .data$feature, y = .data[[xvar]],
                        fill = .data[[xvar]])) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white",
                                    high = "#d6604d", midpoint = 0) +
      ggplot2::labs(x = "Feature", y = paste(xvar, "loading"),
                    title = paste("Feature Loadings:", xvar),
                    fill = "Loading") +
      .lomda_theme(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                    hjust = 1, size = 7))
  } else {
    # Biplot scatter of loadings
    mag <- sqrt(rot[[xvar]]^2 + rot[[yvar]]^2)
    top_feats <- head(rot$feature[order(-mag)], n_top)
    rot$top   <- rot$feature %in% top_feats

    p <- ggplot2::ggplot(rot,
           ggplot2::aes(x = .data[[xvar]], y = .data[[yvar]])) +
      ggplot2::geom_point(ggplot2::aes(color = .data$top), size = 1.8,
                          alpha = 0.7) +
      ggrepel::geom_text_repel(
        data = rot[rot$top, ],
        ggplot2::aes(label = .data$feature),
        size = 3, max.overlaps = 30
      ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "grey60") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                          color = "grey60") +
      ggplot2::scale_color_manual(values = c("FALSE" = "grey70",
                                             "TRUE"  = "#d6604d"),
                                  guide  = "none") +
      ggplot2::labs(x = paste(xvar, "loading"),
                    y = paste(yvar, "loading"),
                    title = "Feature Loadings Biplot") +
      .lomda_theme()
  }
  p
}


#' Trajectory Plot: Mean PC Score Over Time
#'
#' Plots the mean PC score at each visit (+- standard error), with individual
#' subject trajectories overlaid as thin lines. Useful for visualising the
#' time trend estimated by the Stage 2 LMM.
#'
#' @param x A \code{lomda} object.
#' @param pcs Integer vector. Which PCs to plot. Default \code{1:3}.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' dat <- simulate_lomda_data(seed = 10)
#' fit <- lomda(dat, n_pc = 3)
#' plot_trajectory(fit)
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_bw facet_wrap
#' @importFrom dplyr group_by summarise mutate
#' @importFrom tidyr pivot_longer
#' @export
plot_trajectory <- function(x, pcs = 1:3) {
  if (!inherits(x, "lomda"))
    stop("x must be a 'lomda' object.")

  pc_names <- paste0("PC", pcs[pcs <= x$n_pc])
  if (length(pc_names) == 0)
    stop("None of the requested PCs exist in the model.")

  scores_df <- x$scores

  long <- tidyr::pivot_longer(
    scores_df[, c("ID", "visit", pc_names)],
    cols      = dplyr::all_of(pc_names),
    names_to  = "PC",
    values_to = "score"
  )
  long$visit <- as.integer(long$visit)

  # Mean +- SE per time point per PC
  summ <- dplyr::group_by(long, .data$PC, .data$visit)
  summ <- dplyr::summarise(summ,
    mean_score = mean(.data$score),
    se_score   = sd(.data$score) / sqrt(dplyr::n()),
    .groups    = "drop"
  )
  visit_breaks <- sort(unique(long$visit))

  ggplot2::ggplot() +
    # Individual trajectories
    ggplot2::geom_line(
      data = long,
      ggplot2::aes(x = .data$visit, y = .data$score,
                   group = .data$ID),
      color = "grey75", linewidth = 0.35, alpha = 0.6
    ) +
    # Mean trajectory
    ggplot2::geom_line(
      data = summ,
      ggplot2::aes(x = .data$visit, y = .data$mean_score),
      color = "#2166ac", linewidth = 1.2
    ) +
    ggplot2::geom_point(
      data = summ,
      ggplot2::aes(x = .data$visit, y = .data$mean_score),
      color = "#2166ac", size = 3
    ) +
    ggplot2::geom_line(
      data = summ,
      ggplot2::aes(x = .data$visit,
                   y = .data$mean_score + 1.96 * .data$se_score),
      linetype = "dashed", color = "#2166ac", linewidth = 0.7
    ) +
    ggplot2::geom_line(
      data = summ,
      ggplot2::aes(x = .data$visit,
                   y = .data$mean_score - 1.96 * .data$se_score),
      linetype = "dashed", color = "#2166ac", linewidth = 0.7
    ) +
    ggplot2::facet_wrap(~ .data$PC, scales = "free_y") +
    ggplot2::scale_x_continuous(breaks = visit_breaks) +
    ggplot2::labs(
      x = "Visit", y = "PC Score",
      title = "PC Score Trajectory Over Time",
      subtitle = "Blue = group mean +- 95% CI; grey = individual trajectories"
    ) +
    .lomda_theme()
}


#' Variance Explained Plot (Scree Plot)
#'
#' Displays the proportion (and cumulative proportion) of variance explained
#' by each PC extracted in Stage 1.
#'
#' @param x A \code{lomda} object.
#' @param max_pc Integer. Maximum number of PCs to show. Default: all PCs
#'   in the \code{prcomp} object (up to 20).
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' dat <- simulate_lomda_data(seed = 3)
#' fit <- lomda(dat, n_pc = 5)
#' plot_variance_explained(fit)
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_line geom_point labs theme_bw
#' @export
plot_variance_explained <- function(x, max_pc = NULL) {
  if (!inherits(x, "lomda"))
    stop("x must be a 'lomda' object.")

  sdev <- x$fit_pca$sdev
  max_pc <- min(max_pc %||% length(sdev), 20L, length(sdev))
  ve   <- (sdev^2 / sum(sdev^2))[seq_len(max_pc)]
  cum_ve <- cumsum(ve)

  df <- data.frame(
    PC       = seq_len(max_pc),
    var_exp  = ve,
    cum_var  = cum_ve
  )
  df$is_modelled <- df$PC <= x$n_pc

  ggplot2::ggplot(df, ggplot2::aes(x = .data$PC)) +
    ggplot2::geom_col(
      ggplot2::aes(y = .data$var_exp,
                   fill = .data$is_modelled),
      width = 0.7
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$cum_var), color = "#d6604d",
      linewidth = 1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$cum_var), color = "#d6604d", size = 2
    ) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#2166ac", "FALSE" = "grey70"),
      labels = c("TRUE" = "Modelled", "FALSE" = "Not modelled"),
      name   = NULL
    ) +
    ggplot2::labs(
      x     = "Principal Component",
      y     = "Proportion of Variance",
      title = "Variance Explained by PCA",
      subtitle = "Red line = cumulative; blue bars = PCs used in Stage 2 LMM"
    ) +
    .lomda_theme()
}

# Utility: null-coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a

# Shared visual style for package plots.
.lomda_theme <- function(base_size = 13) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", color = "#243447"),
      plot.subtitle = ggplot2::element_text(color = "#52616B"),
      panel.grid.major = ggplot2::element_line(color = "#E7ECEF", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#F4F7F8", color = "#D8E0E3"),
      strip.text = ggplot2::element_text(face = "bold", color = "#243447")
    )
}

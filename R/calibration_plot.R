#' Calibration plots for cbenefit objects
#'
#' @param x A cbenefit object
#' @param limits Optional: a vector of length two with y axis limits
#'
#' @return A list with two ggplot2 objects
#' @export
#'
#' @examples
#' m1 <- model_fit(Surv(time, status) ~  treat + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
#' dat_toy, model = "coxph")
#
#' bf <- benefit_cv(m1, "treat", folds = 5, times = 3)
#' calibration_plot(bf)
calibration_plot <- function(x, limits) {

  #browser()
  if(!inherits(x, "cbenefit")) warning("x must be a cbenefit object, see ?benefit_cv")

  if(missing(limits))
    limits <- c(min(x$obs_pred_nocv$mean_predicted, x$obs_pred_nocv$mean_y, -0.2),
                     max(x$obs_pred_nocv$mean_predicted, x$obs_pred_nocv$mean_y, 0.3))

  # if(missing(data)) {
  #   stop("must have data")
  # }
  #browser()

  calplot_nocv <-
    x$obs_pred_nocv %>%
    ggplot(aes(x = .data$mean_predicted, y = .data$mean_y)) +
    geom_point() +
    geom_linerange(aes(ymin = .data$mean_y - 1.96 * .data$sd_y, ymax = .data$mean_y + 1.96 * .data$sd_y)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(xlim = limits, ylim = limits) +
    labs(x = "Predicted benefit", y = "Observed benefit", title = "No CV") +
    theme_classic()

  calplot_cv <- NA
  if(inherits(x$obs_pred_cv, "data.frame"))
  calplot_cv <-
    x$obs_pred_cv %>%
    ggplot(aes(x = .data$mean_predicted, y = .data$mean_y)) +
    geom_point() +
    geom_linerange(aes(ymin = .data$mean_y - 1.96 * .data$sd_y, ymax = .data$mean_y + 1.96 * .data$sd_y)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_abline(intercept = 0, slope = 1) +
    coord_cartesian(xlim = limits, ylim = limits) +
    labs(x = "Predicted benefit", y = "Observed benefit", title = "CV") +
    theme_classic()

  list(caplot_nocv = calplot_nocv, caplot_cv = calplot_cv)
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

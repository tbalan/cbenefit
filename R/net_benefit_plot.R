#' Net benefit plots for cbenefit objects
#'
#' @param x A cbenefit object
#' @param limits Optional: a vector of length two with y axis limits
#'
#' @return A list with two plots: one before cross validation, one after cross validation
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_abline labs coord_cartesian geom_hline geom_linerange geom_vline theme_classic
#'
#' @examples
#' m2 <- model_fit(status ~  treat + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
#' dat_toy, model = "glm")
#
#' bf <- benefit_cv(m2, "treat", folds = 5, times = 3)
#' net_benefit_plot(bf)
net_benefit_plot <- function(x, limits) {

  # browser()

  if(missing(limits))
    ylim <-    c((x$event_rate_treated - x$event_rate_controls)[1],
                   (x$event_rate_controls - x$event_rate_treated)[1] - min(x$threshold))


  if(missing(x)) {
    stop("must have x")
  }

  net_benefit_plot_nocv <-
    data.frame(threshold = x$threshold, net_benefit = x$nb_nocv) %>%
    ggplot(aes(x = .data$threshold, y = .data$net_benefit)) +
    geom_point() +
    geom_line() +
    geom_line(aes(y = x$event_rate_controls - x$event_rate_treated - x$threshold)) +
    geom_abline(slope = 0) +
    labs(x = "Benefit threshold",
         y = "Net benefit",
         title = "No CV") +
    coord_cartesian(ylim = ylim) +
    theme_classic()

  net_benefit_plot_cv <- NA
  if(inherits(x$nb_cv, "numeric"))
  net_benefit_plot_cv <-
    data.frame(threshold = x$threshold, net_benefit = x$nb_cv) %>%
    ggplot(aes(x = .data$threshold, y = .data$net_benefit)) +
    geom_point() +
    geom_line() +
    geom_line(aes(y = x$event_rate_controls - x$event_rate_treated - x$threshold)) +
    geom_abline(slope = 0) +
    labs(x = "Benefit threshold",
         y = "Net benefit",
         title = "CV") +
    coord_cartesian(ylim = ylim) +
    theme_classic()

  list(nbplot_nocv = net_benefit_plot_nocv, nbplot_cv = net_benefit_plot_cv)

}

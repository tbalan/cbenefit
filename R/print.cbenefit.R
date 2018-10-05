#' Print cbenefit objects
#'
#' @param x An object of class cbenefit
#' @param ... Does nothing at this point
#'
#' @return Nothing
#' @export
#'
#' @examples
#' m1 <- model_fit(Surv(time, status) ~  treat + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
#' dat_toy, model = "coxph")
#
#' cb <- benefit_cv(m1, "treat", folds = 5, times = 3)
#' print(cb)
print.cbenefit <- function(x, ...) {
  # browser()

  cat("Object of type cbenefit")
  cat("\nEvent rate control group:", round(x$event_rate_controls, digits = 3))
  cat("\nEvent rate treatment group:", round(x$event_rate_treated, digits = 3))

  cat("\nc-statistics:\n")
  print(x$cindex_nocv)

  cat("\naverage c-statistics (from", nrow(x$c_index_cv),"CV):\n")
  print(apply(x$c_index_cv, 2, mean))

    cat("\nc-statistics on averaged probabilities (from",nrow(x$c_index_cv),"CV):\n")
  print(x$cindex_cv_ave)


  invisible(x)
}

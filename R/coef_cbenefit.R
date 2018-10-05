#' Printing coefficients for c_glmnet fits
#'
#' @param mod A model as fitted by model_fit
#' @param ... Further argumnets passed to tidy()
#'
#' @return A data frame with the estimates and confidence intervals if available
#' @export
#'
#' @importFrom tibble rownames_to_column
#' @importFrom broom tidy
#' @importFrom dplyr rename
#'
#' @examples
#' m1 <- model_fit(Surv(time, status) ~  treat + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
#' dat_toy, model = "coxph")
#' coef_cbenefit(m1)

coef_cbenefit <- function(mod, ...) {
  if(inherits(mod, "coxph") | inherits(mod, "glm"))
    return(broom::tidy(mod, ...))

  # browser()
  if(inherits(mod, "c_glmnet")) {
    est <- coef(mod$glmnetfit, s = "lambda.min")
    return(as.data.frame(as.matrix(est)) %>% rownames_to_column(var = "term") %>% rename(estimate = .data$`1`))

  }
}

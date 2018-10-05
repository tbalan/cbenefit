#' Model fit for cbenefit
#'
#' @param formula A formula to be evaluated for the model
#' @param data A data frame or the name of the data frame where the formula may be evaluated
#' @param model One of "coxph" for a Cox model or "glm" for logistic regression
#' @param ... Further argumnets that are passed on to `glm` or `coxph`
#'
#' @return A model fit of type coxph or glm
#' @export
#'
#' @seealso model_fit_glmnet
#'
#' @importFrom survival coxph Surv
#'
#' @details The function does not spend time checking its input. In particular, the function has only been tested with complete case data frames, and where the treatment is coded as 0 (control) and 1 (treated).
#' @examples
#' model_fit(Surv(time, status) ~  treat + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
#' dat_toy, model = "coxph")
#'
#' model_fit(status ~ treat +  x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
#'           dat_toy, model = "glm")
model_fit <- function(formula, data, model, ...) {

  # idea: you can pass a data.frame or the name of a data.frame that lies in the parent environemnt
  if(!inherits(data, "data.frame"))
    data = get(data, parent.frame())

  # browser()
  if(model == "coxph") {

    fit <- survival::coxph(formula, data = data, model = TRUE, ...)

  }

  if(model == "glm") {
    fit <- glm(formula = formula, family = binomial(), data = data, ...)
  }

  fit
}


#' Toy data with treatment column and outcome
#'
"dat_toy"


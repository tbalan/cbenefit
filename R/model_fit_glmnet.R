#' Fit glmnet models for the cbenefit package
#'
#' @param formula A formula to be evaluated for the model
#' @param data A data frame or the name of the data frame where the formula may be evaluated
#' @param model One of "coxph" for a Cox model or "glm" for logistic regression
#' @param pen One of "ridge" or "lasso"
#' @param treat_name Name of the treatment column as string
#'
#' @return An object of type c_glmnet
#' @export
#'
#' @seealso model_fit
#'
#' @details A model of class c_glmnet is a named list where: x is the model matrix, x_other is the model matrix with the treatment allocation reversed,
#' y is the outcome, alpha is the glmnet penalization parameter (0 for ridge, 1 for lasso) and glmnetfit is the result of glmnet::cv.glmnet.
#' This augmented object is needed for further calculations in the cbenefit package.
#'
#' The function does not spend time checking its input. In particular, the function has only been tested with complete case data frames, and where the treatment is coded as 0 (control) and 1 (treated).
#' @importFrom dplyr mutate
#' @importFrom glmnet cv.glmnet
#' @importFrom rlang sym
#' @examples
#' model_fit_glmnet(Surv(time, status) ~ treat + x1 + x2 + x3,
#' dat_toy, model = "coxph", pen = "ridge", treat_name = "treat")
#'
model_fit_glmnet <- function(formula, data, model, pen = NULL, treat_name) {

  if(!inherits(data, "data.frame"))
    data = get(data, parent.frame())

  x <- model.matrix(formula, data)[ , -1,drop=FALSE]
  y <- model.frame(update(formula, "~1"), data)[ , 1, drop =TRUE]

  treat <- rlang::sym(treat_name)
  x_other <- model.matrix(formula, mutate(data, treat = 1 - treat))[,-1,drop=FALSE]

  # no penalization for the treatment main effect
  penalty <- rep(1, ncol(x))
  penalty[which(colnames(x) == treat_name)] <- 0

  if(pen == "ridge") alpha <- 0
  if(pen == "lasso") alpha <- 1

  if(model == "glm") {

    fit <- list(x = x,
                x_other = x_other,
                y = y, alpha = alpha,
                glmnetfit = cv.glmnet(x, y, family = "binomial", alpha = alpha, penalty.factor = penalty))

  }

  if(model == "coxph") {

    fit <- list(x = x, x_other = x_other, y = y, alpha = alpha, glmnetfit = cv.glmnet(x, y, family = "cox", alpha = alpha, penalty.factor = penalty))

  }

  attr(fit, "class") <- "c_glmnet"
  fit

}


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_fit_glmnet.R
\name{model_fit_glmnet}
\alias{model_fit_glmnet}
\title{Fit glmnet models for the cbenefit package}
\usage{
model_fit_glmnet(formula, data, model, pen = NULL, treat_name)
}
\arguments{
\item{formula}{A formula to be evaluated for the model}

\item{data}{A data frame or the name of the data frame where the formula may be evaluated}

\item{model}{One of "coxph" for a Cox model or "glm" for logistic regression}

\item{pen}{One of "ridge" or "lasso"}

\item{treat_name}{Name of the treatment column as string}
}
\value{
An object of type c_glmnet
}
\description{
Fit glmnet models for the cbenefit package
}
\details{
A model of class c_glmnet is a named list where: x is the model matrix, x_other is the model matrix with the treatment allocation reversed,
y is the outcome, alpha is the glmnet penalization parameter (0 for ridge, 1 for lasso) and glmnetfit is the result of glmnet::cv.glmnet.
This augmented object is needed for further calculations in the cbenefit package.

The function does not spend time checking its input. In particular, the function has only been tested with complete case data frames, and where the treatment is coded as 0 (control) and 1 (treated).
}
\examples{
model_fit_glmnet(Surv(time, status) ~ treat + x1 + x2 + x3,
dat_toy, model = "coxph", pen = "ridge", treat_name = "treat")

}
\seealso{
model_fit
}

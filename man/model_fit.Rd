% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_fit.R
\name{model_fit}
\alias{model_fit}
\title{Model fit for cbenefit}
\usage{
model_fit(formula, data, model, ...)
}
\arguments{
\item{formula}{A formula to be evaluated for the model}

\item{data}{A data frame or the name of the data frame where the formula may be evaluated}

\item{model}{One of "coxph" for a Cox model or "glm" for logistic regression}

\item{...}{Further argumnets that are passed on to `glm` or `coxph`}
}
\value{
A model fit of type coxph or glm
}
\description{
Model fit for cbenefit
}
\details{
The function does not spend time checking its input. In particular, the function has only been tested with complete case data frames, and where the treatment is coded as 0 (control) and 1 (treated).
}
\examples{
model_fit(Surv(time, status) ~  treat + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
dat_toy, model = "coxph")

model_fit(status ~ treat +  x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
          dat_toy, model = "glm")
}
\seealso{
model_fit_glmnet
}

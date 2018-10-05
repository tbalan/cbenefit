## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ">"
)

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github("tbalan/cbenefit")

## ------------------------------------------------------------------------
library(cbenefit)

## ------------------------------------------------------------------------
data("dat_toy")
head(dat_toy)

## ------------------------------------------------------------------------
# Cox model, wrapper around coxph()
model_fit(Surv(time, status) ~  treat + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
dat_toy, model = "coxph")

# Logistic regression, wrapper around glm()
model_fit(status ~ treat +  x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
          dat_toy, model = "glm")

## ---- results="hide"-----------------------------------------------------
# survival, ridge
model_fit_glmnet(Surv(time, status) ~  treat + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
dat_toy, model = "coxph", pen = "ridge", treat_name = "treat")

# survival, lasso
model_fit_glmnet(Surv(time, status) ~  treat + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
dat_toy, model = "coxph", pen = "lasso", treat_name = "treat")

# logistic, ridge
model_fit_glmnet(status ~  treat + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
dat_toy, model = "glm", pen = "ridge", treat_name = "treat")

# logistic, lasso
model_fit_glmnet(status ~  treat + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
dat_toy, model = "glm", pen = "lasso", treat_name = "treat")

## ------------------------------------------------------------------------
m1 <- model_fit(Surv(time, status) ~  treat + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
dat_toy, model = "coxph")

m2 <- model_fit_glmnet(status ~  treat + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
dat_toy, model = "glm", pen = "lasso", treat_name = "treat")

coef_cbenefit(m1)
coef_cbenefit(m2)

## ------------------------------------------------------------------------
mod <- model_fit(Surv(time, status) ~  treat + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
dat_toy, model = "coxph")

## ------------------------------------------------------------------------
bb <- benefit_cv(mod, "treat", folds = 10, times = 5)

bb

## ------------------------------------------------------------------------
bb <- benefit_cv(mod, "treat", folds = 10, times = 1)

bb

## ------------------------------------------------------------------------
bb <- benefit_cv(mod, "treat", folds = 10, times = 1, reproducible = TRUE)

bb

## ------------------------------------------------------------------------
mcv <- benefit_cv(mod, "treat", folds = 10, times = 1, match_by = "covariates", data = dat_toy)
mcv

## ------------------------------------------------------------------------
matchb <- benefit_cv(mod, "treat", folds = 10, times = 1, matchit_args = list(distance = "mahalanobis"))

matchc <- benefit_cv(mod, "treat", folds = 10, times = 1, matchit_args = list(distance = "cloglog"))

## ------------------------------------------------------------------------
cp <- calibration_plot(bb)
nb <- net_benefit_plot(bb)

library(ggplot2)
cp$caplot_nocv 
cp$caplot_cv

nb$nbplot_nocv 
nb$nbplot_cv

## ------------------------------------------------------------------------
library(egg)
ggarrange(plots = cp, nrow = 1)
ggarrange(plots = nb, nrow = 1)

## ------------------------------------------------------------------------
model <- "coxph"
# at the formula, just specify it without interactions, just the variables you want in
formula <- Surv(time, status) ~ treat + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
dataset <- "dat_toy"
treat_name <- "treat"
model = "coxph"

## ------------------------------------------------------------------------
library(rms) # need it for rcs

# get the names of all the variables
allvars <- attr(terms(formula), "term.labels") 

# model for prognostic index, without treatment
model_pi <- model_fit(formula = as.formula(paste0(formula[2], "~", paste0(setdiff(allvars, treat_name), collapse = "+"))),
                      model = model,
                      data = dataset)

# get prognostic index
pi <- scale(predict(model_pi, type = "lp"), center = TRUE, scale = FALSE)
pi_rcs <- as.vector(rcs(pi, 3)[,2])

# constant relative treatment effect
model_rte <- model_fit(formula = formula,
                       model = model, 
                       data = dataset)

# treatment interaction with prognostic index
model_treat_pi <- model_fit(model = model,
                            formula = as.formula(paste0(formula[2], "~", paste0(treat_name, "*", "pi"))),
                            data = dataset)

# treatment interaction with cubic splines
model_treat_pi_rcs <- model_fit(model = model,
                            formula = as.formula(paste0(formula[2], "~", paste0(treat_name, "*", "(pi + pi_rcs)"))),
                            data = dataset)

# all interactions
model_treat_int <- model_fit(model = model, 
                             formula = as.formula(paste0(formula[2], "~", treat_name, " * (",paste0(setdiff(allvars, treat_name), collapse = "+"),")")),
                             data = dataset)

# ridge, all interactions
model_ridge <- model_fit_glmnet(formula = as.formula(paste0(formula[2], "~", treat_name, " * (",paste0(setdiff(allvars, treat_name), collapse = "+"),")")),
                                model = model, 
                             pen = "ridge",
                             treat_name = treat_name,
                             data = dataset)

# lasso, all interactions
model_lasso <- model_fit_glmnet(
                             formula = as.formula(paste0(formula[2], "~", treat_name, " * (",paste0(setdiff(allvars, treat_name), collapse = "+"),")")),
                             model = model, 
                             pen = "lasso",
                             treat_name = treat_name,
                             data = dataset)


## ------------------------------------------------------------------------
models <- list(#pi = model_pi, 
               rte = model_rte, 
               treat_pi = model_treat_pi, 
               treat_pi_rcs = model_treat_pi_rcs, 
               treat_allint = model_treat_int,
               ridge = model_ridge,
               lasso = model_lasso)

# repeat only 3 times cross validation, to be quick here
res <- lapply(models, benefit_cv, folds = 10, times = 3, reproducible = TRUE, treat_name = "treat")

## ------------------------------------------------------------------------
library(pander)
library(egg)
library(broom)

cat_asis <- function(x) knitr::asis_output(capture.output(x))


for(i in 1:length(models)) {
  
  cat_asis(pandoc.header(names(models)[i], level = 2))
  #knitr::asis_outputpandoc.header(names(models)[i], level = 2)
  
  cat_asis(pandoc.header("Estimates", level = 3))
  
  cat("\n")
  
  print(knitr::kable(coef_cbenefit(models[[1]])))
  
  cat("\n")
  
  cat_asis(pandoc.header("C for benefit", level = 3))
  
  cat("\n")
  print(res[[i]])
  
  cat("\n")
  
  cat_asis(pandoc.verbatim(print(res[[i]])))
  
  cat("\n")
  
  cat_asis(pandoc.header("Plots", level = 3))

  ggarrange(plots = calibration_plot(res[[i]]), nrow = 1)
  ggarrange(plots = net_benefit_plot(res[[i]]), nrow = 1)
  
  cat("\n")
}



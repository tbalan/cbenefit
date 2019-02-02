#' Cross validated c index and c for benefit
#'
#' @param model_fit An object as returned by model_fit or model_fit_glmnet
#' @param treat_name The name of the treatment variable as a string
#' @param folds Number of folds for k-fold cross validation
#' @param times Number of repetitions of the cross validation
#' @param reproducible Logical. If `TRUE`, then at each repetition of the cross validation and every time `matchit` is called the seed is fixed.
#' @param match_by One of "benefit" or "covariates" - what should the matching be based on?
#' @param data Must only be specified if match_by = "covariates", the data frame.
#' @param replace Logical, whether matching should be done with replacement. Default is `FALSE`. This argument is passed to MatchIt::matchit and overrides the parameters given in the `matchit_args` argument.
#' @param matchit_args A list with other arguments that are passed on to matchit. See ?MatchIt::matchit for options like distance or method of matching.
#' @return An object of class cbenefit
#' @export
#'
#' @importFrom Hmisc rcorr.cens
#' @importFrom MatchIt matchit
#' @importFrom magrittr "%>%"
#' @importFrom tibble as_tibble
#' @importFrom stats model.frame terms binomial glm lowess model.matrix quantile sd terms update predict as.formula coef
#' @importFrom survival survfit
#' @importFrom utils tail data
#' @importFrom rlang .data
#' @importFrom dplyr n
#' @examples
#' m1 <- model_fit(Surv(time, status) ~  treat + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
#' dat_toy, model = "coxph")
#
#' benefit_cv(m1, "treat", folds = 5, times = 3)

benefit_cv <- function(model_fit, treat_name, folds = 10, times = 10, reproducible = FALSE, match_by = "benefit", data = NULL, replace = FALSE, matchit_args = list()) {

  # 1. Introduction
  if(folds < 2) stop("folds must be > 1")

  if(match_by == "covariates" & is.null(data)) stop("matching by covariates not possible unless you also specify the data frame")

  # try to get the data frame if you pass it as a string
  if(match_by == "covariates" & !inherits(data, "data.frame"))
      data = get(data, parent.frame())

  # threshold for net benefit
  thr <- seq(from = -0.2,
             to = 0.3,
             length.out = 20)

  matchit_args$replace <- replace


  # 2. Calculate predicted benefit
  if(inherits(model_fit, "glm")) {

    treat <- model.frame(model_fit)[[treat_name]]
    mf <- model.frame(model_fit)
    allvars <- attr(terms(model_fit), "term.labels")
    folds_column <- integer(nrow(mf))

    p <- predict(model_fit, type = "response", newdata = mf)

    mf[[treat_name]] <- 1 - mf[[treat_name]]
    p_other <- predict(model_fit, type = "response", newdata = mf)

    pred_ben <- (p - p_other) * (-1)^(1 - mf[[treat_name]])
    mf[[treat_name]] <- 1 - mf[[treat_name]]
  }

  if(inherits(model_fit, "coxph")) {
    # some stuff that will be useful further on..
    treat <- model.frame(model_fit)[[treat_name]]

    # use the mf to re-fit the model
    mf <- model.frame(model_fit)

    mf$time <- mf[[1]][,1]
    mf$status <- mf[[1]][,2]
    mf <- mf[-1]


    allvars <- attr(terms(model_fit), "term.labels")
    folds_column <- integer(nrow(mf))

    # get c for benefit and c
    # model_fit
    H_tau <- tail(survfit(model_fit, se.fit = FALSE)$cumhaz, n = 1)
    # p <- 1 - exp(- predict(model_fit, type = "risk", newdata = mf)  * H_tau)

    p <- exp(- predict(model_fit, type = "risk", newdata = mf)  * H_tau)

    mf[[treat_name]] <- 1 - mf[[treat_name]]
    p_other <- exp(- predict(model_fit, type = "risk",
                             newdata = mf ) * H_tau)
    # model_fit
    mf[[treat_name]] <- 1 - mf[[treat_name]]

    pred_ben <- (p - p_other) * (-1)^(1 - mf[[treat_name]])


  }

  if(inherits(model_fit$glmnetfit, "cv.glmnet")) {
    folds_column <- integer(nrow(model_fit$x))
    treat <- as.numeric(model_fit$x[,colnames(model_fit$x) == treat_name])


    allvars <- colnames(model_fit$x)[!grepl(":", colnames(model_fit$x))]

    # for later - which columns to not penalize
    penalty <- rep(1, ncol(model_fit$x))
    penalty[which(colnames(model_fit$x) == treat_name)] <- 0

    # model_fit



    # for logistic regression:
    if(inherits(model_fit$y, "numeric")) {
      p <- predict(model_fit$glmnetfit, newx = model_fit$x, s = "lambda.min", type = "response")
      p_other <- predict(model_fit$glmnetfit, newx = model_fit$x_other, s = "lambda.min", type = "response")

      pred_ben <- (p - p_other) * (-1)^treat
    }


    # for cox models:
    if(inherits(model_fit$y, "Surv")) {
      lp <- predict(model_fit$glmnetfit, newx = model_fit$x, s = "lambda.min")
      lp_other <- predict(model_fit$glmnetfit, newx = model_fit$x_other, s = "lambda.min")

      mcox <- coxph(model_fit$y ~ offset(lp))
      H_tau <- tail(survfit(mcox, se.fit = FALSE)$cumhaz, n = 1)

      p <- exp(-predict(mcox,
                        type = "risk",
                        newdata = data.frame(lp = as.numeric(lp)) * H_tau))
      p_other <- exp(- predict(mcox, type = "risk",
                               newdata = data.frame(lp = as.numeric(lp_other) ) * H_tau))

      # check the shit out of this
      # like this it's negative but wth
      pred_ben <- (p - p_other) * (-1)^(1 - treat)
    }
  }

  if(isTRUE(reproducible)) set.seed(1)

  # 3. Matchi
  if(abs(table(treat)[2] - table(treat)[1]) > 0 & replace == FALSE) {
    warning(paste0("Matching without replacement, but some pairs are incomplete. Discarding ", abs(table(treat)[2] - table(treat)[1])," row(s) from the data for which a match could not be found.
See the m.order argument in the MatchIt::matchit() function."))
  }


  # 3a. Matching by predicted benefit
  if(match_by == "benefit") {
    matchit_args$formula <- treat~pred_ben
    matchit_args$data <- data.frame(treat = treat, pred_ben = pred_ben, row.names = NULL)
    match <- do.call(function(...) matchit(...)$match.matrix, matchit_args)
  }

  # 3b. Matching by covariates
  if(match_by == "covariates") {
    matchit_args$formula <-  as.formula(paste0(treat_name, "~", paste0(setdiff(allvars, treat_name), collapse = "+")))
    matchit_args$data <- data
    match <- do.call(function(...) matchit(...)$match.matrix, matchit_args)
  }
  # match <- MatchIt::matchit(as.formula(paste0(treat_name, "~", paste0(setdiff(allvars, treat_name), collapse = "+"))), data = data)$match.matrix

  # 3c. Calculate average predicted benefit for each pair
  treated <- as.numeric(rownames(match))
  control <- as.numeric(match[,1])
  predicted <- (pred_ben[treated] + pred_ben[control])/2

  # browser()
  # 3d. Check if there are NAs. This happens if some individuals don't have a match!
  # browser

  # browser()

  if(any(is.na(predicted))) {


    no_pair <- which(is.na(predicted))
    predicted <- predicted[-no_pair]
    # observed <- observed[-no_pair]
    treated <- treated[-no_pair]
    control <- control[-no_pair]

  }

  # browser()
  # 4. Calculate observed benefit
  intervals <- cut(predicted,
                   breaks = unique(quantile(predicted, probs = seq(from = 0, to = 1, length.out = 5))),
                   include.lowest = TRUE,
                   ordered_result = TRUE)

  if(inherits(model_fit$y, "Surv")) {

    # define observed benefit
    observed <- 1 * (model_fit$y[control, 2] == 1 & model_fit$y[control, 1] < model_fit$y[treated, 1]) -
      1 * (model_fit$y[treated, 2] == 1 & model_fit$y[treated, 1] < model_fit$y[control, 1])


    # observed benefit / benefit quartile
    obs_pred_nocv <- split(
      data.frame(treated, control, predicted, observed),
      intervals) %>%
      lapply(function(x) list(treated = survfit(model_fit$y[x$treated,] ~ 1),
                              control = survfit(model_fit$y[x$control,] ~ 1))) %>%
      lapply(function(x) lapply(x, function(sf) c(tail(sf$surv, n = 1),
                                                  tail(sf$upper - sf$lower, n = 1)/(4 * 1.96)))) %>%
      lapply(function(x) c(mean_y = x$treated[1] - x$control[1],
                           sd_y = sqrt(x$treated[2]^2 + x$control[2]^2)))  %>%
      do.call(rbind, .) %>%
      tibble::as_tibble(rownames = "interval")

    pred_benefit_nocv <- split(predicted, intervals) %>%
      lapply(mean) %>%
      do.call(c, .)

    obs_pred_nocv$mean_predicted <- pred_benefit_nocv

    # net benefit

    event_rate_controls <- 1 -
      tail(survfit(model_fit$y[treat == 0] ~ 1)$surv, n = 1)

    event_rate_treated <- 1 -
      tail(survfit(model_fit$y[treat == 1] ~ 1)$surv, n = 1)

    # for net benefit we do not need pairs
    nb <- vapply(thr, function(threshold) {
      congruent <- (treat == 1 & pred_ben >= threshold) |
        (treat == 0 & pred_ben < threshold)

      event_rate_congruents <- 1 -
        tail(survfit(model_fit$y[congruent]~1)$surv, n = 1)

      reduction <- event_rate_controls - event_rate_congruents

      reduction - threshold * mean(treat[congruent])
    }, FUN.VALUE = 0.0)

    nb_nocv <- nb


  }

  if(inherits(model_fit$y, "numeric")) {
    observed <- model_fit$y[control] - model_fit$y[treated]

    obs_pred_cv <- data.frame(treated, control, predicted, observed, interval = intervals) %>%
      dplyr::group_by(.data$interval) %>%
      dplyr::summarize(mean_pred = mean(predicted),
                       mean_y = mean(observed),
                       sd_y = sd(observed) / sqrt(n()))


    event_rate_controls <- mean(model_fit$y[treat == 0])
    event_rate_treated <- mean(model_fit$y[treat == 1])

    nb_nocv<- vapply(thr, function(threshold) {
      congruent <- (treat == 1 & pred_ben >= threshold) |
        (treat == 0 & pred_ben < threshold)

      event_rate_congruents <- mean(model_fit$y[congruent])
      reduction <- event_rate_controls - event_rate_congruents

      reduction - threshold * mean(treat[congruent])

    }, FUN.VALUE = 0.0)

  }

  Sm <- lowess(predicted, observed, iter=0)

  cindex_nocv <- c(c = Hmisc::rcorr.cens(p, model_fit$y)[[1]],
                   c_ben = Hmisc::rcorr.cens(predicted, observed)[[1]],
                   e_mean_ben = mean(abs(Sm$x - Sm$y)),
                   e_90_ben = quantile(abs(Sm$x - Sm$y), .9, names = FALSE))


  # 5. Cross-validation
  pred_risk <- numeric(length(folds_column))
  pred_ben <- numeric(length(folds_column))
  ave_pred_risk <- numeric(length(folds_column))
  ave_pred_ben <- numeric(length(folds_column))



  if(times < 1) {

    # If there is no cross-validation

    predicted_benefit <- predicted
    cindex_cv_ave <- NA
    cindex_cv <- NA
    obs_pred_cv <- NA
    nb_cv <- NA

  } else {

    # If there is cross validation

   #  browser()
    cindex_cv <- matrix(0, times, 4, dimnames = list(c(1:times), c("c", "c_ben", "e_mean_ben", "e_90_ben")))

    for(iter in 1:times) {

      if(isTRUE(reproducible)) set.seed(iter)

      # Split in folds
      folds_treated <- sample(x = 1:folds, size = sum(treat==1), replace = TRUE)
      folds_control <- sample(x = 1:folds, size = sum(treat==0), replace = TRUE)

      folds_column[treat == 1] <- folds_treated
      folds_column[treat == 0] <- folds_control


      for(i in 1:folds) {

        # Predict benefit in each fold

        if(inherits(model_fit, "glm")) {


          mfit <- glm(model_fit$formula, family = "binomial", data = mf[folds_column != i,])

          pred_df <- mf[folds_column == i,]

          p <- predict(mfit, type = "response", newdata = pred_df)

          pred_df[[treat_name]] <- 1 - pred_df[[treat_name]]
          p_other <- predict(mfit, type = "response", newdata = pred_df)

          pred_ben[folds_column == i] <- (p - p_other) * (-1)^(1 - pred_df[[treat_name]])
          pred_risk[folds_column == i] <- p

        }

        if(inherits(model_fit, "coxph")) {
          mph <- coxph(model_fit$formula, model = TRUE, data =  mf[folds_column != i,])
          H_tau <- tail(survfit(mph, se.fit = FALSE)$cumhaz, n = 1)
          pred_df <- mf[folds_column == i,]
          p <- exp(- predict(mph, type = "risk", newdata = pred_df)  * H_tau)

          pred_df[[treat_name]] <- 1 - pred_df[[treat_name]]
          p_other <- exp(- predict(model_fit, type = "risk",
                                   newdata = pred_df ) * H_tau)

          pred_ben[folds_column == i] <- (p - p_other) * (-1)^(pred_df[[treat_name]])
          pred_risk[folds_column == i] <- p
        }

        if(inherits(model_fit$glmnetfit, "cv.glmnet")) {
          # no penalization for the treatment main effect
          penalty <- rep(1, ncol(model_fit$x))
          penalty[which(colnames(model_fit$x) == treat_name)] <- 0




          if(inherits(model_fit$y, "numeric")) {
            mglmnet <- cv.glmnet(model_fit$x[folds_column != i,],
                                 model_fit$y[folds_column != i],
                                 family = "binomial", alpha = model_fit$alpha, penalty.factor = penalty)

            p <- predict(mglmnet, newx = model_fit$x[folds_column == i, ], s = "lambda.min", type = "response")
            p_other <- predict(mglmnet, newx = model_fit$x_other[folds_column == i, ], s = "lambda.min", type = "response")


            pred_ben[folds_column == i] <- (p - p_other) * (-1)^(treat[folds_column == i])
            pred_risk[folds_column == i] <- p
          }

          if(inherits(model_fit$y, "Surv")) {
            # HERE and if ridge then alpha = 0 blabla
            mglmnet <- cv.glmnet(model_fit$x[folds_column != i,],
                                 model_fit$y[folds_column != i],
                                 family = "cox", alpha = model_fit$alpha, penalty.factor = penalty)

            lp <- predict(mglmnet, newx = model_fit$x[folds_column != i,], s = "lambda.min")
            mcox <- coxph(model_fit$y[folds_column != i] ~ offset(lp))
            H_tau <- tail(survfit(mcox, se.fit = FALSE)$cumhaz, n = 1)

            lp_rest <-
              predict(mglmnet, newx = model_fit$x[folds_column == i, ], s = "lambda.min")
            lp_other_rest <-
              predict(mglmnet, newx = model_fit$x_other[folds_column == i, ], s = "lambda.min")

            p <- exp(-predict(mcox,
                              type = "risk",
                              newdata = data.frame(lp = as.numeric(lp_rest))) * H_tau)
            p_other <- exp(- predict(mcox, type = "risk",
                                     newdata = data.frame(lp = as.numeric(lp_other_rest)) ) * H_tau)


            pred_ben[folds_column == i] <- (p - p_other) * (-1)^(1 - treat[folds_column == i])
            pred_risk[folds_column == i] <- p
          }
        }


      }

      # Re-match the pairs
      if(isTRUE(reproducible)) set.seed(iter)

      if(match_by == "benefit") {
        # matchit_args$formula <- treat~pred_ben
        matchit_args$data <- data.frame(treat = treat, pred_ben = pred_ben, row.names = NULL)
        match <- suppressWarnings(do.call(function(...) matchit(...)$match.matrix, matchit_args))
      }

      # if(match_by == "benefit")
      #   match <- MatchIt::matchit(treat ~ pred_ben, data = data.frame(treat, pred_ben))$match.matrix

      treated <- as.numeric(rownames(match))
      control <- as.numeric(match[,1])

      predicted <- (pred_ben[treated] + pred_ben[control])/2

      # Again remove pairs with NA for predicted, not warning again
      if(any(is.na(predicted))) {
        # warning(paste0("Removed ", sum(is.na(predicted)), " row(s) from the data for which a match could not be found"))
        no_pair <- which(is.na(predicted))
        predicted <- predicted[-no_pair]
        # observed <- observed[-no_pair]
        treated <- treated[-no_pair]
        control <- control[-no_pair]

      }


      if(inherits(model_fit$y, "Surv")) {
        observed <- 1 * (model_fit$y[control, 2] == 1 & model_fit$y[control, 1] < model_fit$y[treated, 1]) -
          1 * (model_fit$y[treated, 2] == 1 & model_fit$y[treated, 1] < model_fit$y[control, 1])
      }

      if(inherits(model_fit$y, "numeric")) {
        observed <- model_fit$y[control] - model_fit$y[treated]
      }




      Sm <- lowess(predicted, observed, iter=0)

      cindex <- c(c = rcorr.cens(pred_risk, model_fit$y)[[1]],
                  c_ben = rcorr.cens(predicted, observed)[[1]],
                  e_mean_ben = mean(abs(Sm$x - Sm$y)),
                  e_90_ben = quantile(abs(Sm$x - Sm$y), .9, names = FALSE))

      # cindex <- c(c = rcorr.cens(-pred_risk, model_fit$y)[[1]],
      #                  c_benefit = rcorr.cens(predicted, observed)[[1]],
      #                  e_mean_benefit = mean(abs(Sm$x - Sm$y)),
      #                  e_90_benefit = quantile(abs(Sm$x - Sm$y), .9, names = FALSE))
      cindex_cv[iter, ] <- cindex

      # again for this stuff we don't really need the pairs

      ave_pred_ben <- ave_pred_ben + pred_ben / times
      ave_pred_risk <- ave_pred_risk + pred_risk / times

    }



  if(isTRUE(reproducible)) set.seed(1)
    if(match_by == "benefit") {
      matchit_args$formula <- treat~ave_pred_ben
      matchit_args$data <- data.frame(treat = treat, ave_pred_ben = ave_pred_ben, row.names = NULL)
      match <- suppressWarnings(do.call(function(...) matchit(...)$match.matrix, matchit_args))
    }

    # if(match_by == "benefit")
    #  match <- MatchIt::matchit(treat ~ ave_pred_ben, data = data.frame(treat, ave_pred_ben))$match.matrix

  treated <- as.numeric(rownames(match))
  control <- as.numeric(match[,1])


  predicted_ave <- (ave_pred_ben[treated] + ave_pred_ben[control])/2

  # View(data.frame(treated, control))
  if(any(is.na(predicted_ave))) {
    # warning(paste0("Removed ", sum(is.na(predicted)), " row(s) from the data for which a match could not be found"))
    no_pair <- which(is.na(predicted_ave))
    predicted_ave <- predicted_ave[-no_pair]

    # observed <- observed[-no_pair]
    treated <- treated[-no_pair]
    control <- control[-no_pair]

  }

  intervals <- cut(predicted_ave,
                   breaks = unique(quantile(predicted, probs = seq(from = 0, to = 1, length.out = 5))),
                   include.lowest = TRUE,
                   ordered_result = TRUE)


  if(inherits(model_fit$y, "Surv")) {
    observed <- 1 * (model_fit$y[control, 2] == 1 & model_fit$y[control, 1] < model_fit$y[treated, 1]) -
      1 * (model_fit$y[treated, 2] == 1 & model_fit$y[treated, 1] < model_fit$y[control, 1])



    obs_pred_cv <- split(
      data.frame(treated, control, predicted, observed),
      intervals) %>%
      lapply(function(x) list(treated = survfit(model_fit$y[x$treated,] ~ 1),
                              control = survfit(model_fit$y[x$control,] ~ 1))) %>%
      lapply(function(x) lapply(x, function(sf) c(tail(sf$surv, n = 1),
                                                  tail(sf$upper - sf$lower, n = 1)/(4 * 1.96)))) %>%
      lapply(function(x) c(mean_y = x$treated[1] - x$control[1],
                           sd_y = sqrt(x$treated[2]^2 + x$control[2]^2)))  %>%
      do.call(rbind, .) %>%
      tibble::as_tibble(rownames = "interval")

    pred_benefit_cv <- split(predicted, intervals) %>%
      lapply(mean) %>%
      do.call(c, .)

    obs_pred_cv$mean_predicted <- pred_benefit_cv


    nb_cv <- vapply(thr, function(threshold) {
      congruent <- (treat == 1 & ave_pred_ben >= threshold) |
        (treat == 0 & ave_pred_ben < threshold)

      event_rate_congruents <- 1 -
        tail(survfit(model_fit$y[congruent]~1)$surv, n = 1)

      reduction <- event_rate_controls - event_rate_congruents

      reduction - threshold * mean(treat[congruent])
    }, FUN.VALUE = 0.0)




  }

  if(inherits(model_fit$y, "numeric")) {
    observed <- model_fit$y[control] - model_fit$y[treated]


    obs_pred_nocv <- data.frame(treated, control, predicted, observed, interval = intervals) %>%
      dplyr::group_by(.data$interval) %>%
      dplyr::summarize(mean_pred = mean(predicted),
                       mean_y = mean(observed),
                       sd_y = sd(observed) / sqrt(n()))

    nb_cv<- vapply(thr, function(threshold) {
      congruent <- (treat == 1 & ave_pred_ben >= threshold) |
        (treat == 0 & ave_pred_ben < threshold)

      event_rate_congruents <- mean(model_fit$y[congruent])
      reduction <- event_rate_controls - event_rate_congruents

      reduction - threshold * mean(treat[congruent])

    }, FUN.VALUE = 0.0)
    }

  Sm <- lowess(predicted, observed, iter=0)
  # browser()
  #
  # ave_pred_risk
  #

  cindex_cv_ave <- c(c = rcorr.cens(ave_pred_risk, model_fit$y)[[1]],
                     c_ben = rcorr.cens(predicted, observed)[[1]],
                     e_mean_ben = mean(abs(Sm$x - Sm$y)),
                     e_90_ben = quantile(abs(Sm$x - Sm$y), .9, names = FALSE))



}

  # browser()


  # browser()


  res <- list(obs_pred_nocv = obs_pred_nocv,
              obs_pred_cv = obs_pred_cv,
              nb_nocv = nb_nocv,
              nb_cv = nb_cv,
              threshold = thr,
              event_rate_controls = event_rate_controls,
              event_rate_treated = event_rate_treated,
              cindex_nocv = cindex_nocv,
              cindex_cv_ave = cindex_cv_ave,
              c_index_cv = cindex_cv,
              call = match.call())

  attr(res, "class") <- "cbenefit"


  return(res)

}


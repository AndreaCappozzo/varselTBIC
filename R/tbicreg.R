# Computes BIC for the regression of y on x, after doing subset selection
# on the predictors
# require(BMA)

BICreg <- function(mod)
{
  sigma <- sqrt(mean(residuals(mod) ^ 2))
  n <- length(mod$residuals)
  p <- n - df.residual(mod) + 1
  # calculate the BIC for the regression
  -n * log(2 * pi) - 2 * n * log(sigma) - n - log(n) * p
}

fit_reg <- function(x, y)
{
  x <- as.matrix(x)
  y <- as.vector(y)
  mod <- bicreg(y = y, x = x, nbest = 1)
  subset <- which(mod$which[1, ])
  lm.fit(y = y, x = cbind(1, x[, subset, drop = FALSE]))
}

TBIC_computation <- function(X_p,
                             X_c,
                             cltrain,
                             union_trimmed_units,
                             fitted_model,
                             model_type) {
  Xtrain <- switch (model_type,
                    "GR" = cbind(X_c, X_p),
                    "NG" = X_c)

  untrimmed_units <-
    setdiff(1:fitted_model$Ntrain, union_trimmed_units)

  # Mixture part

  component_density <-
    cdens(
      modelName = fitted_model$modelName,
      data = Xtrain[untrimmed_units, , drop = FALSE],
      logarithm = TRUE,
      parameters = list(
        mean = fitted_model$parameters$mean,
        variance = fitted_model$parameters$variance
      )
    )

  log_lik_mixture <-
    sum(
      mclust::unmap(cltrain[untrimmed_units]) * sweep(
        x = component_density,
        MARGIN = 2,
        STATS = log(fitted_model$parameters$pro),
        FUN = "+"
      )
    )

  tbic_mixture <-
    2 * log_lik_mixture - mclust::nMclustParams(modelName = fitted_model$modelName,
                                                d = fitted_model$d,
                                                G = fitted_model$G) * log(length(untrimmed_units))
  if (model_type == "GR") {
    return(tbic_mixture)
  } else {
    ll_reg <-
      sum(
        dnorm(
          X_p[untrimmed_units],
          mean = fitted_model$parameters_regression$alpha + Xtrain[untrimmed_units, names(fitted_model$parameters_regression$beta), drop = FALSE] %*% fitted_model$parameters_regression$beta,
          sd = fitted_model$parameters_regression$sigma,
          log = T
        )
      )
    tbic_regression <-
      2 * ll_reg - fitted_model$parameters_regression$number_of_param * log(length(untrimmed_units))
    tbic_mixture + tbic_regression
  }
}


# General function for computing BIC once models are robustly estimated --------


BIC_computation <- function(X_p,
                             X_c,
                             cltrain,
                             fitted_model,
                             model_type) {
  Xtrain <- switch (model_type,
                    "GR" = cbind(X_c, X_p),
                    "NG" = X_c,
                    "only_mixture" = X_c)
  
  # Mixture part
  
  component_density <-
    cdens(
      modelName = fitted_model$modelName,
      data = Xtrain,
      logarithm = TRUE,
      parameters = list(
        mean = fitted_model$parameters$mean,
        variance = fitted_model$parameters$variance
      )
    )
  
  log_lik_mixture <-
    sum(
      mclust::unmap(cltrain) * sweep(
        x = component_density,
        MARGIN = 2,
        STATS = log(fitted_model$parameters$pro),
        FUN = "+"
      )
    )
  
  bic_mixture <-
    2 * log_lik_mixture - mclust::nMclustParams(modelName = fitted_model$modelName,
                                                d = fitted_model$d,
                                                G = fitted_model$G) * log(nrow(Xtrain))
  if (model_type == "GR"| model_type == "only_mixture" ) {
    return(bic_mixture)
  } else {
    ll_reg <-
      sum(
        dnorm(
          X_p,
          mean = fitted_model$parameters_regression$alpha + Xtrain[, names(fitted_model$parameters_regression$beta), drop = FALSE] %*% fitted_model$parameters_regression$beta,
          sd = fitted_model$parameters_regression$sigma,
          log = T
        )
      )
    bic_regression <-
      2 * ll_reg - fitted_model$parameters_regression$number_of_param * log(nrow(Xtrain))
    bic_mixture + bic_regression
  }
}

redda_model_varsel <- function(X_p,
                               #the variable (univariate) proposed for inclusion
                               X_c,
                               # the variables currently included in the model
                               cltrain,
                               alpha_Xtrain = 0,
                               # The proportion of obs to be trimmed
                               modelName = "EEE",
                               nsamp = 500,
                               # the number of initial samples of the patterned MCD-like procedure
                               model = c("GR", "NG"),
                               #Grouping or the No Grouping model to be estimated?
                               ...) {
  X_p <-
    if (is.null(X_p)) {
      NULL
    } else {
      data.matrix(X_p)
    } #used just for the single variable group model
  X_c <- data.matrix(X_c)
  Xtrain <- switch (model,
                    "GR" = cbind(X_c, X_p),
                    "NG" = X_c)
  Ntrain <- nrow(Xtrain)
  ltrain <- mclust::unmap(cltrain)
  d <- ncol(Xtrain)
  G <- ncol(ltrain)
  cltrain <- as.factor(cltrain)
  classLabel <- levels(cltrain)

  # Core algorithm ----------------------------------------------------------


  if (alpha_Xtrain != 0) {
    Ntrain_trim <- Ntrain - ceiling(Ntrain * alpha_Xtrain)
    robust_result <-
      patterned_MCD_varsel(nsamp = nsamp,
                    # I perform the estimation starting from nsamp J_g subsets of sample size (d+1), inspired to what done in \cite{Hubert2018}
                    Xtrain,
                    X_p,
                    #the variable (univariate) proposed for inclusion
                    cltrain,
                    ltrain,
                    alpha_Xtrain,
                    modelName,
                    model)
    ll <- robust_result$ll
    fitm <- robust_result$fitm
    if (model == "NG") {
      ll <-
        ll - fitm$ll_reg # this is the likelihood of the normal density on Xtrain, that is X_c
      lm_fit <- fitm$lm_fit
      alpha <- lm_fit$coefficients[1]
      beta <- lm_fit$coefficients[-1]
      sigma <- sqrt(mean(residuals(lm_fit) ^ 2))
    }
  } else {
    Ntrain_trim <- NULL
    fitm <- mclust::mstep(modelName = modelName,
                          data = Xtrain,
                          z = ltrain)
    mTau.train <-
      matrix(log(fitm$parameters$pro),
             nrow(Xtrain),
             fitm$G,
             byrow = TRUE)
    lDensity.train <- do.call(mclust::cdens, c(list(
      data = Xtrain,
      logarithm = TRUE
    ), fitm))

    sum.train <- mTau.train + lDensity.train
    mat.train <- ltrain * sum.train
    ll <- sum(mat.train)
    if (model == "NG") {
      lm_fit <-
        fit_reg(x = Xtrain, y = X_p) #linear regression with stepwise variable selection
      alpha <- lm_fit$coefficients[1]
      beta <- lm_fit$coefficients[-1]
      sigma <- sqrt(mean(residuals(lm_fit) ^ 2))
      ll_reg <-
        sum(dnorm(
          X_p,
          mean = alpha + Xtrain[, names(beta), drop = FALSE] %*% beta,
          sd = sigma,
          log = T
        ))
    }
  }


  # Checking if errors in the procedure -------------------------------------

  fitetrain <-
    tryCatch(
      do.call(mclust::estep, c(list(data = Xtrain), fitm)),
      error = function(e) {
        list(z = NA)
      }
    )
  emptyz <- ifelse(all(!is.na(fitetrain$z)), yes = FALSE, no = TRUE)

  # Results Collection ------------------------------------------------------

  if (!emptyz) {
    res <- list()
    res$Ntrain <- Ntrain
    res$Ntrain_after_trimming <- Ntrain_trim
    res$alpha_Xtrain <- alpha_Xtrain
    res$d <- d
    res$G <- G
    res$modelName <- modelName
    res$parameters <- fitm$parameters
    #Name the different groups
    names(res$parameters$pro) <- classLabel
    if (d == 1) {
      names(res$parameters$mean) <- classLabel
    } else {
      colnames(res$parameters$mean) <- classLabel
    }
    if (model == "NG") {
      res$parameters_regression <-
        list(alpha = alpha,
             beta = beta,
             sigma = sigma,
             number_of_param = length(lm_fit$residuals) - df.residual(lm_fit) + 1)
    }
    ztrain <- fitetrain$z
    cltrain <-
      factor(sapply(map(ztrain), function(i)
        classLabel[i]), levels = classLabel) # I classify a posteriori also the trimmed units
    pos_trimmed_train <- NULL
    cltrain_after_trimming <- NULL
    if (alpha_Xtrain != 0) {
      D_Xtrain_cond <- do.call(mclust::cdens, c(list(
        data = Xtrain, # computing the component density, this is done because I am interested in trimming also obs that might have
        logarithm = T
      ), fitm)) # been WRONGLY assinged to a class
      ind_D_Xtrain_cdens <-
        cbind(1:Ntrain, mclust::map(ltrain)) # matrix with obs and respective group from ltrain
      D_Xtrain <-
        D_Xtrain_cond[ind_D_Xtrain_cdens] # I compute D_g conditioning ON the fact that I know the supposed true class
      if (model == "NG") {
        D_reg <-
          as.vector(dnorm(
            X_p,
            mean = alpha + Xtrain[, names(beta), drop = FALSE] %*% beta,
            sd = sigma,
            log = T
          ))
        D_Xtrain <- D_Xtrain + D_reg #D_{No Grouping}
      }
      pos_trimmed_train <-
        # I trim the conditional density \phi(x_n; \mu_g, \Sigma_g) when x_n comes from group g
        which(D_Xtrain <= (sort(D_Xtrain, decreasing = F)
                           [[ceiling(Ntrain * alpha_Xtrain)]]))
      cltrain_after_trimming <- cltrain
      cltrain_after_trimming <-
        factor(cltrain, levels = c(classLabel, "0"))
      cltrain_after_trimming[pos_trimmed_train] <- "0"
    }

    res$train <- list()
    res$train$z <- ztrain
    res$train$cl <- cltrain
    res$train$cl_after_trimming <- cltrain_after_trimming
    res$train$obs_trimmed <- pos_trimmed_train
    res$train$alpha_Xtrain <- alpha_Xtrain
    bic_noreg <-
      2 * ll - mclust::nMclustParams(
        modelName = modelName,
        d = d,
        G = G,
        noise = FALSE,
        equalPro = FALSE
      ) * log(fitm$n)
    # this is the TBIC on Xtrain (if model == "NG" it is X_c), without considering the regression contribution
    if (model == "NG") {
      bic.all <- bic_noreg + BICreg(lm_fit)
    } else {
      bic.all <- bic_noreg
    }
    res$ll <- ll
    res$bic <- bic.all
    res$bic_noreg <- bic_noreg
  }
  else {
    res <- list()
    res$error <-
      "Either training groups too small or trimming level too large for this model type"
    res$ll <- NA
    res$bic <- NA
    res$bic_noreg <- NA
  }
  res
}

# Function used in random initialization ---------------------------------


patterned_MCD <-
  # main function that will be repeated nsamp times, run only if alpha_Xtrain!=0
  function(nsamp,
           Xtrain,
           X_p,
           cltrain,
           ltrain,
           alpha_Xtrain,
           modelName, model) {

    robust_start <-
      function(Xtrain,
               X_p,
               cltrain,
               ltrain,
               alpha_Xtrain,
               modelName, model) {
        Ntrain <- nrow(Xtrain)
        J_ind <-
          tryCatch(
            c(sapply(
              split(seq(cltrain), cltrain),
              sample,
              size = ncol(Xtrain) + 1,
              replace = FALSE
            )),
            error = function(e)
              seq(cltrain)
          )

        #Init with the random (p + 1)-subset J_g
        fitm <- tryCatch(
          mclust::mstep(
            modelName = modelName,
            data = Xtrain[J_ind,,drop=F],
            z = ltrain[J_ind,,drop=F]
          ),
          error = function(e) {
            list(parameters = NA)
          }
        )

        if (model == "NG") {
          lm_fit <- fit_reg(x=Xtrain[J_ind,,drop=F], y=X_p[J_ind,,drop=F])
        }
        emptyz <- TRUE
        iter <- 0
        pos_old <-
          1:ceiling(Ntrain * alpha_Xtrain)  # used in the while loop
        # initial random value just using the first ceiling(Ntrain * alpha_Xtrain) obs
        # the estimation will continue running until for two consecutive iterations exactly
        # the same obs are trimmed
        criterion <- TRUE
        while (criterion) {
          iter <- iter + 1
          # in performing robust estimation I discard those obs with smallest value of
          # component density given their label
          fitetrain <-
            tryCatch(
              do.call(mclust::estep, c(list(data = Xtrain), fitm)),
              error = function(e) {
                list(z = NA)
              }
            )
          Xtrain_fit <- NULL # I initialize its value for preventing erros when the EM algorithm fails
          emptyz <- TRUE
          if (all(!is.na(fitetrain$z))) {
            emptyz <- FALSE
            D_Xtrain_cond <-
              do.call(mclust::cdens, c(list(data = Xtrain, # computing the component density, this is done because I am interested in trimming also obs that might have
                                            logarithm = T), fitm)) # been WRONGLY assinged to a class
            ind_D_Xtrain_cdens <-
              cbind(1:Ntrain, mclust::map(ltrain)) # matrix with obs and respective group from ltrain
            D_Xtrain <-
              D_Xtrain_cond[ind_D_Xtrain_cdens] #if model==GR Xtrain=cbind(X_c,X_p) so D_Xtrain is D_grouping
            if (model == "NG") {
              alpha <- lm_fit$coefficients[1]
              beta <- lm_fit$coefficients[-1]
              sigma <- sqrt(mean(residuals(lm_fit) ^ 2))
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
              which(D_Xtrain <= (sort(D_Xtrain, decreasing = F)
                                 [[ceiling(Ntrain * alpha_Xtrain)]]))
            if (length(pos_trimmed_train) != ceiling(Ntrain * alpha_Xtrain)) {
              #condition due to the fact that for some models (the simplest ones usually) it might happen that 2 obs have exactly the same comp prob, leading to errors
              pos_trimmed_train <-
                pos_trimmed_train[1:ceiling(Ntrain * alpha_Xtrain)]
            }
            ltrain_fit <- ltrain[-pos_trimmed_train, , drop = F]
            Xtrain_fit <- Xtrain[-pos_trimmed_train, , drop = F]
            fitm <-
              tryCatch(
                mclust::mstep(
                  modelName = modelName,
                  data = Xtrain_fit,
                  z = ltrain_fit
                ),
                error = function(e) {
                  list(parameters = NA)
                }
              )
            if (model == "NG") {
              X_p_fit <- X_p[-pos_trimmed_train, , drop =F]
              lm_fit <-
                fit_reg(x = Xtrain_fit, y = X_p_fit)
            }
            if (all(pos_old == pos_trimmed_train)) {
              criterion <- FALSE
            } else {
              pos_old <- pos_trimmed_train
            }
          }
          else {
            criterion <- FALSE
          }
        }
        if (!emptyz) {
          mTau.train <-
            matrix(log(fitm$parameters$pro),
                   nrow(Xtrain_fit),
                   fitm$G,
                   byrow = TRUE)
          lDensity.train <- do.call(mclust::cdens, c(list(
            data = Xtrain_fit,
            logarithm = TRUE
          ), fitm))

          sum.train <- mTau.train + lDensity.train
          mat.train <- ltrain_fit * sum.train
          ll <- sum(mat.train)
          if (model == "NG") {
            alpha <- lm_fit$coefficients[1]
            beta <- lm_fit$coefficients[-1]
            sigma <- sqrt(mean(residuals(lm_fit) ^ 2))
            ll_reg <-
              sum(dnorm(
                X_p_fit,
                mean = alpha + Xtrain_fit[, names(beta), drop = FALSE] %*% beta,
                sd = sigma,
                log = T
              ))
            ll <- ll + ll_reg
            #Adding the \hat{\boldsymbol{\beta}^{(0)}} and \hat{\boldsymbol{\sigma}^{(0)}} at the returned estimates
            fitm$lm_fit <- lm_fit
            fitm$ll_reg <- ll_reg
          }
        } else {
          ll <- NA
        }
        fitm$X <- Xtrain_fit # this is used in the initialization of raedda_transductive_model: for VVE and EVE models for which the constraint is not satisfied the data from which the M-step is computed are needed for the MM algorithm
        return(list(ll = ll, fitm = fitm))
      }
    n_init <-
      replicate(
        nsamp,
        expr = robust_start(Xtrain,
                            X_p,
                            cltrain,
                            ltrain,
                            alpha_Xtrain,
                            modelName, model),
        simplify = T
      )
    if(!all(is.na(n_init[1,]))){
      ind <- which.max(n_init[1, ])
    } else {
      ind <- 1 #if no initial subsample worked, get the first one and then the alg will break down further in the code
    }
    return(n_init[, ind])

  }

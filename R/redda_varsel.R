redda_varsel <-
  function(X_p,
           X_c,
           cltrain,
           alpha_Xtrain = 0,
           modelNames = NULL,
           nsamp = 500,
           model,
           ...) {
    X_p <- if(is.null(X_p)) {NULL} else {data.matrix(X_p)} #used just for the single variable group model
    X_c <- data.matrix(X_c)
    Xtrain <- switch (model,
                      "GR" = cbind(X_p, X_c),
                      "NG" = X_c

    )
    Xtrain <- data.matrix(Xtrain)
    if (is.null(modelNames)) {
      if(ncol(Xtrain)==1){
        modelNames <- c("E", "V")
      } else {
        modelNames <- mclust::mclust.options("emModelNames")
      }
    }
    RES <- list()
    bestBIC <- -Inf
    RES[["Best"]] <- list()
    for (modelName in modelNames) {
      RES[[modelName]] <- list()
      RES[[modelName]] <- redda_model_varsel(X_p,
                                             X_c,
                                           cltrain,
                                           alpha_Xtrain,
                                           modelName,
                                           nsamp,
                                           model,
                                           ...)
      if (!is.na(RES[[modelName]]$bic)) {
        if (RES[[modelName]]$bic > bestBIC) {
          RES[["Best"]] <- RES[[modelName]]
          bestBIC <- RES[[modelName]]$bic
        }
      }
    }
    RES
  }

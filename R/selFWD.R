#############################################################################
## Sequential & parallel robust forward greedy search REDDA
#############################################################################

redda_varsel_gr_fwd_TBIC_no_union <- function(Xtrain,
                                cltrain,
                          emModels1 = c("E","V"),
                          emModels2 = mclust.options("emModelNames"),
                          alpha_Xtrain,
                          forcetwo = TRUE,
                          BIC.diff = 0,
                          itermax = 100,
                          nsamp= 100,
                          parallel = FALSE,
                          fit= TRUE,
                          verbose = interactive())
{
  Xtrain <- data.matrix(Xtrain)
  Ntrain <- nrow(Xtrain) # number of rows = number of observations
  Ntrain_trim <- if(alpha_Xtrain!=0) {Ntrain - ceiling(Ntrain * alpha_Xtrain)} else {NULL}
  d <- ncol(Xtrain) # number of columns = number of variables
  if(is.null(colnames(Xtrain))|any(duplicated(colnames(Xtrain)))){
    colnames(Xtrain) <- paste("X", 1:d, sep = "")
  }
  G <- length(unique(cltrain))

  # Start "parallel backend" if needed
  if(is.logical(parallel)) {
    if (parallel) {
      parallel <- clustvarsel::startParallel(parallel)
      stopCluster <-
        TRUE
    } else {
      parallel <- stopCluster <- FALSE
    }
  } else {
    stopCluster <- if (inherits(parallel, "cluster"))
      FALSE
    else
      TRUE
    parallel <- clustvarsel::startParallel(parallel)
  }
  on.exit(if (parallel & stopCluster)
    parallel::stopCluster(attr(parallel, "cluster")))

  # define operator to use depending on parallel being TRUE or FALSE
  `%DO%` <- if(parallel) `%dopar%` else `%do%`
  i <- NULL # dummy to trick R CMD check

  # First Step - selecting single variable
  if(verbose) cat(paste("iter 1\n+ adding step\n"))
  out <- foreach(i = 1:d) %DO% {
    xBIC <- NULL
    # Fit the single variable grouping models
    try(xBIC <- redda_varsel(
      X_p = NULL,
      X_c = Xtrain[, i],
      cltrain = cltrain,
      alpha_Xtrain = alpha_Xtrain,
      modelNames = emModels1,
      nsamp = nsamp,
      model = "GR"
    )$Best,
    silent = TRUE)

    # maxBIC is the maximum BIC over all grouping models fit
    maxBIC <- if(is.null(xBIC$bic)) NA else xBIC$bic
    # Fit and get BIC for a single component no-group normal model
    if (alpha_Xtrain != 0) {
      best_subset <-
        MASS::cov.mcd(x = Xtrain[, i], quantile.used = floor((1 - alpha_Xtrain) * Ntrain))$best
      # I recover the best subset on which the mean and the variance are computed.
      # Afterwards, I perform a robust single component no-grouping normal model using Mclust,
      # this is done because cov.mcd (as well as covMcd in robustbase)
      # use consistency correction that I do not want to include

      try(oneBIC <-
            mclust::Mclust(
              Xtrain[best_subset, i],
              G = 1,
              modelNames = emModels1,
              verbose = FALSE
            )$BIC[1],
          silent = TRUE)
    } else {
      try(oneBIC <- mclust::Mclust(
        Xtrain[, i],
        G = 1,
        modelNames = emModels1,
        verbose = FALSE
      )$BIC[1],
      silent = TRUE)
    }
    return(list(maxBIC, oneBIC, xBIC$modelName, xBIC$G))
  }
  maxBIC <- sapply(out, "[[", 1)
  oneBIC <- sapply(out, "[[", 2)
  # Difference between maximum BIC for grouping and BIC for no grouping
  maxdiff <- maxBIC - oneBIC
  # Find the single variable with the biggest difference between
  # grouping and no grouping
  m <- max(maxdiff[is.finite(maxdiff)])
  arg <- which(maxdiff==m,arr.ind=TRUE)[1]
  # This is our first selected variable/S is the matrix of currently selected
  # grouping variables
  S <- Xtrain[,arg,drop=FALSE]
  # BICS is the BIC value for the grouping model with the variable(s) in S
  BICS <- maxBIC[arg]
  # NS is the matrix of currently not selected variables
  NS <- Xtrain[,-arg,drop=FALSE]
  # info records the proposed variable, BIC for the S matrix and difference
  # in BIC for grouping versus no grouping on S, whether it was an
  # addition step and if it was accepted
  info <- data.frame(Var = colnames(S),
                     BIC = BICS, BICdiff = maxdiff[arg],
                     Step = "Add", Decision = "Accepted",
                     Model = out[[arg]][[3]],
                     G = out[[arg]][[4]],
                     stringsAsFactors = FALSE)

  if(verbose)
    { print(info[,c(1,3:5),drop=FALSE])
      cat("iter 2\n+ adding step\n") }

  # Second Step - selecting second variable
  out <- foreach(i = 1:ncol(NS)) %DO%
  {
    # cindepBIC is the BIC for the grouping model on S and the
    # regression model of the new variable on S
    try(cindepBIC <- redda_varsel(
      X_p = NS[,i],
      X_c = S,
      cltrain = cltrain,
      alpha_Xtrain = alpha_Xtrain,
      modelNames = emModels1,
      nsamp = nsamp,
      model = "NG"
    )$Best,
    silent = TRUE)
    cindepBIC <- if(is.null(cindepBIC$bic)) NA else cindepBIC$bic
    # Fit the cluster model on the two variables
    sBIC <- NULL
    try(sBIC <- redda_varsel(
      X_p = NS[,i],
      X_c = S,
      cltrain = cltrain,
      alpha_Xtrain = alpha_Xtrain,
      modelNames = emModels2,
      nsamp = nsamp,
      model = "GR"
    )$Best,
        silent = TRUE)
    # depBIC is the BIC for the grouping model with both variables
    depBIC <- if(is.null(sBIC$bic)) NA else sBIC$bic

    return(list(cindepBIC, depBIC, sBIC$modelName, sBIC$G))
  }
  cindepBIC <- sapply(out, "[[", 1)
  depBIC <- sapply(out, "[[", 2)
  # cdiff is the difference between BIC for the models with variables' being
  # grouping variables versus them being conditionally independent of
  # the grouping
  cdiff <- depBIC - cindepBIC
  # Choose the variable with the largest difference
  m <- max(cdiff[is.finite(cdiff)])
  arg <- which(cdiff==m,arr.ind=TRUE)[1]

  # if forcetwo is true automatically add the best second variable,
  # otherwise only add it if its difference is large than BIC.diff
  if(forcetwo || cdiff[arg] > BIC.diff) {
    k <- c(colnames(S), colnames(NS)[arg])
    nks <- c(colnames(NS)[-arg])
    BICS <- depBIC[arg]
    info <- rbind(info, c(colnames(NS)[arg], BICS, cdiff[arg],
                          "Add", "Accepted",
                          out[[arg]][3:4]))

    S <- cbind(S, NS[, arg])
    NS <- as.matrix(NS[, -arg])
    colnames(S) <- k
    colnames(NS) <- nks
  } else {
    info <- rbind(info, c(colnames(NS)[arg], BICS, cdiff[arg],
                          "Add", "Rejected",
                          out[[arg]][3:4]))
  }
  info$BIC <- as.numeric(info$BIC)
  info$BICdiff <- as.numeric(info$BICdiff)

  if (verbose){
    print(info[2, c(1, 3:5), drop = FALSE])
  }
  criterion <- 1
  iter <- 2
  while((criterion == 1) & (iter < itermax))
  {
    iter <- iter+1
    check1 <- colnames(S)

    if(verbose) cat(paste("iter", iter, "\n"))

    # Adding step
    if(verbose) cat("+ adding step\n")
    # For the special case where we have removed all the grouping
    # variables/S is empty
    if(ncol(S)==0 || is.null(ncol(S))){
        # We simply choose the same variable as in the first step and check
        # whether the difference between the BIC for grouping versus not
        # grouping is positive or not
        m <- max(maxdiff[is.finite(maxdiff)])
        arg <- which(maxdiff==m,arr.ind=TRUE)[1]
        if(maxdiff[arg] > BIC.diff)
          {
            # if the difference is positive this variable is selected
            # as a grouping variable
            S <- matrix(c(Xtrain[,arg]), Ntrain, 1)
            BICS <- maxBIC[arg]
            colnames(S) <- colnames(Xtrain)[arg]
            NS <- as.matrix(Xtrain[,-arg])
            colnames(NS) <- colnames(Xtrain)[-arg]
            info <- rbind(info, c(colnames(S), BICS, maxdiff[arg],
                                  "Add","Accepted", NA, NA))
          } else { # if the difference is not > BIC.diff no grouping variables exist
            BICS <- NA
            info <- rbind(info,
                          c(colnames(Xtrain)[arg], BICS, maxdiff[arg],
                            "Add", "Rejected", NA, NA))
          }
      } else { # Addition Step in general (for all cases except when S is empty)
        if(ncol(NS) != 0 & !is.null(ncol(NS))) {
          
          if(ncol(S)==1) { 
            modelNames_NG <- emModels1 # if S is one-dimensional, univariate models need to be fitted in the No-Grouping structure
          } else{
            modelNames_NG <- emModels2
          }
            out <- foreach(i = 1:ncol(NS)) %DO%
            {
              # cindepBIC is the BIC for the grouping model on S and the
              # regression model of the new variable on S
              cindepBIC <- NULL
              try(cindepBIC <- redda_varsel(
                X_p = NS[,i],
                X_c = S,
                cltrain = cltrain,
                alpha_Xtrain = alpha_Xtrain,
                modelNames = modelNames_NG,
                nsamp = nsamp,
                model = "NG"
              )$Best,
              silent = TRUE)

              cindepBIC <- if(is.null(cindepBIC$bic)) NA else cindepBIC$bic
              # Fit the cluster model on the two variables
              sBIC <- NULL
              try(sBIC <- redda_varsel(
                X_p = NS[,i],
                X_c = S,
                cltrain = cltrain,
                alpha_Xtrain = alpha_Xtrain,
                modelNames = emModels2,
                nsamp = nsamp,
                model = "GR"
              )$Best,
              silent = TRUE)
              # depBIC is the BIC for the grouping model with both variables
              depBIC <- if(is.null(sBIC$bic)) NA else sBIC$bic

              return(list(cindepBIC, depBIC, sBIC$modelName, sBIC$G))
            }
            cindepBIC <- sapply(out, "[[", 1)
            depBIC <- sapply(out, "[[", 2)
            # cdiff is the difference between BIC for the models with
            # variables' being grouping variables versus them being
            # conditionally independent of the grouping
            cdiff <- depBIC - cindepBIC
            # Choose the variable with the largest difference
            m <- max(cdiff[is.finite(cdiff)])
            arg <- which(cdiff==m,arr.ind=TRUE)[1]
            if(cdiff[arg] > BIC.diff)
              { # if this difference is positive add this variable to S
                # and update the grouping model's BICS
                BICS <- depBIC[arg]
                k <- c(colnames(S),colnames(NS)[arg])
                nks <- c(colnames(NS)[-arg])
                info <- rbind(info,
                              c(colnames(NS)[arg], BICS, cdiff[arg],
                                "Add", "Accepted", out[[arg]][3:4]))
                S <- cbind(S,NS[,arg])
                NS <- as.matrix(NS[,-arg])
                colnames(S) <- k
                colnames(NS) <- nks
              } else {
                info <- rbind(info,
                              c(colnames(NS)[arg], depBIC[arg], cdiff[arg],
                                "Add", "Rejected", out[[arg]][3:4]))
              }
          }
      }
    if(verbose) cat("- removing step\n")
    # Removal Step for the special case where S contains only a single variable
    if(ncol(S) == 1){
        cdiff <- 0
        oneBIC <- NA
        if (alpha_Xtrain != 0) {
          best_subset <-
            MASS::cov.mcd(x = S, quantile.used = floor((1 - alpha_Xtrain) * Ntrain))$best
          # I recover the best subset on which the mean and the variance are computed.
          # Afterwards, I perform a robust single component no-grouping normal model using Mclust,
          # this is done because cov.mcd (as well as covMcd in robustbase)
          # use consistency correction that I do not want to include

          try(oneBIC <-
                Mclust(
                  S[best_subset,,drop=F],
                  G = 1,
                  modelNames = emModels1,
                  verbose = FALSE
                )$BIC[1],
              silent = TRUE)
        } else {
          try(oneBIC <- Mclust(
            as.matrix(S),
            G = 1,
            modelNames = emModels1,
            verbose = FALSE
          )$BIC[1],
          silent = TRUE)
        }
        # Difference between maximum BIC for grouping and BIC
        # for no grouping
        cdiff <- c(BICS - oneBIC)
        if(is.na(cdiff)) cdiff <- 0
        # Check if difference is less than BIC.diff
        if(cdiff <= BIC.diff){ # if negative remove the variable from S and set the BIC
            # for the model to NA
            BICS <- NA
            info <- rbind(info, c(colnames(S), BICS, cdiff,
                          "Remove", "Accepted", NA, NA))
            k <- c(colnames(NS),colnames(S))
            NS <- cbind(NS,S)
            S <- NULL
            colnames(NS) <- k
          } else { # Otherwise leave S and BICS alone
            info <- rbind(info, c(colnames(S),
                          info[nrow(info),2], cdiff,
                          "Remove", "Rejected",
                          unlist(info[nrow(info),c(6,7)])))
          }
      } else { # Removal step in general (for all cases except when S is a single
        # variable or empty)
        if(ncol(S) >= 2) {
            # Check if the data is at least 3 dimensional
          name <- if(ncol(S) > 2) emModels2 else emModels1
            out <- foreach(i = 1:ncol(S)) %DO%
            {
              try(sBIC <- redda_varsel(
                X_p = S[,i],
                X_c = S[,-i, drop=FALSE],
                cltrain = cltrain,
                alpha_Xtrain = alpha_Xtrain,
                modelNames = name,
                nsamp = nsamp,
                model = "NG"
              )$Best,
              silent = TRUE)
              # cindepBIC is the BIC for the clustering model on the other
              # variables in S and the regression model of the proposed
              #  variable on the other variables in S
              cindepBIC <- if(is.null(sBIC$bic)) NA else sBIC$bic
              # rdep is the the BIC for the clustering model on the other
              # variables in S (excluding the proposed variable)
              rdep <- if(is.null(sBIC$bic_noreg)) NA else sBIC$bic_noreg
              return(list(cindepBIC, rdep, sBIC$modelName, sBIC$G))
            }
            cindepBIC <- sapply(out, "[[", 1)
            rdep <- sapply(out, "[[", 2)
            # depBIC is the BIC for the grouping model with all
            # variables in S
            depBIC <- BICS
            # cdiff is the difference between BIC for the models with
            # variables' being grouping variables versus them being
            # conditionally independent of the grouping
            cdiff <- depBIC - cindepBIC
            # Choose the variable with the smallest difference
            m <- min(cdiff[is.finite(cdiff)])
            arg <- which(cdiff==m,arr.ind=TRUE)[1]
            if(cdiff[arg] <= BIC.diff)
              { # if this difference is less than BIC.diff remove this
                # variable from S and update the grouping model's BICS
                BICS <- rdep[arg]
                k <- c(colnames(NS),colnames(S)[arg])
                nks <- c(colnames(S)[-arg])
                info <- rbind(info,
                              c(colnames(S)[arg], BICS, cdiff[arg],
                                "Remove", "Accepted", out[[arg]][3:4]))
                NS <- cbind(NS,S[,arg])
                S <- as.matrix(S[,-arg])
                colnames(S) <- nks
                colnames(NS) <- k
              } else { info <- rbind(info,
                           c(colnames(S)[arg], rdep[arg], cdiff[arg],
                             "Remove", "Rejected", out[[arg]][3:4]))
              }
          }
      }
    info$BIC <- as.numeric(info$BIC)
    info$BICdiff <- as.numeric(info$BICdiff)

    if(verbose)
      print(info[seq(nrow(info)-1,nrow(info)),c(1,3:5),drop=FALSE])
    # Check if the variables in S have changed or not
    check2 <- colnames(S)
    if(is.null(check2)) # all variables have been removed
      { criterion <- 0 }
    else
      # if they have changed (either added one or removed one or changed one)
      # then continue the algorithm (criterion is 1) otherwise stop
      # (criterion is 0)
      { if(length(check2) != length(check1))
          { criterion <- 1 }
        else
          { criterion <- if(sum(check1==check2) != length(check1)) 1 else 0 }
      }
  }

  if(iter >= itermax)
    warning("Algorithm stopped because maximum number of iterations was reached")

  # List the selected variables and the matrix of steps' information
  info$BIC <- as.numeric(info$BIC)
  info$BICdiff <- as.numeric(info$BICdiff)
  # reorder steps.info
  info <- info[,c(1,4,2,6,7,3,5),drop=FALSE]
  colnames(info) <- c("Variable proposed", "Type of step",
                      "BICclust", "Model", "G", "BICdiff", "Decision")
  varnames <- colnames(Xtrain)
  subset <- if(is.null(S)) NULL
            else sapply(colnames(S), function(x) which(x == varnames))

  out <- list(variables = varnames,
              subset = subset,
              steps.info = info,
              search = "greedy",
              direction = "forward")

  class(out) <- "clustvarsel"
  return(out)
}

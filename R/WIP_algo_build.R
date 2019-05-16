#
# # Testing Raedda on real datasets -----------------------------------------
#
# library(caret)
# library(mclust)
# library(pgmm)
# library(dplyr)
# library(raedda) #package_under_development
# library(BMA)
# library(reddavarsel)
#
# # Data Loading ------------------------------------------------------------
# #
# data("iris")
# df <- iris[, -5]
# class <- iris[, 5]
#
# # data(mclust::banknote)
# # df <- banknote[, -1]
# # class <- banknote[, 1]
#
# # data("wine")
# # df <- (wine[,-1])
# # class <- wine$Type
# #
# # data("thyroid")
# # df <- thyroid[,-1]
# # class <- thyroid[,1]
#
#
# # Case 0: Tr Te split and model fitting original dataset -------------------------------------------
# perc <- .7
# inTrain <- createDataPartition(y = class, times = 1, p = perc)$Resample1
#
# X <- df[inTrain,]
# Y <- df[-inTrain,]
# class_x <- class[inTrain]
# class_y <- class[-inTrain]
#
# X_p <- XX[,4]
# X_c <- XX[,-4]
# cltrain = class_x
# alpha_Xtrain = 0.05
# modelName = "EEE"
# nsamp = 50
# verbose = F
# model <- "NG"
#
# #check whether redda_varsel_works
#
# #0 trimming
# edda_var_sel <- redda_varsel(
#   X_p = X_p,
#   X_c = X_c,
#   cltrain = class_x,
#   alpha_Xtrain = 0,
#   modelNames = NULL,
#   nsamp = 1,
#   model = "GR"
# )
#
# edda <- RAEDDA_l(
#   Xtrain =  X,
#   cltrain = class_x,
#   alpha_Xtrain = 0,
#   modelNames = NULL,
#   nsamp = 1
# )
# edda_var_sel$Best$bic
# edda_var_sel$Best$bic_noreg
# edda$Best$bic
#
# #Can I do it with just one variable?
# edda_var_sel <- redda_varsel(
#   X_p = NULL,
#   X_c = XX[,1],
#   cltrain = class_x,
#   alpha_Xtrain = 0,
#   modelNames = NULL,
#   nsamp = 1,
#   model = "GR"
# )
#
# edda <- RAEDDA_l(
#   Xtrain =  XX[,1],
#   cltrain = class_x,
#   alpha_Xtrain = 0,
#   modelNames = NULL,
#   nsamp = 1
# )
# edda_var_sel$Best$bic
# edda$Best$bic
#
# #No trimming and No grouping
#
# edda_var_sel <- redda_varsel(
#   X_p = X_p,
#   X_c = X_c,
#   cltrain = class_x,
#   alpha_Xtrain = 0,
#   modelNames = NULL,
#   nsamp = 1,
#   model = "NG"
# )
#
# edda <- RAEDDA_l(
#   Xtrain =  X_c,
#   cltrain = class_x,
#   alpha_Xtrain = 0,
#   modelNames = NULL,
#   nsamp = 1
# )
#
# edda_var_sel$Best$bic-edda$Best$bic
# #BIC regression
# reddavarsel::BICreg(reddavarsel::fit_reg(x=X_c, y = X_p))
# edda_var_sel$Best$bic_noreg
# edda$Best$bic
#
# #With trimming
# set.seed(33)
# edda_var_sel <- redda_varsel(
#   X_p = X_p,
#   X_c = X_c,
#   cltrain = class_x,
#   alpha_Xtrain = 0.05,
#   modelNames = NULL,
#   nsamp = 10,
#   model = "GR"
# )
#
# set.seed(33)
# edda <- RAEDDA_l(
#   Xtrain =  X,
#   cltrain = class_x,
#   alpha_Xtrain = 0.05,
#   modelNames = NULL,
#   nsamp = 100
# )
# edda_var_sel$Best$bic
# edda$Best$bic
#
# # model==NG
#
# edda_var_sel <- redda_varsel(
#   X_p = X_p,
#   X_c = X_c,
#   cltrain = class_x,
#   alpha_Xtrain = 0.05,
#   modelName = NULL,
#   nsamp = 100,
#   model = "NG"
# )
#
# edda_var_sel$Best$parameters$variance
#
#
# # First Step- Selecting Single Variable (robust) --------------------------
#
#
# i <- 2
# oneX <- XX[,2]
# oneModel <- Mclust(oneX, G = 1)
#
# oneModel$parameters$mean
# mean(oneX)
#
# oneModel$parameters$variance$sigmasq
# var(oneX)*(length(oneX)-1)/length(oneX)
#
# oneModel$loglik
# (ll_oneModel <- sum(dnorm(
#   oneX,
#   mean = mean(oneX),
#   sd = sqrt(oneModel$parameters$variance$sigmasq),
#   log = T
# ))
# )
#
# oneModel$BIC
# (oneBIC <- 2*ll_oneModel-2*log(length(oneX)))
#
# #Robust_Version
# alpha_Xtrain <- 0.1
# Ntrain <- length(XX[,i])
# Ntrain_trim <- if(alpha_Xtrain!=0) {Ntrain - ceiling(Ntrain * alpha_Xtrain)} else {NULL}
#
# if(alpha_Xtrain!=0){
#   best_subset <-
#     MASS::cov.mcd(x = XX[, i], quantile.used = floor((1 - alpha_Xtrain) * Ntrain))$best
#   # I recover the best subset on which the mean and the variance are computed.
#   # Afterwards, I perform a robust single component no-cluster normal model using Mclust,
#   # this is done because cov.mcd (as well as covMcd in robustbase)
#   # use consistency correction that I do not want to include
#
#   try(oneBIC <- Mclust(XX[best_subset,i], G = 1, modelNames = emModels1,
#                        initialization = list(subset = sub),
#                        verbose = FALSE)$BIC[1],
#       silent = TRUE)
# } else {
#   try(oneBIC <- Mclust(XX[,i], G = 1, modelNames = emModels1,
#                        initialization = list(subset = sub),
#                        verbose = FALSE)$BIC[1],
#       silent = TRUE)
# }
#
#
# # Var_sel_forward_building ------------------------------------------------
#
# library(MASS)
# library(mclust)
# library(reddavarsel)
# # library(parallel)
# # library(foreach)
#
# n <- 200
# pro <- 0.5
# mu1 <- c(0,0)
# mu2 <- c(3,3)
# sigma1 <- matrix(c(1,0.5,0.5,1),2,2,byrow=TRUE)
# sigma2 <- matrix(c(1.5,-0.7,-0.7,1.5),2,2,byrow=TRUE)
# XX <- matrix(0, n, 5)
# colnames(XX) <- paste("XX", 1:ncol(XX), sep ="")
# # generate the grouping variables
# u <- runif(n)
# Class <- ifelse(u < pro, 1, 2)
# XX[u < pro, 1:2]  <- mvrnorm(sum(u < pro), mu = mu1, Sigma = sigma1)
# XX[u >= pro, 1:2] <- mvrnorm(sum(u >= pro), mu = mu2, Sigma = sigma2)
# # generate the non-grouping variables
# XX[,3] <- XX[,1] + rnorm(n)
# XX[,4] <- rnorm(n, mean = 1.5, sd = 2)
# XX[,5] <- rnorm(n, mean = 2, sd = 1)
#
# # plot the data
# clPairs(XX, Class, gap = 0)
#
#
# #Try to see whether it works alpha_Xtrain=0
# group_var_sel_no_trim <-
#   redda_varsel_gr_fwd(Xtrain = XX,
#                       cltrain = Class,
#                       alpha_Xtrain = 0)
#
# #Try to see whether it works alpha_Xtrain=0.1
# t_var_sel <- system.time(
# group_var_sel_trim <-
#   redda_varsel_gr_fwd(Xtrain = XX,
#                       cltrain = Class,
#                       alpha_Xtrain = 0.1,
#                       emModels2 = "EVE")
# )[3]
# #Function arguments
# Xtrain <- X
# cltrain <- Class
# emModels1 = c("E", "V")
# emModels2 = mclust.options("emModelNames")
# alpha_Xtrain <- .05
# forcetwo = TRUE
# BIC.diff = 0
# itermax = 100
# nsamp = 100
# parallel = FALSE
# fit = TRUE
# verbose = TRUE

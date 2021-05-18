# # # Var_sel_forward_building ------------------------------------------------
# # 
# library(MASS)
# library(mclust)
# library(varselTBIC)
# library(parallel)
# library(foreach)
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
#   redda_varsel_greedy_forward(Xtrain = XX,
#                       cltrain = Class,
#                       alpha_Xtrain = 0)
# 
# #Try to see whether it works alpha_Xtrain=0.1
# t_var_sel <- system.time(
#   group_var_sel_trim <-
#     redda_varsel_greedy_forward(
#       Xtrain = XX,
#       cltrain = Class,
#       alpha_Xtrain = 0.1,
#       emModels2 = NULL,
#       nsamp = 5
#     )
# )[3]
# 
# 
# # group_var_sel_no_trim <-
# #   redda_varsel_gr_fwd_union_trimming(Xtrain = XX,
# #                       cltrain = Class,
# #                       alpha_Xtrain = 0)
# #
# # #Try to see whether it works alpha_Xtrain=0.1
# # t_var_sel <- system.time(
# #   group_var_sel_trim <-
# #     redda_varsel_gr_fwd_union_trimming(
# #       Xtrain = XX,
# #       cltrain = Class,
# #       alpha_Xtrain = 0.05,
# #       emModels2 = NULL,
# #       nsamp = 2
# #     )
# # )[3]
# #
# #
# # #Try to see whether it works alpha_Xtrain=0.1 and trimming only for grouping model
# #
# # t_var_sel <- system.time(
# #   group_var_sel_trim <-
# #     redda_varsel_gr_fwd_trimming_full_model(
# #       Xtrain = XX,
# #       cltrain = Class,
# #       alpha_Xtrain = 0.05,
# #       emModels2 = NULL,
# #       nsamp = 2
# #     )
# # )[3]
# 
# #Function arguments
# 
# Xtrain <- XX
# cltrain <- Class
# emModels1 = c("E", "V")
# emModels2 = mclust.options("emModelNames")
# alpha_Xtrain <- .005
# forcetwo = TRUE
# BIC.diff = 0
# itermax = 100
# nsamp = 50
# parallel = FALSE
# fit = TRUE
# verbose = TRUE
# # 

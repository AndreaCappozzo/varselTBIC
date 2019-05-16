# Computes BIC for the regression of y on x, after doing subset selection
# on the predictors
# require(BMA)

BICreg <- function(mod)
{
  sigma <- sqrt(mean(residuals(mod)^2))
  n <- length(mod$residuals)
  p <- n - df.residual(mod) + 1
  # calculate the BIC for the regression
  -n*log(2*pi) -2*n*log(sigma) -n -log(n)*p
}

fit_reg <- function(x, y)
{
  x <- as.matrix(x)
  y <- as.vector(y)
  mod <- bicreg(y = y, x = x, nbest = 1)
  subset <- which(mod$which[1,])
  lm.fit(y = y, x = cbind(1,x[,subset, drop=FALSE]))
}

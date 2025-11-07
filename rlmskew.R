source("Basic Functions.R")

#---------------------------------------------------------------------------*;
# Efficiency factor
#---------------------------------------------------------------------------*;
# Input:
  # x: error
  # loss: type of loss
  # ...: value of tuning parameter, 
    # one for Huber and Tukey losses, two for Asymmetric Huber and Tukey losses

# Output:
  # tau.hat: estimate of efficiency factor

tau <- function(x, loss, ...){
  if (loss == "H") {
    args <- list(...)
    if (length(args) != 1) {
      stop("When loss = 'H', exactly one additional parameter (k) is required")
    }
    k <- args[[1]]
    A = mean(psi2.AH(x, k, k))
    B = mean(psi.d.AH(x, k, k))
    return(tau.hat = A / B^2)
  } else if (loss == "AH") {
    args <- list(...)
    if (length(args) != 2) {
      stop("When loss = 'AH', exactly two additional parameters (k1, k2) are required")
    }
    k1 <- args[[1]]
    k2 <- args[[2]]
    A = mean(psi2.AH(x, k1, k2))
    B = mean(psi.d.AH(x, k1, k2))
    return(tau.hat = A / B^2)
  } else if (loss == "T"){
    args <- list(...)
    if (length(args) != 1) {
      stop("When loss = 'H', exactly one additional parameter (k) is required")
    }
    k <- args[[1]]
    A = mean(psi2.AT(x, k, k))
    B = mean(psi.d.AT(x, k, k))
    return(tau.hat = A / B^2)
  } else if (loss == "AT"){
    args <- list(...)
    if (length(args) != 2) {
      stop("When loss = 'AH', exactly two additional parameters (k1, k2) are required")
    }
    k1 <- args[[1]]
    k2 <- args[[2]]
    A = mean(psi2.AT(x, k1, k2))
    B = mean(psi.d.AT(x, k1, k2))
    return(tau.hat = A / B^2)
  } else {
    stop("Invalid loss type. Must be 'H', 'AH', 'T' or 'AT'")
  }
}

#---------------------------------------------------------------------------*;
# Select tuning parameters
#---------------------------------------------------------------------------*;
# Input:
#  err: error
#  loss: type of loss

# Output:
#  k.sel: selected tuning parameter for the Huber or Tukey losses
#  k1.sel: selected tuning parameter corresponding to negative residuals in AH or AT losses
#  k2.sel: selected tuning parameter corresponding to positive residuals in AH or AT losses
#  tau.hat: estimate of efficiency factor

sel.tun <- function(err, loss = c("H", "AH", "T", "AT")){
  sig.mad <- median(abs(err - median(err))) * 1.4826
  err <- err / sig.mad
  
  if(loss %in% c("H", "AH")){
    ini = 0.1
    end = 3
    if(loss == "H"){
      k <- seq(ini, end, 0.1)
      tau.value <- sapply(k, function(k) tau(x = err, loss = "H", k))
      tau.sel <- tau.value[which.min(tau.value)]
      k.sel <- k[which.min(tau.value)]
      return(list(k.sel = k.sel, tau.hat = tau.sel))
    }else{
      k1 <- seq(ini, end, by = 0.1)
      k2 <- k1
      grid <- expand.grid(k1 = k1, k2 = k2)
      tau.value <- mapply(function(k1, k2) tau(x = err, loss = "AH", k1, k2), grid$k1, grid$k2)
      tau.sel <- tau.value[which.min(tau.value)]
      k1.sel <- grid$k1[which.min(tau.value)]
      k2.sel <- grid$k2[which.min(tau.value)]
      return(list(tau.hat = tau.sel, k1.sel = k1.sel, k2.sel = k2.sel))
    }
  }else if(loss %in% c("T", "AT") ){
    ini = 0.1
    end = 6
    if(loss == "T"){
      k <- seq(ini, end, 0.1)
      tau.value <- sapply(k, function(k) tau(x = err, loss = "T", k))
      tau.sel <- tau.value[which.min(tau.value)]
      k.sel <- k[which.min(tau.value)]
      return(list(k.sel = k.sel, tau.hat = tau.sel))
    }else{
      k1 <- seq(ini, end, by = 0.1)
      k2 <- k1
      grid <- expand.grid(k1 = k1, k2 = k2)
      tau.value <- mapply(function(k1, k2) tau(x = err, loss = "AT", k1, k2), grid$k1, grid$k2)
      tau.sel <- tau.value[which.min(tau.value)]
      k1.sel <- grid$k1[which.min(tau.value)]
      k2.sel <- grid$k2[which.min(tau.value)]
      return(list(tau.hat = tau.sel, k1.sel = k1.sel, k2.sel = k2.sel))
    }
  }
}

#---------------------------------------------------------------------------*;
# Weight function for iteratively reweighted least squares (IRLS)
#---------------------------------------------------------------------------*;
# Input:
#  x: error
#  loss: type of loss

# Output:
#  w.AH(x, k, k): weight function when the loss is the Huber loss
#    x: error
#    k: tuning parameter
#  w.AH(x, k1, k2): weight function when the loss is the asymmetric Huber loss
#    x: error
#    k1: tuning parameter for the negative errors
#    k2: tuning parameter for the positive errors
#  w.AT(x, k, k): weight function when the loss is the Tukey loss
#    x: error
#    k: tuning parameter
#  w.AT(x, k1, k2): weight function when the loss is the asymmetric Tukey loss
#    x: error
#    k1: tuning parameter for the negative errors
#    k2: tuning parameter for the positive errors

weight <- function(x, loss, ...) {
  if (loss == "H") {
    args <- list(...)
    if (length(args) != 1) {
      stop("When loss = 'H', exactly one additional parameter (k) is required")
    }
    k <- args[[1]]
    w.AH(x, k, k)
  } else if (loss == "AH") {
    args <- list(...)
    if (length(args) != 2) {
      stop("When loss = 'AH', exactly two additional parameters (k1, k2) are required")
    }
    k1 <- args[[1]]
    k2 <- args[[2]]
    w.AH(x, k1, k2)
  } else if (loss == "T"){
    args <- list(...)
    if (length(args) != 1) {
      stop("When loss = 'T', exactly one additional parameter (k) is required")
    }
    k <- args[[1]]
    w.AT(x, k, k)
  } else if (loss == "AT"){
    args <- list(...)
    if (length(args) != 2) {
      stop("When loss = 'AT', exactly two additional parameters (k1, k2) are required")
    }
    k1 <- args[[1]]
    k2 <- args[[2]]
    w.AT(x, k1, k2)
  } else {
    stop("Invalid loss type. Must be 'H', 'AH', 'T' or 'AT'")
  }
}

#---------------------------------------------------------------------------*;
# Robust linear model for skew data
#---------------------------------------------------------------------------*;
# Input:
#  formula: regression formula
#  data: data.frame containing response and covariance variables
#  loss: type of loss
#  tuning: selection method of tuning parameters
#  scale.est: estimation method of the scale
#  itmax: maximum number of iterations
#  acc: accuracy of the coefficients stopping iteration
#  ...: k, k1, k2 are given? 

# Output:
#  coeff: estimates of the coefficients
#  coefficients: estimates of the coefficients and its standard error
#  k.obj: used tuning parameters
#  tau.hat: estimated efficiency factor

rlmskew <- function(formula, data, loss = c("H", "T", "AH", "AT"), tuning = c("fixed", "asy"),
                    scale.est = c("MAD"), itmax = 20, acc = 1e-5, ...){
  y.call <- all.vars(terms(formula, response = TRUE))[1]
  y <- data[[y.call]]
  X <- model.matrix(formula, data)
  
  require(MASS)
  lad <- rlm(formula, data = data)
  res0 <- resid(lad)
  coeff0 <- coef(lad)
  
  all.args <- as.list(match.call())[-1]
  tuning.test <- "tuning" %in% names(all.args)
  user.args <- list(...)
  k.test <- "k" %in% names(user.args)
  k1.test <- "k1" %in% names(user.args)
  k2.test <- "k2" %in% names(user.args)
  
  if(loss == "H"){
    if(k.test){
      k = user.args$k
    }else if(!tuning.test){
      stop("Please provide 'k' for Huber loss or specify either 'fixed' or 'asy' for the 'tuning' argument.")
    }else if(tuning == "fixed"){
      k <- 1.345
      if(any(k1.test, k2.test)) warning("Tuning parameters are defined but unused. Please provide 'k' for Huber loss.")
    }else if(tuning == "asy"){
      tun.res <- sel.tun(res0, loss = "H")
      k <- tun.res$k.sel
      if(any(k1.test, k2.test)) warning("Tuning parameters are defined but unused. Please provide 'k' for Huber loss.")
    }
    w <- function(x) weight(x, loss = "H", k = k)
    tau.est <- function(x) tau(x, loss = "H", k = k)
    k.obj <- c(k = k)
  }else if(loss == "AH"){
    if(k1.test & k2.test){
      k1 = user.args$k1
      k2 = user.args$k2
    }else if(!tuning.test){
      stop("Please provide 'k1' and 'k2' for asymmetric Huber loss or specify either 'fixed' or 'asy' for the 'tuning' argument.")
    }else if(tuning == "fixed"){
      k1 <- 1.345
      k2 <- k1
      if(any(k1.test, k2.test, k.test)) warning("Tuning parameters are defined but unused. Please provide 'k1' and 'k2' for asymmetric Huber loss.")
    }else if(tuning == "asy"){
      tun.res <- sel.tun(res0, loss = "AH")
      k1 <- tun.res$k1.sel
      k2 <- tun.res$k2.sel
      if(any(k1.test, k2.test, k.test)) warning("Tuning parameters are defined but unused. Please provide 'k1' and 'k2' for asymmetric Huber loss.")
    }
    w <- function(x) weight(x, loss = "AH", k1 = k1, k2 = k2)
    tau.est <- function(x) tau(x, loss = "AH", k1 = k1, k2 = k2)
    k.obj <- c(k1 = k1, k2 = k2)
  }else if(loss == "T"){
    if(k.test){
      k = user.args$k
    }else if(!tuning.test){
      stop("Please provide 'k' for Tukey loss or specify either 'fixed' or 'asy' for the 'tuning' argument.")
    }else if(tuning == "fixed"){
      k <- 4.685
      if(any(k1.test, k2.test)) warning("Tuning parameters are defined but unused. Please provide 'k' for Tukey loss.")
    }else if(tuning == "asy"){
      k <- sel.tun(res0, loss = "T")$k.sel
      if(any(k1.test, k2.test)) warning("Tuning parameters are defined but unused. Please provide 'k' for Tukey loss.")
    }
    w <- function(x) weight(x, loss = "T", k = k)
    tau.est <- function(x) tau(x, loss = "T", k = k)
    k.obj <- c(k = k)
  }else if(loss == "AT"){
    if(k1.test & k2.test){
      k1 = user.args$k1
      k2 = user.args$k2
    }else if(!tuning.test){
      stop("Please provide 'k1' and 'k2' for asymmetric Tukey loss or specify either 'fixed' or 'asy' for the 'tuning' argument.")
    }else if(tuning == "fixed"){
      k1 <- 4.685
      k2 <- k1
      if(any(k1.test, k2.test, k.test)) warning("Tuning parameters are defined but unused. Please provide 'k1' and 'k2' for asymmetric Tukey loss.")
    }else if(tuning == "asy"){
      tun.res <- sel.tun(res0, loss = "AT")
      k1 <- tun.res$k1.sel
      k2 <- tun.res$k2.sel
      if(any(k1.test, k2.test, k.test)) warning("Tuning parameters are defined but unused. Please provide 'k1' and 'k2' for asymmetric Tukey loss.")
    }
    w <- function(x) weight(x, loss = "AT", k1 = k1, k2 = k2)
    tau.est <- function(x) tau(x, loss = "AT", k1 = k1, k2 = k2)
    k.obj <- c(k1 = k1, k2 = k2)
  }
  
  coeff <- coeff0
  for (s in 1:itmax) {
    oldcoeff <- coeff
    res <- y - X %*% coeff
    sig.mad <- median(abs(res - median(res))) * 1.4826
    W <- w(res / sig.mad)
    coeff <- solve(t(X) %*% diag(W) %*% X) %*% t(X) %*% diag(W) %*% y
    if(any(is.na(coeff))){
      coeff <- oldcoeff
      break
    }else if(sum((coeff - oldcoeff)^2) < acc){
      break
    }
  }
  
  fit <- X %*% coeff
  res <- y - X %*% coeff
  sig.mad <- median(abs(res - median(res))) * 1.4826
  tau.hat <- tau.est(res / sig.mad)
  cov.mat <- solve(t(X) %*% diag(W) %*% X / sum(W)) * sig.mad^2 * tau.hat / n
  Std.Error <- sqrt(diag(cov.mat))
  
  Coefficients <- cbind(coeff, Std.Error)
  coeff <- c(coeff)
  names(coeff) <- colnames(X)
  
  return(list(coeff = unlist(coeff), coefficients = Coefficients, k.obj = k.obj, tau.hat = tau.hat))
}

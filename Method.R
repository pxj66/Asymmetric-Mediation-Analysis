setwd("E:/课题/2025 Asymmetric Mediation Analysis/Rcode")
require(quantreg)
source("Basic Functions.R")
source("LocalInfRobustMediation.R")
require(robmed)
require(boot)
source("rlmskew.R")
source("DGPs.R")


boot.f <- function(data, ind, met, ...){
  if(met == "LS"){
    ls1 <- lm(m ~ x, data = data[ind, ])
    ls2 <- lm(y ~ x + m, data = data[ind, ])
    a.hat <- coef(ls1)["x"]
    b.hat <- coef(ls2)["m"]
    a.hat * b.hat
  }else if(met == "LAD"){
    ls1 <- rq(m ~ x, data = data[ind, ])
    ls2 <- rq(y ~ x + m, data = data[ind, ])
    a.hat <- coef(ls1)["x"]
    b.hat <- coef(ls2)["m"]
    a.hat * b.hat
  }else if(met == "HF"){
    rlm1 <- rlmskew(m ~ x, data = data[ind, ], loss = "H", tuning = "fixed")
    rlm2 <- rlmskew(y ~ x + m, data = data[ind, ], loss = "H", tuning = "fixed")
    a.hat <- rlm1$coeff["x"]
    b.hat <- rlm2$coeff["m"]
    a.hat * b.hat
  }else if(met == "HD"){
    args <- list(...)
    k.1 <- args$k.1
    k.2 <- args$k.2
    rlm1 <- rlmskew(m ~ x, data = data[ind, ], loss = "H", tuning = "asy", k = k.1)
    rlm2 <- rlmskew(y ~ x + m, data = data[ind, ], loss = "H", tuning = "asy", k = k.2)
    a.hat <- rlm1$coeff["x"]
    b.hat <- rlm2$coeff["m"]
    a.hat * b.hat
  }else if(met == "AHD"){
    args <- list(...)
    k11 <- args$k11
    k12 <- args$k12
    k21 <- args$k21
    k22 <- args$k22
    rlm1 <- rlmskew(m ~ x, data = data[ind, ], loss = "AH", k1 = k11, k2 = k12, tuning = "asy")
    rlm2 <- rlmskew(y ~ x + m, data = data[ind, ], loss = "AH", k1 = k21, k2 = k22, tuning = "asy")
    a.hat <- rlm1$coeff["x"]
    b.hat <- rlm2$coeff["m"]
    a.hat * b.hat
  }else if(met == "ATD"){
    args <- list(...)
    k11 <- args$k11
    k12 <- args$k12
    k21 <- args$k21
    k22 <- args$k22
    rlm1 <- rlmskew(m ~ x, data = data[ind, ], loss = "AT", k1 = k11, k2 = k12, tuning = "asy")
    rlm2 <- rlmskew(y ~ x + m, data = data[ind, ], loss = "AT", k1 = k11, k2 = k12, tuning = "asy")
    a.hat <- rlm1$coeff["x"]
    b.hat <- rlm2$coeff["m"]
    a.hat * b.hat
  }else if(met == "LI"){
    Z <- as.matrix(data[ind, ])
    muSigRes <- MeanCov(Z)
    MLERes <- MLEst(muSigRes$S)
    MLERes[1] * MLERes[2]
  }else if(met == "RT"){
    Z <- as.matrix(data[ind, ])
    TunRes <- HuberTun(0.05, p = ncol(Z))
    EstRes <- robEst(Z, r = TunRes$r, tau = TunRes$tau, ep = 1e-5)
    EstRes$theta[1] * EstRes$theta[2]
  }
}


Met <- function(met, dgp){
  data = dgp
  n <- nrow(data)
  x <- data$x
  m <- data$m
  y <- data$y
  
  if (met == "LS"){
    ls1 <- lm(m ~ x, data = data)
    ls2 <- lm(y ~ x + m, data = data)
    a.hat <- coef(ls1)["x"]
    b.hat <- coef(ls2)["m"]
    a.se <- summary(ls1)$coefficients[2, 2]
    b.se <- summary(ls2)$coefficients[3, 2]
    IDE <- a.hat * b.hat
    SE.sobel <- sqrt((a.hat * b.se)^2 + (b.hat * a.se)^2)
    test.sobel <- (IDE - qnorm(0.975) * SE.sobel) * (IDE + qnorm(0.975) * SE.sobel) > 0
    ab.boot <- boot(data = data, function(data, ind) boot.f(data, ind, met = "LS"), R = 1000)
    test.perc <- prod(boot.ci(ab.boot, type = "perc")$percent[c(4:5)]) > 0
  }else if (met == "LAD"){
    lad1 <- rq(m ~ x, data = data)
    lad2 <- rq(y ~ x + m, data = data)
    a.hat <- coef(lad1)["x"]
    b.hat <- coef(lad2)["m"]
    a.se <- summary(lad1, se = "iid")$coefficients[2, 2]
    b.se <- summary(lad2, se = "iid")$coefficients[3, 2]
    IDE <- a.hat * b.hat
    SE.sobel <- sqrt((a.hat * b.se)^2 + (b.hat * a.se)^2)
    test.sobel <- (IDE - qnorm(0.975) * SE.sobel) * (IDE + qnorm(0.975) * SE.sobel) > 0
    ab.boot <- boot(data = data, function(data, ind) boot.f(data, ind, met = "LAD"), R = 1000)
    test.perc <- prod(boot.ci(ab.boot, type = "perc")$percent[c(4:5)]) > 0
  }else if (met == "HF"){
    rlm1 <- rlmskew(m ~ x, data = data, loss = "H", tuning = "fixed")
    rlm2 <- rlmskew(y ~ x + m, data = data, loss = "H", tuning = "fixed")
    a.hat <- rlm1$coeff["x"]
    b.hat <- rlm2$coeff["m"]
    a.se <- rlm1$coefficients[2, 2]
    b.se <- rlm2$coefficients[3, 2]
    IDE <- a.hat * b.hat
    SE.sobel <- sqrt((a.hat * b.se)^2 + (b.hat * a.se)^2)
    test.sobel <- (IDE - qnorm(0.975) * SE.sobel) * (IDE + qnorm(0.975) * SE.sobel) > 0
    ab.boot <- boot(data = data, function(data, ind) boot.f(data, ind, met = "HF"), R = 1000)
    test.perc <- prod(boot.ci(ab.boot, type = "perc")$percent[c(4:5)]) > 0
  }else if (met == "HD"){
    rlm1 <- rlmskew(m ~ x, data = data, loss = "H", tuning = "asy")
    rlm2 <- rlmskew(y ~ x + m, data = data, loss = "H", tuning = "asy")
    k.1 = rlm1$k.obj
    k.2 = rlm2$k.obj
    a.hat <- rlm1$coeff["x"]
    b.hat <- rlm2$coeff["m"]
    a.se <- rlm1$coefficients[2, 2]
    b.se <- rlm2$coefficients[3, 2]
    IDE <- a.hat * b.hat
    SE.sobel <- sqrt((a.hat * b.se)^2 + (b.hat * a.se)^2)
    test.sobel <- (IDE - qnorm(0.975) * SE.sobel) * (IDE + qnorm(0.975) * SE.sobel) > 0
    ab.boot <- boot(data = data, function(data, ind) boot.f(data, ind, met = "HD", k.1 = k.1, k.2 = k.2), R = 1000)
    test.perc <- prod(boot.ci(ab.boot, type = "perc")$percent[c(4:5)]) > 0
  }else if (met == "AHD"){
    rlm1 <- rlmskew(m ~ x, data = data, loss = "AH", tuning = "asy")
    rlm2 <- rlmskew(y ~ x + m, data = data, loss = "AH", tuning = "asy")
    a.hat <- rlm1$coeff["x"]
    b.hat <- rlm2$coeff["m"]
    a.se <- rlm1$coefficients[2, 2]
    b.se <- rlm2$coefficients[3, 2]
    IDE <- a.hat * b.hat
    SE.sobel <- sqrt((a.hat * b.se)^2 + (b.hat * a.se)^2)
    test.sobel <- (IDE - qnorm(0.975) * SE.sobel) * (IDE + qnorm(0.975) * SE.sobel) > 0
    ab.boot <- boot(data = data, function(data, ind){
      boot.f(data, ind, met = "AHD", k11 = rlm1$k.obj[1], k12 = rlm1$k.obj[2], 
             k21 = rlm2$k.obj[1], k22 = rlm2$k.obj[2])}, R = 1000)
    test.perc <- prod(boot.ci(ab.boot, type = "perc")$percent[c(4:5)]) > 0
  }else if (met == "ATD"){
    rlm1 <- rlmskew(m ~ x, data = data, loss = "AT", tuning = "asy")
    rlm2 <- rlmskew(y ~ x + m, data = data, loss = "AT", tuning = "asy")
    a.hat <- rlm1$coeff["x"]
    b.hat <- rlm2$coeff["m"]
    a.se <- rlm1$coefficients[2, 2]
    b.se <- rlm2$coefficients[3, 2]
    IDE <- a.hat * b.hat
    SE.sobel <- sqrt((a.hat * b.se)^2 + (b.hat * a.se)^2)
    test.sobel <- (IDE - qnorm(0.975) * SE.sobel) * (IDE + qnorm(0.975) * SE.sobel) > 0
    ab.boot <- boot(data = data, function(data, ind){
      boot.f(data, ind, met = "ATD", k11 = rlm1$k.obj[1], k12 = rlm1$k.obj[2], 
             k21 = rlm2$k.obj[1], k22 = rlm2$k.obj[2])}, R = 1000)
    test.perc <- prod(boot.ci(ab.boot, type = "perc")$percent[c(4:5)]) > 0
  }else if (met == "LI"){
    Z <- as.matrix(data)
    LIRes <- LocalInf(Z)
    data.LI <- data[LIRes$B.v <= 0.8, ]
    
    Z <- as.matrix(data.LI)
    muSigRes <- MeanCov(Z)
    MLERes <- MLEst(muSigRes$S)
    IDE <- MLERes[1] * MLERes[2]
    SE.Res <- SEML(Z, MLERes)
    SE.sobel <- SE.Res$sand[7]
    test.sobel <- (IDE - qnorm(0.975) * SE.sobel) * (IDE + qnorm(0.975) * SE.sobel) > 0
    ab.boot <- boot(data = data.LI, function(data, ind) boot.f(data, ind, met = "LI"), R = 1000)
    test.perc <- prod(boot.ci(ab.boot, type = "perc")$percent[c(4:5)]) > 0
  }else if (met == "RT"){
    Z <- as.matrix(data)
    TunRes <- HuberTun(0.05, p = ncol(Z))
    EstRes <- robEst(Z, r = TunRes$r, tau = TunRes$tau, ep = 1e-5)
    IDE <- EstRes$theta[1] * EstRes$theta[2]
    SE.Res <- SErob(Z = Z, mu = EstRes$mu, Sigma = EstRes$Sigma,
                   theta = EstRes$theta, d = EstRes$d,
                   r = TunRes$r, tau = TunRes$tau)
    data.RT <- data.frame(SE.Res$Zr)
    SE.sobel <- SE.Res$sand[7]
    test.sobel <- (IDE - qnorm(0.975) * SE.sobel) * (IDE + qnorm(0.975) * SE.sobel) > 0
    ab.boot <- boot(data = data.RT, function(data, ind) boot.f(data, ind, met = "RT"), R = 1000)
    test.perc <- prod(boot.ci(ab.boot, type = "perc")$percent[c(4:5)]) > 0
  }else if (met == "FAR"){
    Res <- test_mediation(y ~ x + m(m), data = data, robust = T, type = "perc")
    IDE <- Res$indirect
    test.sobel <- NA
    test.perc <- prod(Res$ci) > 0
  }
  return(list(IDE = IDE, test.sobel = test.sobel, test.perc = test.perc))
}


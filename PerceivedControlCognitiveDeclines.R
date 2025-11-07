# Mediation Analysis Based on the data.frame: data
load("LPACOG.rda")
library(boot)
set.seed(123)

#-------------------------------------------------------------------------------
# LS mediation
#-------------------------------------------------------------------------------
ls1 <- lm(lpa ~ . -cog, data = LPACOG)
ls2 <- lm(cog ~ ., data = LPACOG)
a <- coef(ls1)['pc']
b <- coef(ls2)['lpa']
a * b
se.a <- summary(ls1)$coefficients["pc", 2]
se.b <- summary(ls2)$coefficients["lpa", 2]
sobel.ab <- sqrt(a^2 * se.b^2 + b^2 * se.a^2)
paste0("[", round(a * b - qnorm(0.975) * sobel.ab, digits = 5), ",",
       round(a * b + qnorm(0.975) * sobel.ab, digits = 5), "]")

ls.boot <- function(d, i){
  d0 <- d[i, ]
  ls1 <- lm(lpa ~ . - cog, data = d0)
  ls2 <- lm(cog ~ ., data = d0)
  a <- coef(ls1)['pc']
  b <- coef(ls2)['lpa']
  a * b
}
ls.res <- boot(data = LPACOG, ls.boot, R = 1000)
boot.ci(ls.res, type = "perc")

#-------------------------------------------------------------------------------
# LAD mediation
#-------------------------------------------------------------------------------
library(quantreg)
lad1 <- rq(lpa ~ . -cog, data = LPACOG)
lad2 <- rq(cog ~ ., data = LPACOG)
a <- coef(lad1)['pc']
b <- coef(lad2)['lpa']
a * b
se.a <- summary(lad1, se = "iid")$coefficients["pc", 2]
se.b <- summary(lad2, se = "iid")$coefficients["lpa", 2]
sobel.ab <- sqrt(a^2 * se.b^2 + b^2 * se.a^2)
paste0("[", round(a * b - qnorm(0.975) * sobel.ab, digits = 5), ",",
       round(a * b + qnorm(0.975) * sobel.ab, digits = 5), "]")

lad.boot <- function(d, i){
  d0 <- d[i, ]
  lad1 <- rq(lpa ~ . - cog, data = d0)
  lad2 <- rq(cog ~ ., data = d0)
  a <- coef(lad1)['pc']
  b <- coef(lad2)['lpa']
  # a1 * b
  a * b
}

lad.res <- boot(data = LPACOG, lad.boot, R = 1000)
boot.ci(lad.res, type = "perc")

#-------------------------------------------------------------------------------
# Huber with fixed tuning parameter mediation
#-------------------------------------------------------------------------------
setwd("/Real Data Analysis")
load("LPACOG.rda")
source('Basic Functions.R')
source("rlmskew.R")

hf1 <- rlmskew(lpa ~ .- cog, data = LPACOG, loss = "H", tuning = "fixed")
hf2 <- rlmskew(cog ~ ., loss = "H", tuning = "fixed", data = LPACOG)
a <- hf1$coeff['pc']
b <- hf2$coeff['lpa']
a * b
se.a <- hf1$coefficients["pc", 2]
se.b <- hf2$coefficients["lpa", 2]
sobel.ab <- sqrt(a^2 * se.b^2 + b^2 * se.a^2)
paste0("[", round(a * b - qnorm(0.975) * sobel.ab, digits = 5), ",",
       round(a * b + qnorm(0.975) * sobel.ab, digits = 5), "]")


hf.boot <- function(d, i){
  d0 <- d[i, ]
  hf1 <- rlmskew(lpa ~ .- cog, loss = "H", tuning = "fixed", data = d0)
  hf2 <- rlmskew(cog ~ ., loss = "H", tuning = "fixed", data = d0)
  a <- hf1$coeff['pc']
  b <- hf2$coeff['lpa']
  # a1 * b
  a * b
}
hf.res <- boot(data = LPACOG, hf.boot, R = 1000)
boot.ci(hf.res, type = "perc")

#-------------------------------------------------------------------------------
#  Huber with data-dependent tuning parameter mediation
#-------------------------------------------------------------------------------
setwd("/Real Data Analysis")
load("LPACOG.rda")
library(boot)
source('Basic Functions.R')
source("rlmskew.R")

hd1 <- rlmskew(lpa ~ .- cog, loss = "H", tuning = "asy", data = LPACOG)
hd2 <- rlmskew(cog ~ ., loss = "H", tuning = "asy", data = LPACOG)
a <- hd1$coeff['pc']
b <- hd2$coeff['lpa']
a * b
se.a <- hd1$coefficients["pc", 2]
se.b <- hd2$coefficients["lpa", 2]
sobel.ab <- sqrt(a^2 * se.b^2 + b^2 * se.a^2)
paste0("[", round(a * b - qnorm(0.975) * sobel.ab, digits = 5), ",",
       round(a * b + qnorm(0.975) * sobel.ab, digits = 5), "]")

hd.boot <- function(d, i){
  d0 <- d[i, ]
  hd1 <- rlmskew(lpa ~ .- cog, loss = "H", tuning = "asy", data = d0)
  hd2 <- rlmskew(cog ~ ., loss = "H", tuning = "asy", data = d0)
  a <- hd1$coeff['pc']
  b <- hd2$coeff['lpa']
  a * b
}
hd.res <- boot(data = LPACOG, hd.boot, R = 1000)
boot.ci(hd.res, type = "perc")

#-------------------------------------------------------------------------------
#  Asymmetric Huber with data-dependent tuning parameter mediation
#-------------------------------------------------------------------------------
setwd("GitHub/Real Data Analysis")
load("LPACOG.rda")
source('Basic Functions.R')
source("rlmskew.R")
library(boot)

ahd1 <- rlmskew(lpa ~ .- cog, loss = "AH", tuning = "asy", data = LPACOG)
ahd2 <- rlmskew(cog ~ ., loss = "AH", tuning = "asy", data = LPACOG)
a <- ahd1$coeff['pc']
b <- ahd2$coeff['lpa']
a * b
se.a <- ahd1$coefficients["pc", 2]
se.b <- ahd2$coefficients["lpa", 2]
sobel.ab <- sqrt(a^2 * se.b^2 + b^2 * se.a^2)
paste0("[", round(a * b - qnorm(0.975) * sobel.ab, digits = 5), ",",
       round(a * b + qnorm(0.975) * sobel.ab, digits = 5), "]")

ahd.boot <- function(d, i){
  d0 <- d[i, ]
  ahd1 <- rlmskew(lpa ~ .- cog, loss = "AH", data = d0, k1 = 2.2, k2 = 1)
  ahd2 <- rlmskew(cog ~ ., loss = "AH", data = d0, k1 = 2.1, k2 = 1.2)
  a <- ahd1$coeff['pc']
  b <- ahd2$coeff['lpa']
  a * b
}
ahd.res <- boot(data = LPACOG, ahd.boot, R = 1000)
boot.ci(ahd.res, type = "perc")

#-------------------------------------------------------------------------------
#  Asymmetric Tukey with data-dependent tuning parameter mediation
#-------------------------------------------------------------------------------
setwd("E:/课题/Completed/2025 Asymmetric Mediation Analysis/Rcode/GitHub/Real Data Analysis")
load("LPACOG.rda")
source('Basic Functions.R')
source("rlmskew.R")
library(boot)

atd1 <- rlmskew(lpa ~ .- cog, loss = "AT", tuning = "asy", data = LPACOG)
atd2 <- rlmskew(cog ~ ., loss = "AT", tuning = "asy", data = LPACOG)
a <- atd1$coeff['pc']
b <- atd2$coeff['lpa']
a * b
se.a <- atd1$coefficients["pc", 2]
se.b <- atd2$coefficients["lpa", 2]
sobel.ab <- sqrt(a^2 * se.b^2 + b^2 * se.a^2)
paste0("[", round(a * b - qnorm(0.975) * sobel.ab, digits = 7), ",",
       round(a * b + qnorm(0.975) * sobel.ab, digits = 5), "]")

atd.boot <- function(d, i){
  d0 <- d[i, ]
  atd1 <- rlmskew(lpa ~ .- cog, loss = "AT", k1 = 6, k2 = 5.1, data = d0)
  atd2 <- rlmskew(cog ~ ., loss = "AT", k1 = 6, k2 = 5.2, data = d0)
  a <- atd1$coeff['pc']
  b <- atd2$coeff['lpa']
  a * b
}

atd.res <- boot(data = LPACOG, atd.boot, R = 1000)
boot.ci(ahd.res, type = "perc")

#-------------------------------------------------------------------------------
#  Fast and robust mediation
#-------------------------------------------------------------------------------
library(robmed)
load("LPACOG.rda")
res <- test_mediation(cog ~ m(lpa) + pc + covariates(age, sex, race, educ, inc, pses, limit),
                      data = LPACOG,
                      robust = "MM", type = "perc")
res$indirect

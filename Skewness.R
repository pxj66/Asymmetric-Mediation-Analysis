#setwd("/Real Data Analysis")

load("LPACOG.rda")
source('Basic Functions.R')
source("rlmskew.R")
library(e1071)
library(symmetry)

ahd1 <- rlmskew(lpa ~ .- cog, loss = "AH", tuning = "asy", data = LPACOG)
ahd2 <- rlmskew(cog ~ ., loss = "AH", tuning = "asy", data = LPACOG)

e1 <- ahd1$res
e2 <- ahd2$res

skewness(e1)
kurtosis(e1)
ks.test(e1, "pnorm")
symmetry_test(e1, "MOI", k = 3, mu = 0)
ahd1$k.obj

skewness(e2)
kurtosis(e2)
ks.test(e2, "pnorm")
symmetry_test(e2, "MOI", k = 3, mu = 0)
ahd2$k.obj

skewness(e3)
kurtosis(e3)
ks.test(e3, "pnorm")
symmetry_test(e3, "MOI", k = 3, mu = 0)
ahd3$k.obj
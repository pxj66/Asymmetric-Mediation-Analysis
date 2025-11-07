setwd("E:/课题/2025 Asymmetric Mediation Analysis/Rcode")
source("DGPs.R")
source("Method.R")

Eff = c("s.pow", "m.pow", "l.pow", "s.t1e", "m.t1e", "l.t1e") 
Method = c("LS", "LAD", "HF", "HD", "AHD", "ATD", "LI", "RT", "FAR")
Dist = c("norm", "mix-norm", "laplace", "mix-exp,r=0.1","mix-exp,r=0.3","mix-exp,r=0.5", "LN", "chisq5")
N <- c(50, 200, 1000)
library(data.table)
library(foreach)
library(doParallel)

cl <- makeCluster(12)
registerDoParallel(cl)

Res <- foreach(eff = Eff, .combine = 'rbind') %:%
  foreach(err = Dist, .combine = 'rbind') %:%
  foreach(n = N, .combine = 'rbind') %:%
  foreach(i = 1000, .combine = 'rbind', .packages = c("data.table", "quantreg", "robmed", "boot"))%dopar%{
    Data <- DGP(eff, err, n)
    do.call(rbind, lapply(Method, function(met){
      ResMet <- Met(met, Data)
      data.table(eff = eff, err = err, n = n, i = i, met = met, 
                 IDE = ResMet$IDE, 
                 test.sobel = ResMet$test.sobel, 
                 test.perc = ResMet$test.perc)
    }))
  }


save(Res, file = "Simulation.RData")
stopCluster(cl)

Res[eff == "s.pow" & err == "chisq5", mean((IDE - 0.14 * 0.14)^2), by = .(n, met)]


library(data.table)
#, "mix-t", "laplace", "t2", "mix-exp", "mix-chi", "LN", "chisq"
Results[eff == "s" & err == "chisq", mean((IDE - 0.14 * 0.14)^2), by = .(n, met)]
Results[eff == "m" & err == "chisq", mean((IDE - 0.39 * 0.39)^2), by = .(n, met)]
Results[eff == "l" & err == "chisq", mean((IDE - 0.59 * 0.59)^2), by = .(n, met)]

n = 1000
err = "norm"
eff = "s.pow"
library(foreach)
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)
system.time(
  Res <- foreach(i = 1:1000, .combine = 'rbind', .packages = c("data.table", "quantreg", "robmed", "boot"))%dopar%{
  Data <- DGP(eff, err, n)
  do.call(rbind, lapply(Method, function(met){
    ResMet <- Met(met, Data)
    data.table(eff = eff, err = err, n = n, i = i, met = met, 
               IDE = ResMet$IDE, 
               test.sobel = ResMet$test.sobel, 
               test.perc = ResMet$test.perc)
  }))
})

































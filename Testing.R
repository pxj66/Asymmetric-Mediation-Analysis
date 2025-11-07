library(foreach)
library(doSNOW)

source("DGPs.R")
source("Method.R")

Eff = c("s.pow", "m.pow", "l.pow", "s.t1e", "m.t1e", "l.t1e") 
Method = c("LS", "LAD", "HF", "HD", "AHD", "ATD", "LI", "RT", "FAR")
Dist = c("norm", "mix-norm", "laplace", "mix-exp,r=0.1","mix-exp,r=0.3","mix-exp,r=0.5", "LN", "chisq5")
N <- c(50, 200, 1000)

for (eff in Eff) {
  for (err in Dist) {
    time.start <- Sys.time()
    for (n in 1000) {
      time0 <- Sys.time()
      #-------------------------------------------------------------------------
      cl <- makeSOCKcluster(12)
      registerDoSNOW(cl)
      
      Data <- DGP(eff, err, n)
      
      pb <- txtProgressBar(max = 1000, style = 3)
      progress = function(n) setTxtProgressBar(pb, n)
      
      Res <- foreach(i = 1:1000, 
                     .combine = "rbind", 
                     .options.snow = list(progress = progress), 
                     .packages = c("quantreg", "data.table", "robmed", "boot")
                     )%dopar%{
                       R <- c()
                       for (met in Method) {
                         R0 <- Met(met, Data)
                         R1 <- data.table(eff = eff, err = err, n = n, i = i, met = met, 
                                          IDE = R0$IDE, 
                                          test.sobel = R0$test.sobel, 
                                          test.perc = R0$test.perc)
                         R <- rbind(R, R1)
                       }
                       R
                     }
      close(pb)
      stopCluster(cl)
      file.name <- sprintf("eff=%s,err=%s,n=%d.RData", eff, err, n)
      save(Res, file = file.name)
      #-------------------------------------------------------------------------------
      time1 <- Sys.time()
      print(time1 - time0)
      cat("Eff =", eff, "Err =", err, "N =", n, "\n", sep = " ")
      print(time1 - time.start)
    }
  }
}


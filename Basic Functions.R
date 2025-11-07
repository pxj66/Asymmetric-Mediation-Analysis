## Basic functions for asymmetric Huber loss
psi.AH <- function(x, k1, k2){
  sapply(x, function(t){
    y <- numeric(0)
    if (t < -k1){
      y <- - k1
    }else if (t < k2){
      y <- t
    }else{
      y <- k2
    }
    return(y)
  })
}

w.AH <- function(x, k1, k2){
  sapply(x, function(t){
    y <- numeric(0)
    if (t < -k1){
      y <- - k1 / t
    }else if (t < k2){
      y <- 1
    }else{
      y <- k2 / t
    }
    return(y)
  })
}

psi2.AH <- function(x, k1, k2){
  sapply(x, function(t){
    y <- numeric(0)
    if (t < -k1){
      y <- k1^2
    }else if (t < k2){
      y <- t^2
    }else{
      y <- k2^2
    }
    return(y)
  })
}

psi.d.AH <- function(x, k1, k2){
  sapply(x, function(t){
    y <- numeric(0)
    if (t < -k1){
      y <- 0
    }else if (t < k2){
      y <- 1
    }else{
      y <- 0
    }
    return(y)
  })
}

## Basic functions for asymmetric Tukey loss
psi.AT <- function(x, k1, k2){
  sapply(x, function(t){
    y <- numeric(0)
    if (t < -k1){
      y <- 0
    }else if (t < 0){
      y <- t - 2 * t^3 / k1^2 + t^5 / k1^4
    }else if (t < k2){
      y <- t - 2 * t^3 / k2^2 + t^5 / k2^4
    }else{
      y <- 0
    }
    return(y)
  })
}

w.AT <- function(x, k1, k2){
  sapply(x, function(t){
    y <- numeric(0)
    if (t < -k1){
      y <- 0
    }else if (t < 0){
      y <- 1 - 2 * t^2 / k1^2 + t^4 / k1^4
    }else if (t < k2){
      y <- 1 - 2 * t^2 / k2^2 + t^4 / k2^4
    }else{
      y <- 0
    }
    return(y)
  })
}

psi2.AT <- function(x, k1, k2){
  sapply(x, function(t){
    y <- numeric(0)
    if (t < -k1){
      y <- 0
    }else if (t < 0){
      y <- (t - 2 * t^3 / k1^2 + t^5 / k1^4)^2
    }else if (t < k2){
      y <- (t - 2 * t^3 / k2^2 + t^5 / k2^4)^2
    }else{
      y <- 0
    }
    return(y)
  })
}

psi.d.AT <- function(x, k1, k2){
  sapply(x, function(t){
    y <- numeric(0)
    if (t < -k1){
      y <- 0
    }else if (t < 0){
      y <- 1 - 6 * t^2 / k1^2 + 5 * t^4 / k1^4
    }else if (t < k2){
      y <- 1 - 6 * t^2 / k2^2 + 5 * t^4 / k2^4
    }else{
      y <- 0
    }
    return(y)
  })
}

psi.H <- function(x, k){
  x[x > k]  <- k
  x[x < -k] <- -k
  return(x = c(x))
}

w.H <- function(x, k){
  ind1 <- x > k
  ind2 <- x < -k
  ind <- abs(x) <= k
  x[ind] <- 1
  x[ind1] <- k / x[ind1]
  x[ind2] <- -k / x[ind2]
  return(x = c(x))
}

psi.dot.H <- function(x, k){
  ind <- abs(x) <= k
  x[ind] <- 1
  x[!ind] <- 0
  return(x = c(x))
}

psi2.H <- function(x, k){
  ind <- abs(x) > k
  x[ind] <- k^2
  x[!ind] <- x[!ind]^2
  return(x = c(x))
}

plaplace <- function(q, mu = 0, b = 1){
  if(q <= mu){
    0.5 * exp((q - mu) / b)
  }else{
    1 - 0.5 * exp(- (q - mu) / b)
  }
}

inv.laplace <- function(p, mu = 0, b = 1){
  if(p <= 0.5){
    mu + b * log(2 * p)
  }else{
    mu - b * log(2 * (1 - p))
  }
}

rlaplace <- function(n, mu = 0, b = 1){
  u <- runif(n)
  sapply(u, function(p) inv.laplace(p, mu = mu, b = b))
}













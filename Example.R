source("rlmskew.R")

n <- 100
x <- rnorm(n)
m <- 0.14 * x + rnorm(n)
y <- x + 0.14 * m + rnorm(n)
data <- data.frame(x = x, m = m, y = y)

rlmskew(formula = y ~ m + x, data = data, loss = "H", tuning = "fixed")
rlmskew(formula = y ~ m + x, data = data, loss = "AH", tuning = "fixed")

rlmskew(formula = y ~ m + x, data = data, loss = "H", tuning = "asy")
rlmskew(formula = y ~ m + x, data = data, loss = "AH", tuning = "asy")


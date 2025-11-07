
X <- matrix(c(rep(1, 6), -4, -3, -2, -1, 0, 10), nrow = 6, byrow = F)
y <- c(2.48, 0.73, -0.04, -1.44, -1.32, 0)
H <- X %*% solve(t(X) %*% X) %*% t(X)

y - H %*% y
diag(H)[diag(H) > sum(diag(H)) / 6]


# mediation
library(VGAM)
n = 100
x <- rnorm(n)
m <- 0.14 * x + rlaplace(n)
y <- x + 0.14 * m + rlaplace(n)

Z <- matrix(c(rep(1, n), x, m), nrow = n, byrow = F)
Z_add <- rbind(Z, c(1, x[n], m[n]))

H <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
H_tilda <- Z_add %*% solve(t(Z_add) %*% Z_add) %*% t(Z_add)


1 / tail(diag(H), 1) - 1 / tail(diag(H_tilda), 1)

tail(diag(H), 1) / (1 + tail(diag(H), 1))
tail(diag(H_tilda), 1)




rank(diag(H) * (1 - diag(H)))
rank(diag(H))
boxplot(y - H %*% y)


x <- seq(0, 1.5, 0.01)
y <- sqrt(x)
plot(x, y, type = "l", xlab = "h")
abline(h = 1, col = "gray", lty = "dashed")
points(x = .2, y = sqrt(.2), col = "blue")













# irdec examples ----
## continuous case ----
n = 1000
x <- matrix(runif(n * 3), nrow = n)
y <- (x[, 1] + x[, 2]) %% 1
irdec(y, x[, 1])
irdec(y, x[, 2])
irdec(y, x[, 3])

## discrete case
n <- 10000
s <- 0.1
x1 <- c(rep(0,n*s), runif(n*(1-s)))
x2 <- runif(n)
y <- x1

irdec(y,x1,dist.type.X = "discrete")
irdec(y,x2)

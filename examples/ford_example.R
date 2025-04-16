# ford example ----
n = 1000
p = 100
x <- matrix(rnorm(n * p), nrow = n)
colnames(x) = paste0(rep("x", p), seq(1, p))
y <- x[, 1] * x[, 10] + x[, 20]^2
# with num_features equal to 3 and stop equal to FALSE, ford will give a list of
# three selected features
result1 = ford(y, x, num_features = 3, stop = FALSE, numCores = 1)
result1



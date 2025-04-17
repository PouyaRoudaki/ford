# irdec examples ----
## continuous case ----
n = 1000
x <- matrix(runif(n * 3), nrow = n)
y <- (x[, 1] + x[, 2]) %% 1
irdec(y, x[, 1])
irdec(y, x[, 2])
irdec(y, x[, 3])

## discrete case ----

### example 1 ----
n <- 10000
s <- 0.1
x1 <- c(rep(0,n*s), runif(n*(1-s)))
x2 <- runif(n)
y <- x1

irdec(y,x1,dist.type.X = "discrete")
irdec(y,x2)

### example 2 ----
n <- 10000
x1 <- runif(n)

y1 <- rbinom(n,1,0.5)

y2 <- as.numeric(x1 >= 0.5)

irdec(y1,x1,dist.type.X = "discrete")
irdec(y2,x1,dist.type.X = "discrete")

FOCI::codec(y1,x1)
FOCI::codec(y2,x1)

### example 3 ----
library(ggplot2)

r_hurdle_poisson <- function(n, p_zero = 0.3, lambda = 2) {
  # n: number of samples
  # p_zero: probability of zero (Bernoulli part)
  # lambda: Poisson parameter for the positive counts

  # Step 1: Decide which values are zero (Bernoulli trial)
  is_zero <- rbinom(n, size = 1, prob = p_zero)

  # Step 2: For non-zero outcomes, draw from zero-truncated Poisson
  rztpois <- function(m, lambda) {
    # Rejection sampling for zero-truncated Poisson
    samples <- numeric(m)
    for (i in 1:m) {
      repeat {
        x <- rpois(1, lambda)
        if (x > 0) {
          samples[i] <- x
          break
        }
      }
    }
    return(samples)
  }

  # Fill in the result
  result <- numeric(n)
  num_nonzero <- sum(is_zero == 0)
  result[is_zero == 0] <- rztpois(num_nonzero, lambda)
  return(result)
}

# Generate samples
set.seed(123)
n <- 1000
p_zero <- 0.4
lambda <- 10

hurdle <- r_hurdle_poisson(n = n, p_zero = p_zero, lambda = lambda)
gamma_mix <- c(rep(0, round(p_zero * n)), rgamma(round((1 - p_zero) * n), shape = lambda, rate = 1))

# Combine into one data frame
df <- data.frame(
  value = c(hurdle, gamma_mix),
  source = rep(c("Hurdle Poisson", "Gamma Mixture"), each = n)
)

# Plot
ggplot(df, aes(x = value, fill = source)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 40) +
  labs(title = "Comparison: Hurdle Poisson vs Gamma Mixture",
       x = "Value", y = "Count", fill = "Distribution") +
  theme_bw()

x1 <- sort(gamma_mix)

y1 <- rbinom(n,1,0.5)

y2 <- sort(hurdle)



irdec(y1,x1,dist.type.X = "discrete")
irdec(y2,x1,dist.type.X = "discrete")

FOCI::codec(y1,x1)
FOCI::codec(y2,x1)

### example 4 ----

x1 <- sort(hurdle)

y1 <- rbinom(n,1,0.5)

y2 <- sort(gamma_mix)


irdec(y1,x1,dist.type.X = "discrete")
irdec(y2,x1,dist.type.X = "discrete")

FOCI::codec(y1,x1)
FOCI::codec(y2,x1)


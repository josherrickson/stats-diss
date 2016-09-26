rm(list = ls())
n <- 1000
d <- data.frame(x = rep(0:1, each = n/2),
                t = rep(0:1, times = n/2),
                m = rep(1:40, each = n/40)[1:n])

additive <- TRUE

reps <- 1000
results <- vector(length = reps)

for (i in 1:reps) {
  if (additive) {
    d$y <- rbinom(n, 1, .05 + d$x*.45 + d$t*.05)
    # treatment effect is .05, linear
  } else {
    d$y <- rbinom(n, 1, .05 + d$x*.45 + d$t*.05 + d$t*d$x*.13)
    # treatment effect is ~2.1 on the multiplicative scale
  }

  loss <- function(real, fitted) {
    # squared error
    #max((real - fitted)^2)
    # Fail.

    # log loss
    # mean(-real*log(fitted) + (1-real)*log(1-fitted), na.rm=TRUE)
    # mean log loss, 70% correct in Additive, and 95% correct in
    # Multiplicative.

    # logsitic (different?)
    mean(log(1 + exp(-fitted*real)))

    # boosting loss
    #mean(real * sqrt((1 - fitted)/fitted) + (1 - real) * sqrt(fitted/(1 - fitted)), na.rm=TRUE)
    # mean boosting around the same as log-loss

    # misclassification loss
    #mean(real*(fitted <= .05) + (1 - real)*(fitted > .05))
    # ~50% for both, as it takes only two values.
    # Fail

    # Logistic loss
    # mean(real*log(1 + exp(-2*fitted)) + (1 - real*log(1 + exp(-2*(1 - fitted)))))
    # ~50% for additive, 100% for multiplicative

    # hinge loss
    # mean(real * (1 - fitted) + (1 - real) * fitted)
    # ~50% for additive, 100% for multiplicative
  }


  m1 <- lm(y ~ t + x, data = d)
  m2 <- glm(y ~ t + x, data = d, family = binomial)

  #m1fitted <- pmin(pmax(fitted(m1), 0), 1)
  m1fitted <- fitted(m1)
  m1fitted[m1fitted >= 1] <- pmax(.999, 1 - (1 - max(m1fitted[m1fitted < 1])/2))
  m1fitted[m1fitted <= 0] <- pmin(.001, min(m1fitted[m1fitted > 0])/2)
  #if (any(m1fitted == 1 | m1fitted == 0))
  #  print(sum(m1fitted == 1 | m1fitted == 0))

  m1diff <- loss(d$y, m1fitted)
  m2diff <- loss(d$y, fitted(m2))

  results[i] <- ifelse(m1diff < m2diff, "additive", "multiplicative")

}
a_m <- ifelse(additive, "additive", "multiplicative")
table(results)[a_m]

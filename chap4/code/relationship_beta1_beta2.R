rm(list=ls())
x <- seq(.01, .99, by = .01)

beta2 <- seq(-10, 10, by=.1)
dd <- matrix(nrow = length(beta2), ncol = 2)
colnames(dd) <- c("beta2", "beta1")
for (i in seq_along(beta2)) {
  y <- arm::invlogit(car::logit(x) + beta2[i])
  dd[i, ] <- c(beta2[i], coef(lm(y~x))[2])
}

plot(dd, type = 'l')
d <- dnorm(beta2, )
lines(1/max(d)*d ~ beta2, lty=2)

minmax <- function(data, min = 0, max = 1) {
  pmin(max, pmax(min, data))
}

n <- 1000

x <- seq(-4,4,length = n)
z <- sample(0:1, n, TRUE)

yc <- arm::invlogit(2*x)
plot(yc ~ x, type = 'l', ylim = c(0,1), xlim = c(min(x), max(x)))

# multiplicative:
yt <- arm::invlogit(.1 * z + car::logit(yc))
lines(yt[z == 1] ~ x[z == 1], col = 'blue')

#linear
yt <- minmax(.01*z + yc)
lines(yt[z == 1] ~ x[z == 1], col = 'red')


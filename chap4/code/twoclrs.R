library(survival)
minmax <- function(data, min = 0, max = 1) {
  pmin(max, pmax(min, data))
}
meanseloss <- function(real, fitted) {
  mean((real - fitted) ^ 2)
}

meanlogloss <- function(real, fitted) {
  mean(-real * log(fitted) - (1 - real) * log(1 - fitted), na.rm = TRUE)
}

meanexploss <- function(real, fitted) {
  mean(real * sqrt((1 - fitted)/fitted) + (1 - real) *
         sqrt(fitted/(1 - fitted)), na.rm = TRUE)
}

TRUTH <- "additive"
#TRUTH <- "multiplicative"

n <- 1000
d <- data.frame(x = rnorm(n),
                t = rep(0:1, times = n/2),
                s = factor(rep(1:(n/4), each = 4)))
pi_c <- arm::invlogit(.05 + d$x)

beta_add <- .2
beta_mult <- 3

reps <- 100
save <- vector(length=reps)
for (i in 1:reps) {
  # additive
  pi_obs <- pi_c + beta_add*d$t
  d$y_add <- rbinom(n, 1, minmax(pi_obs))

  # multiplicative
  pi_obs <- arm::invlogit(beta_mult*d$t - car::logit(pi_c))
  d$y_mult <- rbinom(n, 1, pi_obs)

  if (TRUTH == "additive") {
    d$y <- d$y_add
  } else {
    d$y <- d$y_mult
  }


  ybars <- aggregate(d$y, by = list(d$s), FUN=mean)
  delta <- 1/ybars[,2]*(1-ybars[,2])
  delta[!is.finite(delta)] <- 0
  delta <- rep(delta, each = 4)

  additive <- clogit(y ~ x + strata(s) + t, data = d)
  multipli <- clogit(y ~ x + strata(s) + I(t*delta), data = d)

  fadd <- predict(additive, type='risk')
  fadd <- fadd/(1+fadd)

  fmult <- predict(multipli, type='risk')
  fmult <- fmult/(1+fmult)

   save[i] <- ifelse(logLik(additive) < logLik(multipli),
  #save[i] <- ifelse(summary(additive)$rsq[1] < summary(multipli)$rsq[1] ,
  # save[i] <- ifelse(meanlogloss(d$y, fadd) < meanlogloss(d$y, fmult),
  # save[i] <- ifelse(meanexploss(d$y, fadd) < meanexploss(d$y, fmult),
  # save[i] <- ifelse(meanseloss(d$y, fadd) < meanseloss(d$y, fmult),
                    "additive",
                    "multiplicative")
}
table(save)[TRUTH]

cachefile <- "cache/twostagesim.Rdata"
if (!file.exists(cachefile)) {
  meanseloss <- function(real, fitted) {
    mean((real - fitted) ^ 2)
  }

  meanlogloss <- function(real, fitted) {
    mean(-real * log(fitted) - (1 - real) * log(1 - fitted), na.rm = TRUE)
  }

  minmax <- function(data, min = 0, max = 1) {
    pmin(max, pmax(min, data))
  }

  n <- 1000
  d <- data.frame(x = rnorm(n),
                  t = rep(0:1, times = n/2))

  reps <- 10000

  bigsave <- matrix(nrow = 10, ncol = 6)
  colnames(bigsave) <- c("beta additive",
                         "beta multiplicative",
                         "mean se additive",
                         "mean se multiplicative",
                         "mean log additive",
                         "mean log multiplicative")

  beta_additives <- seq(0.01, .4, length = 20)
  beta_multiplicatives <- seq(.25, 2, length = 20)

  for (j in 1:10) {

    results <- matrix(nrow = reps, ncol = 4)
    colnames(results) <- colnames(bigsave)[-(1:2)]
    beta_add <- beta_additives[j]
    beta_mult <- beta_multiplicatives[j]
    for (i in 1:reps) {
      pi_c <- arm::invlogit(.05 + d$x)

      # additive
      pi_obs <- pi_c + beta_add*d$t
      d$y_add <- rbinom(n, 1, minmax(pi_obs))

      # multiplicative
      pi_obs <- arm::invlogit(beta_mult*d$t - car::logit(pi_c))
      d$y_mult <- rbinom(n, 1, pi_obs)

      # Additive
      mod1 <- glm(y_add ~ x, data = d, family = binomial, subset = (t == 0))

      m1 <- lm(y_add ~ 0 + t + offset(predict(mod1, newdata = d,
                                      type = 'response')), data = d)
      m2 <- glm(y_add ~ 0 + t + offset(predict(mod1, newdata = d)),
                data = d, family = binomial)

      m1fitted <- fitted(m1)
      m1fitted[m1fitted >= 1] <- pmax(.999, 1 - (1 - max(m1fitted[m1fitted < 1])/2))
      m1fitted[m1fitted <= 0] <- pmin(.001, min(m1fitted[m1fitted > 0])/2)

      results[i, 1] <- meanseloss(d$y_add, m1fitted)    < meanseloss(d$y_add, fitted(m2))
      results[i, 3] <- meanlogloss(d$y_add, m1fitted)   < meanlogloss(d$y_add, fitted(m2))

      # multiplicative
      mod1 <- glm(y_mult ~ x, data = d, family = binomial, subset = (t == 0))

      m1 <- lm(y_mult ~ t + offset(predict(mod1, newdata = d,
                                      type = 'response')), data = d)
      m2 <- glm(y_mult ~ t + offset(predict(mod1, newdata = d)),
                data = d, family = binomial)

      #m1fitted <- pmin(pmax(fitted(m1), 0), 1)
      m1fitted <- fitted(m1)
      m1fitted[m1fitted >= 1] <- pmax(.999, 1 - (1 - max(m1fitted[m1fitted < 1])/2))
      m1fitted[m1fitted <= 0] <- pmin(.001, min(m1fitted[m1fitted > 0])/2)
      #if (any(m1fitted == 1 | m1fitted == 0))
      #  print(sum(m1fitted == 1 | m1fitted == 0))

      results[i, 2]  <- meanseloss(d$y_mult, m1fitted)    > meanseloss(d$y_mult, fitted(m2))
      results[i, 4]  <- meanlogloss(d$y_mult, m1fitted)   > meanlogloss(d$y_mult, fitted(m2))


    }

    bigsave[j,] <- c(beta_add, beta_mult, apply(results, 2, sum))
  }

  bigsave[, -(1:2)] <- bigsave[, -(1:2)]/reps

  save(bigsave, file = cachefile)
}

load(cachefile)
rm(cachefile)

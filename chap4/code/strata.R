library(survival)
library(glmnet)
#cachefile <- "cache/strata.Rdata"
#if (!file.exists(cachefile)) {
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

  minmax <- function(data, min = 0, max = 1) {
    pmin(max, pmax(min, data))
  }

  n <- 1000
  d <- data.frame(x = rnorm(n),
                  t = rep(0:1, times = n/2),
                  s = factor(rep(1:(n/2), each = 2)))

  reps <- 10

  bigsave <- matrix(nrow = 10, ncol = 8)
  colnames(bigsave) <- c("beta additive",
                         "beta multiplicative",
                         "mean se additive",
                         "mean se multiplicative",
                         "mean log additive",
                         "mean log multiplicative",
                         "mean exp additive",
                         "mean exp multiplicative")

  beta_additives <- seq(0.01, .3, length = 10)
  beta_multiplicatives <- seq(.1, 5, length = 10)

  for (j in 1:10) {
    print(j)
    results <- matrix(nrow = reps, ncol = 6)
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

      run_with <- function(y) {
        d$y <- d[, y]

        # m1 - two stage, lm with strata as fixed
        # m2 - two stage, glm with no strata
        # m3 - one stage, CLR with strata
        # m4 - two stage, CLR with strata
        mod1 <- glm(y ~ x, data = d, family = binomial, subset = (t == 0))
#        m1 <- lm(y ~ t + s + offset(predict(mod1, newdata = d,
#                                                type = 'response')), data = d)
# Adding the strata made the linear model massively overfit
        m1 <- lm(y ~ t + offset(predict(mod1, newdata = d,
                                            type = 'response')), data = d)
        m2 <- glm(y ~ t + offset(predict(mod1, newdata = d)),
                  data = d, family = binomial)
        m3 <- clogit(y ~ t + x + strata(s), data = d)
        m4 <- clogit(y ~ t + offset(predict(mod1, newdata = d)) + strata(s),
                     data = d)
        xridge <- model.matrix(y ~ t + offset(predict(mod1, newdata = d, type = 'response')) + s, data = d)
        m5 <- cv.glmnet(xridge, d$y_add, alpha = 0, penalty.factor = c(0, 0, 0, rep(1, ncol(xridge) - 2)))


        m1fitted <- fitted(m1)
        m1fitted[m1fitted >= 1] <- pmax(.999, 1 - (1 - max(m1fitted[m1fitted < 1])/2))
        m1fitted[m1fitted <= 0] <- pmin(.001, min(m1fitted[m1fitted > 0])/2)

        m5fitted <- predict(m5, xridge, 0)
        m5fitted[m5fitted >= 1] <- pmax(.999, 1 - (1 - max(m5fitted[m5fitted < 1])/2))
        m5fitted[m5fitted <= 0] <- pmin(.001, min(m5fitted[m5fitted > 0])/2)

        fitted <- data.frame(m1 = m1fitted)
        fitted$m2 <- fitted(m2)
        fitted$m3 <- plogis(predict(m3, type = 'lp'))
        fitted$m4 <- plogis(predict(m4, type = 'lp'))
        fitted$m5 <- m5fitted

        return(fitted)
      }
      add <- run_with("y_add")
      mult <- run_with("y_mult")

      compare1 <- "m5"
      compare2 <- "m4"

      results[i, 1] <- meanseloss(d$y_add, add[,compare1])   < meanseloss(d$y_add, add[,compare2])
      results[i, 3] <- meanlogloss(d$y_add, add[,compare1])  < meanlogloss(d$y_add, add[,compare2])
      results[i, 5] <- meanexploss(d$y_add, add[,compare1])  < meanexploss(d$y_add, add[,compare2])

      results[i, 2]  <- meanseloss(d$y_mult, mult[,compare1])  > meanseloss(d$y_mult, mult[,compare2])
      results[i, 4]  <- meanlogloss(d$y_mult, mult[,compare1]) > meanlogloss(d$y_mult, mult[,compare2])
      results[i, 6]  <- meanexploss(d$y_mult, mult[,compare1]) > meanexploss(d$y_mult, mult[,compare2])


    }

    bigsave[j,] <- c(beta_add, beta_mult, apply(results, 2, sum))
  }

  bigsave[, -(1:2)] <- bigsave[, -(1:2)]/reps

#  save(bigsave, file = cachefile)
#}

#load(cachefile)
#rm(cachefile)

cachefile <- "cache/fixtausim.Rdata"
if (!file.exists(cachefile)) {
  n <- 100
  p <- 5
  pc <- .4
  informative <- .4
  true_t <- .5
  sigma2 <- 1
  true_inter <- -.5

  fixtausim.save <- matrix(nrow = 3, ncol = 9)
  seeds <- c(26, 28, 3345)

  for (ss in 1:3) {
    set.seed(seeds[ss])
    covs <- data.frame(matrix(rnorm(n*p), nrow = n))
    truebeta <- rep(0, p)
    truebeta[sample(1:p, round(informative*p))] <- rnorm(round(informative*p),0,1)

    treatment <- rep(0:1, c(n*pc, n*(1 - pc)))

    noise <- rnorm(n)
    yc_un <- as.matrix(covs) %*% truebeta
    yt_un <- yc_un + true_t*treatment + true_inter*treatment*yc_un
    resp <- ifelse(treatment == 1, yt_un, yc_un) + noise
    d <- data.frame(y = resp, covs)

    mod1 <- lm(y ~ ., data = d, subset = treatment == 0)

    object <- pbph(mod1, treatment, d)
    confint(object, returnShape = TRUE)



    bAndM <- pbph:::createBreadAndMeat(object, object$pbph$cluster)
    level <- .95

    bound <- 3

    tosolve <- function(eta) {
      corrected <- pbph:::corrVar(eta, object, bAndM, bound)[2,2]
      mod1 <- object$pbph$mod1
      df <- ifelse(is(mod1, "glm"), mod1$df.null, mod1$df)
      stat <- qt( (1 - level) / 2, df)
      return( (object$coef[2] - eta)^2 - stat^2 * corrected)
    }

    x <- -1:1
    t <- sapply(x,tosolve)
    dmat <- cbind(1, x, (x)^2)
    coefs_mid <- solve(t(dmat) %*% dmat) %*% t(dmat) %*% t

    x <- -28:-30
    t <- sapply(x,tosolve)
    dmat <- cbind(1, x, (x)^2)
    coefs_left <- solve(t(dmat) %*% dmat) %*% t(dmat) %*% t

    x <- 28:30
    t <- sapply(x,tosolve)
    dmat <- cbind(1, x, (x)^2)
    coefs_right <- solve(t(dmat) %*% dmat) %*% t(dmat) %*% t



    fixtausim.save[ss,] <- round(c(coefs_left, coefs_mid, coefs_right), 1)
  }

  save(fixtausim.save, file = cachefile)
}

load(cachefile)
rm(cachefile)

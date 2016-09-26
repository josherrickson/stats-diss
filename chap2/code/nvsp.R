cachefile <- "cache/nvsp.Rdata"
if (!file.exists(cachefile)) {
  pc <- .5
  true_t <- .5
  sigma2 <- 1
  ns <- c(50, 100, 250, 500, 1000)
  cs <- seq(1,10, by = .5)
  ps <- 1:1000
  reps <- 10000
  bigsave <- matrix(ncol = length(cs),
                    nrow = length(ns))
  colnames(bigsave) <- cs
  rownames(bigsave) <- ns
  psave <- bigsave

  for (j in seq_along(ns)) {
    n <- ns[j]
    print(n)
    for (k in seq_along(cs)) {
      c <- cs[k]
      print(c)
      p <- tail(which(ps^2*log(ps)^2/n < c),1)
      psave[j,k] <- p
      if (p > .2*n) next()
      ti <- runif(1,-1,2)
      save <- vector(length = reps)
      for (i in 1:reps) {
        covs <- data.frame(matrix(rnorm(n*p), nrow = n))
        truebeta <- rnorm(p)

        treatment <- rbinom(n, 1, 1 - pc)

        noise <- rnorm(n)
        yc_un <- as.matrix(covs) %*% truebeta
        yt_un <- yc_un + true_t*treatment + ti*treatment*yc_un
        resp <- ifelse(treatment == 1, yt_un, yc_un) + noise
        d <- data.frame(y = resp, covs)

        mod1 <- lm(y ~ ., data = d, subset = treatment == 0)

        e <- pbph(mod1, treatment, d)
        ci <- confint(e, "pred", returnShape = TRUE)
        type <- attr(ci, "shape")
        covered <- ci[1] < ti & ti < ci[2]
        if (type == "disjoint") {
          covered <- TRUE
        }

        save[i] <- covered

      }
      bigsave[j,k] <- mean(save)
    }
  }
  save(bigsave, file = cachefile)
}

load(cachefile)
rm(cachefile)

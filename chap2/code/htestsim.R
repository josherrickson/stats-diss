cachefile <- "cache/htestsim.Rdata"
if (!file.exists(cachefile)) {
  pc <- .4
  informative <- .4
  sigma2 <- 1

  reps <- 1000

  htestsim.save <- makeSaveMatrix(c("truth100",  "estimate100",  "pval100",
                                    "truth1000", "estimate1000", "pval1000"), reps=reps)

  for (n in c(100, 1000)) {
    which <- if(n == 100) 1:3 else 4:6
    p <- if (n == 100) 7 else 17
    for (i in seq_len(reps)) {
      ti <- 0
      tt <- rnorm(1)

      covs <- data.frame(matrix(rnorm(n*p), nrow=n))
      inf <- round(informative * p)
      truebeta <- c(rnorm(inf), rep(0, p - inf))
      treatment <- rbinom(n, 1, 1 - pc)

      noise <- rnorm(n)
      yc_un <- as.matrix(covs)%*%truebeta
      yt_un <- yc_un + tt*treatment + ti*treatment*yc_un
      resp <- ifelse(treatment==1, yt_un, yc_un) + noise
      d <- data.frame(y=resp, covs)

      mod1 <- lm(y ~ . , data=d, subset=treatment==0)

      sm1 <- summary(mod1)

      e <- pbph(mod1, treatment, d)
      sm <- summary(e)

      htestsim.save[i,which] <- c(ti, e$coef[2], sm$coef[2,4])

    }
  }
  htestsim.save <- data.frame(htestsim.save)
  save(htestsim.save, file = cachefile)
}

load(cachefile)
rm(cachefile)

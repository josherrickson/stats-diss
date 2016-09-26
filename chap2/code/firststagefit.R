p <- 7
pc <- .4
n <- 100

cachefile <- "cache/firststagefit.Rdata"
if (!file.exists(cachefile)) {
  informative <- .6
  sigma2 <- 1

  reps <- 1000

  firststagefit.save <- makeSaveMatrix(c("truth", "estimate", "tstat",
                                         "pval", "lb", "ub", "type",
                                         "covered", "F-stat", "R^2"),
                         reps=reps)

  for (i in seq_len(reps)) {
    ti <- runif(1, -1, 2)
    tt <- rnorm(1)

    covs <- data.frame(matrix(rnorm(n*p), nrow=n))
    truebeta <- c(rnorm(round(informative*p)),
                  rep(0, round((1-informative)*p)))
    treatment <- rbinom(n, 1, 1 - pc)

    noise <- rnorm(n)
    yc_un <- as.matrix(covs)%*%truebeta
    yt_un <- yc_un + tt*treatment + ti*treatment*yc_un
    resp <- ifelse(treatment==1, yt_un, yc_un) + noise
    d <- data.frame(y=resp, covs)

    mod1 <- lm(y ~ ., data=d[treatment==0,])
    sm1 <- summary(mod1)
    mod2 <- pbph(mod1, treatment, d)

    ci <- confint(mod2, "pred", returnShape = TRUE)

    if (attr(ci, "type") == "finite") {
      covered <- ci[1] < ti & ti < ci[2]
    } else if (attr(ci, "type") == "infinite") {
      covered <- TRUE
    } else if (attr(ci, "type") == "disjoint") {
      covered <- ti < ci[1] | ci[2] < ti
    } else {
      stop(paste("Problem:", attr(ci, "type")))
    }

    citype <- match(attr(ci, "type"), c("finite", "infinite", "disjoint"))


    firststagefit.save[i,] <- c(ti, mod2$coef[2], summary(mod2)$coef[2,3:4], ci,
                  citype, covered, sm1$fstat[1], sm1$r.s)


  }

  firststagefit.save <- as.data.frame(firststagefit.save)
  firststagefit.save$type <- as.factor(firststagefit.save$type)
  levels(firststagefit.save$type) <- c("finite", "infinite", "disjoint")

  save(firststagefit.save, file = cachefile)
}

load(cachefile)
rm(cachefile)

boxplot(log(firststagefit.save$"F-stat") ~ firststagefit.save$type,
        ylab="log(F)", xaxt="n",)
axis(1, 1:3, c("Finite", "Infinite", "Disjoint"))

abline(h=log(qf(.95, p, pc*n-p-1)), col='blue')
abline(h=log(qf(.999, p, pc*n-p-1)), col='red')
legend("topright", legend=c("95% Signif", "99.9% Signif"),col=c('blue', 'red'), lty=c(1,1))

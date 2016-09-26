cachefile <- "cache/coverage_by_type1000.Rdata"
if (!file.exists(cachefile)) {
  n <- 1000
  p <- 17
  pc <- .4
  informative <- .4
  true_t <- .5
  sigma2 <- 1
  true_inter <- seq(-1,2,by=.5)

  reps <- 1000
  cov.by.type1000.save <- makeSaveMatrix(c("truth", "estimate",
                                            "cont_un", "cont_cov",
                                            "disjoint_un",
                                            "disjoint_cov"),
                                          reps=length(true_inter))
  for (j in 1:length(true_inter)) {
    ti <- true_inter[j]
    save <- makeSaveMatrix(c("estimate","lb","ub", "type", "covered"),
                           reps=reps)
    for (i in 1:reps) {
      covs <- data.frame(matrix(rnorm(n*p), nrow=n))
      truebeta <- c(rnorm(round(informative*p)),
                    rep(0, round((1-informative)*p)))

      treatment <- rbinom(n, 1, 1 - pc)

      noise <- rnorm(n)
      yc_un <- as.matrix(covs)%*%truebeta
      yt_un <- yc_un + true_t*treatment + ti*treatment*yc_un
      resp <- ifelse(treatment==1, yt_un, yc_un) + noise
      d <- data.frame(y=resp, covs)

      mod1 <- lm(y ~ ., data=d[treatment==0,])
      mod2 <- pbph(mod1, treatment, d)

      ci <- confint(mod2, "pred", returnShape=TRUE)
      type <- attr(ci, "type")
      if (type %in% c("finite", "infinite")) {
        covered <- ci[1] < ti & ti < ci[2]
      } else if (type == "disjoint") {
        covered <- ti < ci[1] | ci[2] < ti
      } else {
        stop(paste("Problem:", type))
      }
      type <- switch(type,
                     finite=1,
                     infinite=1,
                     disjoint=2)

      save[i,] <- c(mod2$coef[2], ci[1], ci[2], type, covered)
    }
    save <- data.frame(save)
    save$type <- as.factor(save$type)
    levels(save$type) <- 1:2
    save$covered <- factor(save$covered, levels=0:1)

    #bytype <- aggregate(save$covered, by=list(save$type), FUN=mean)[,2]
    overall <- mean(save$covered==1)

    cov.by.type1000.save[j,] <- c(ti, mean(save$estimate),
                     as.vector(table(save$covered, save$type)))

  }

  cov.by.type1000.save <- as.data.frame(cov.by.type1000.save)
  cov.by.type1000.save$overall_un <- apply(cov.by.type1000.save[names(cov.by.type1000.save)[grepl("_un", names(cov.by.type1000.save))]], 1, sum)
  cov.by.type1000.save$overall_cov <- apply(cov.by.type1000.save[names(cov.by.type1000.save)[grepl("_cov", names(cov.by.type1000.save))]], 1, sum)

  cov.by.type1000.save$overall_per <- cov.by.type1000.save$overall_cov/(cov.by.type1000.save$overall_un + cov.by.type1000.save$overall_cov)
  cov.by.type1000.save$cont_per <- cov.by.type1000.save$cont_cov/(cov.by.type1000.save$cont_un + cov.by.type1000.save$cont_cov)
  cov.by.type1000.save$disjoint_per <- cov.by.type1000.save$disjoint_cov/(cov.by.type1000.save$disjoint_un + cov.by.type1000.save$disjoint_cov)

  cov.by.type1000.save <- cov.by.type1000.save[,
                         c(1, 2, 7:9, 3:4, 10, 5:6, 11)]

  save(cov.by.type1000.save, file = cachefile)
}

load(cachefile)
rm(cachefile)

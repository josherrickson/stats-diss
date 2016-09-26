cachefile <- "cache/clustersim_1000_10.Rdata"
if (!file.exists(cachefile)) {
  n <- 1000
  p <- 13
  C <- 10
  pc <- .4
  informative <- .4
  true_t <- .5
  sigma2 <- 1
  true_inter <- seq(-1,2,by = .5)

  reps <- 1000
  clus_1000_10 <- pbph:::makeSaveMatrix(c("truth", "estimate", "overall_un",
                                          "overall_cov", "cont_un", "cont_cov",
                                          "disjoint_un", "disjoint_cov"),
                                        reps = length(true_inter))
  for (j in 1:length(true_inter)) {
    ti <- true_inter[j]
    print(ti)
    save <- pbph:::makeSaveMatrix(c("estimate", "lb", "ub", "type", "covered"),
                                 reps)
    for (i in 1:reps) {
      covs <- data.frame(matrix(rnorm(n*p), nrow = n))
      truebeta <- rep(0, p)
      truebeta[sample(1:p, round(informative*p))] <- rnorm(round(informative*p),0,1)

      cluster <- sample(1:C, n, TRUE)
      tbyc <- rbinom(C, 1, 1 - pc)
      while (sum(tbyc == 1) <= 1 | sum(tbyc == 0) <= 1)
        tbyc <- rbinom(C, 1, 1 - pc)
      treatment <- tbyc[cluster]


      noise <- rnorm(n)
      yc_un <- as.matrix(covs) %*% truebeta
      yt_un <- yc_un + true_t*treatment + ti*treatment*yc_un
      resp <- ifelse(treatment == 1, yt_un, yc_un) + noise
      d <- data.frame(y = resp, covs)

      mod1 <- lm(y ~ ., data = d, subset = treatment == 0)

      e <- pbph(mod1, treatment, d, cluster = cluster)
      ci <- confint(e, "pred", returnShape = TRUE)
      type <- attr(ci, "type")
      covered <- ci[1] < ti & ti < ci[2]
      if (type == "disjoint") {
        covered <- !covered
      }
      type <- switch(type,
                     finite = 1,
                     infinite = 1,
                     disjoint = 2)

      save[i,] <- c(e$coef[2], ci[1], ci[2], type, covered)

    }
    save <- data.frame(save)
    save$type <- as.factor(save$type)
    levels(save$type) <- 1:2
    save$covered <- factor(save$covered, levels = 0:1)

    overall <- mean(save$covered == 1)

    clus_1000_10[j,] <- c(ti, mean(save$estimate), table(save$covered),
                     as.vector(table(save$covered, save$type)))

  }

  clus_1000_10 <- as.data.frame(clus_1000_10)
  clus_1000_10$perc_cov <- clus_1000_10$overall_cov/reps*100
  save(clus_1000_10, file = cachefile)
}


load(cachefile)
rm(cachefile)

cachefile <- "cache/overunder.Rdata"
reps <- 1000
if (!file.exists(cachefile)) {
  n <- 100
  pc <- .4
  p <- 7


  saveCols <- c("truth", "estimate", "lb", "ub", "type", "covered")

  saveoracle <- makeSaveMatrix(saveCols, reps=reps)
  saveunder <- makeSaveMatrix(saveCols, reps=reps)
  saveover <- makeSaveMatrix(saveCols, reps=reps)

  for (i in seq_len(reps)) {
    ti <- runif(1,-1,2)
    tt <- rnorm(1)

    covs <- data.frame(matrix(rnorm(n*p), nrow=n))
    truebeta <- c(rnorm(4), rep(0, 3))
    treatment <- rbinom(n, 1, 1 - pc)

    noise <- rnorm(n)
    yc_un <- as.matrix(covs)%*%truebeta
    yt_un <- yc_un + tt*treatment + ti*treatment*yc_un
    resp <- ifelse(treatment==1, yt_un, yc_un) + noise
    d <- data.frame(y=resp, covs)

    mod1under  <- lm(y ~ X1               , data=d, subset=treatment==0)
    mod1oracle <- lm(y ~ X1 + X2 + X3 + X4, data=d, subset=treatment==0)
    mod1over   <- lm(y ~ .                , data=d, subset=treatment==0)

    mod2under  <- pbph(mod1under, treatment, d)
    mod2oracle <- pbph(mod1oracle, treatment, d)
    mod2over   <- pbph(mod1over, treatment, d)

    ciunder <- confint(mod2under, "pred", returnShape=TRUE)
    typeunder <- attr(ciunder, "type")
    if (typeunder %in% c("finite", "infinite")) {
      coveredunder <- ciunder[1] < ti & ti < ciunder[2]
    } else if (typeunder == "disjoint") {
      coveredunder <- ti < ciunder[1] | ciunder[2] < ti
    } else {
      stop(paste("Problem:", typeunder))
    }
    typeunder <- switch(typeunder,
                        finite=1,
                        infinite=1,
                        disjoint=2)


    cioracle <- confint(mod2oracle, "pred", returnShape=TRUE)
    typeoracle <- attr(cioracle, "type")
    if (typeoracle %in% c("finite", "infinite")) {
      coveredoracle <- cioracle[1] < ti & ti < cioracle[2]
    } else if (typeoracle == "disjoint") {
      coveredoracle <- ti < cioracle[1] | cioracle[2] < ti
    } else {
      stop(paste("Problem:", typeoracle))
    }
    typeoracle <- switch(typeoracle,
                         finite=1,
                         infinite=1,
                         disjoint=2)


    ciover <- confint(mod2over, "pred", returnShape=TRUE)
    typeover <- attr(ciover, "type")
    if (typeover %in% c("finite", "infinite")) {
      coveredover <- ciover[1] < ti & ti < ciover[2]
    } else if (typeover == "disjoint") {
      coveredover <- ti < ciover[1] | ciover[2] < ti
    } else {
      stop(paste("Problem:", typeover))
    }
    typeover <- switch(typeover,
                       finite=1,
                       infinite=1,
                       disjoint=2)

    saveunder[i,] <- c(ti, mod2under$coef[2], ciunder,
                       typeunder, coveredunder)
    saveoracle[i,] <- c(ti, mod2oracle$coef[2], cioracle,
                        typeoracle, coveredoracle)
    saveover[i,] <- c(ti, mod2over$coef[2], ciover,
                      typeover, coveredover)

  }

  saveunder <- as.data.frame(saveunder)
  saveunder$type <- as.factor(saveunder$type)
  levels(saveunder$type) <- c("Finite", "Disjoint")
  saveoracle <- as.data.frame(saveoracle)
  saveoracle$type <- as.factor(saveoracle$type)
  levels(saveoracle$type) <- c("Finite", "Disjoint")
  saveover <- as.data.frame(saveover)
  saveover$type <- as.factor(saveover$type)
  levels(saveover$type) <- c("Finite", "Disjoint")

  overunder.save <- list(saveunder=saveunder,
                         saveoracle=saveoracle,
                         saveover=saveover)

  save(overunder.save, file = cachefile)
}

load(cachefile)
rm(cachefile)

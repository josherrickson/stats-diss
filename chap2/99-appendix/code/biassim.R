cachefile <- "cache/biassim.Rdata"
if (!file.exists(cachefile)) {

  for ( n in c(100, 1000)) {
    p <- ifelse(n == 100, 7, 17)
    pc <- .5
    informative <- .4
    true_t <- .5
    sigma2 <- 1

    reps <- 1000

    truths <- seq(-1,2, by=.1)

    bigsave <- matrix(nrow=length(truths), ncol=4)
    for (j in 1:length(truths)) {
     # print(paste("Starting",truths[j]))

      save <- matrix(nrow=reps, ncol=4)
      true_inter <- truths[j]

      for(i in seq_len(reps)) {
        x <- matrix(rnorm(n*p), nrow=n)
        truebeta <- rep(0, p)
        truebeta[sample(1:p, round(informative*p))] <- rnorm(round(informative*p))

        t <- rep(0:1, c(n*pc, n*(1-pc)))

        noise <- rnorm(n)
        yc_un<- x%*%truebeta
        yt_un<- yc_un+ true_t*t + true_inter*t*yc_un
        y <- ifelse(t==1, yt_un, yc_un) + noise

        d <- data.frame(cbind(y,x))
        colnames(d) <- c("y", paste("x", 1:p, sep=''))

        # Stage 1:
        mod1 <- lm(y ~ ., data=d[t == 0,])

        d$pred <- predict(mod1, d)

        # Stage 2:
        mod2 <- lm(y - pred ~ pred, data=d[t == 1,])

        # naive est
        predtr <- d$pred[t == 1]
        predct <- d$pred[t == 0]

        bread.naive <- crossprod(cbind(rep(1,length(predtr)), predtr))
        meat.naive <- crossprod(cbind(mod2$res, mod2$res*predtr))

        naive <- sqrt((solve(bread.naive)%*%meat.naive%*%solve(bread.naive))[2,2])

        # corrected

        xt <- d[t == 1,]
        xt <- xt[,names(xt) %in% names(mod1$coef)]
        xt <- cbind(rep(1,nrow(xt)), xt)
        yt <- d$y[t == 1]

        xc <- d[t == 0,]
        xc <- xc[,names(xc) %in% names(mod1$coef)]
        xc <- cbind(rep(1,nrow(xc)), xc)
        yc <- d$y[t == 0]

        bread.22 <- bread.naive


        # model based
        nullres <- d[t==1,]$y - d[t==1,]$pred - mod2$coef[1]
        meat.22.mb <- crossprod(cbind(nullres, nullres*predtr))

        # empirical
        meat.22.emp <- meat.naive

        # model-based
        bread.21.mb <- rbind(apply(-xt, 2, sum),
                             apply((yt - mod2$coef[1] -
                                      2*predtr)*xt, 2, sum))
        # empirical
        bread.21.emp <- rbind(apply(-(1+mod2$coef[2])*xt, 2, sum),
                              apply((yt - mod2$coef[1] -
                                       2*(1 + mod2$coef[2])*predtr)*xt, 2, sum))

        meat.22 <- meat.22.emp
        bread.21 <- bread.21.emp

        bread.11 <- t(as.matrix(xc))%*%as.matrix(xc)

        meat.11 <- t(mod1$res*as.matrix(xc))%*%(mod1$res*as.matrix(xc))

        corrected <- sqrt((solve(bread.22)%*%
                             (meat.22 +
                                bread.21%*%solve(bread.11)%*%meat.11%*%
                                solve(bread.11)%*%t(bread.21))%*%
                             solve(bread.22))[2,2])

        coef <- mod2$coef[2]

        # bias
        Sigmac <- t(t(x[t==0,])%*%x[t==0,])/sum(t==0)
        Sigmat <- t(t(x[t==1,])%*%x[t==1,])/sum(t==1)

        bias <- summary(mod2)$sigma^2*sum(solve(Sigmac)*Sigmat)/(sum(t==0)-1)

        coefcor <- with(d[t==1,],(cov(y - pred, pred) + bias)/(var(pred) + bias))

        save[i,] <- c(coef, coefcor, naive, corrected)
      }

      tstat <- qt(.975, mod1$df)


      bigsave[j,] <- c(true_inter,
                       sum(abs(save[,1]-true_inter) -
                             tstat*save[,3] <= 0)/reps,
                       sum(abs(save[,1]-true_inter) -
                             tstat*save[,4] <= 0)/reps,
                       sum(abs(save[,2]-true_inter) -
                             tstat*save[,4] <= 0)/reps)
    }

    bigsave <- data.frame(bigsave)
    colnames(bigsave) <- c("true", "none", "se", "bias_se")
    library(ggplot2)

    assign(paste0("bs_",n), reshape(bigsave,varying = c("none", "se", "bias_se"),
                                    v.names = "coverage", direction = "long",
                                    times = c("none", "se", "bias_se"),
                                    new.row.names = 1:(5*nrow(bigsave))))

  }

  biassim <- list(bs_100, bs_1000)

  save(biassim, file = cachefile)
}

load(cachefile)
rm(cachefile)

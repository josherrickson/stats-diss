# Some of this stuff is used for plotting, so define it now
n <- 100
t <- rep(0:1, each=50)
p <- 40
p <- floor(p/2)*2 # make sure its even

reps <- 100

cachefile <- "cache/bias_r2.Rdata"
if (!file.exists(cachefile)) {
  bias_r2.save<- matrix(nrow=reps*p, ncol=9)
  for (i in seq_len(reps)) {
    x <- matrix(rnorm(n*p), nrow=n)
    beta <- c(rnorm(p/2), rep(0,p/2))
    y <- x%*%beta + rnorm(n) # no treatment effect
    d <- data.frame(y,x,t)

    for (j in seq_len(p)) {
      mod1 <- lm(y ~ . + 0, data=d[d$t==0, 1:(j+1), drop=FALSE])
      d$pred <- predict(mod1, newdata=d)
      mod2 <- lm(y - pred ~ pred, data=d[d$t==1,])

      bias_r2.save[(i-1)*p+j,] <- c(j, summary(mod1)$r.s, mod2$coef[2], summary(mod1)$fstat, summary(mod2)$fstat)
    }
  }

  bias_r2.save<- as.data.frame(bias_r2.save)
  names(bias_r2.save) <- c("p", "mod1r2", "coef", "fstatmod1", "fdf1mod1", "fdf2mod1", "fstatmod2", "fdf1mod2", "fdf2mod2")

  save(bias_r2.save, file = cachefile)
}
load(cachefile)
rm(cachefile)

par(mfrow=c(2,1), oma=c(0,0,0,0), mar=c(4.1, 4.1, 1, 1))
plot(mean(mod1r2~p, data=bias_r2.save), ylim=c(0,1), type="l", lty=2, ylab=parse(text="R^2"), xlab="")
abline(v=20, lty=2, col='lightgrey')

plot(mean(coef~p, data=bias_r2.save), ylim=c(-.5,0), pch=15, col='white',
     xlab="Number of included variables", ylab="Bias Confidence Bounds")

upper <- mean(coef~p, data=bias_r2.save) + 1.96*sd(coef~p, data=bias_r2.save)/sqrt(reps)
lower <- mean(coef~p, data=bias_r2.save) - 1.96*sd(coef~p, data=bias_r2.save)/sqrt(reps)
polygon(c(1:p, p:1), c(lower, rev(upper)), col=rgb(1,0,0,.5), border=rgb(1,0,0,.5))
abline(h=0, lty=3, col='lightgrey')
abline(v=20, lty=2, col='lightgrey')

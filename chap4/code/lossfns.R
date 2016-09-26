logloss <- function(y, yhat) {
  -y*log(yhat) - (1 - y)*log(1-yhat)
}
plot(function(x) logloss(1, x), 0, 1, ylim=c(0, 5), xlim=c(-1,1), col=2)
plot(function(x) logloss(0, x), 0, 1, add=TRUE, lty=2, col=2)

logisticloss <- function(y, yhat) {
  log(1 + exp(-y*yhat))
}
plot(function(x) logisticloss(1, x), -1, 1, add=TRUE, col=3)
plot(function(x) logisticloss(-1, x), -1, 1, add=TRUE, col=3, lty=2)

boostloss <- function(y, yhat) {
  y * sqrt((1 - yhat)/yhat) + (1 - y) * sqrt(yhat/(1 - yhat))
}
plot(function(x) boostloss(1, x), 0, 1, add=TRUE, col=4)
plot(function(x) boostloss(0, x), 0, 1, add=TRUE, col=4, lty=2)

seloss <- function(y, yhat) {
  (y - yhat)^2
}
plot(function(x) seloss(1, x), 0, 1, add=TRUE, col=5)
plot(function(x) seloss(0, x), 0, 1, add=TRUE, col=5, lty=2)

hingeloss <- function(y, yhat) {
  y * (1 - yhat) + (1 - y) * yhat
}
plot(function(x) hingeloss(1, x), 0, 1, add=TRUE, col=6)
plot(function(x) hingeloss(0, x), 0, 1, add=TRUE, col=6, lty=2)

legend("topleft", legend = c("Y = 0", "Y = 1", "log-loss",
                             "logistic", "boost", "squared error",
                             "hinge"),
       lty = c(1,2,1,1,1,1,1),
       col = c(1,1,2,3,4,5,6))

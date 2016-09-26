cachefile <- "cache/choose_version_of_var.Rdata"
if (!file.exists(cachefile)) {
  n <- 100
  p <- 5
  pc <- .4
  informative <- .4
  true_t <- .5
  sigma2 <- 1
  true_inter <- seq(-1,2,by = .5)

  reps <- 1000

  ## 1) theta_0 in all (model-based)
  ## 2) \hat\theta in all (empirical)
  ## 3) Bread uses hat, meat uses null
  ## 4) Bread uses null, meat uses hat

  ############# Make bread use etahat instead of eta0
  changeBreadToUseHat <- function() {
    bread21 <- function(model, eta) {
      mod1 <- model$pbph$mod1
      data <- model$pbph$data

      treatment <- data$treatment

      resp <- eval(formula(mod1)[[2]], envir = data)[treatment == 1]
      covs <- model.matrix(formula(mod1), data = data)
      covs <- pbph:::addIntercept(covs)[treatment == 1, , drop = FALSE]
      pred <- data$pred[treatment == 1]

      scale <- pbph:::glmScale(mod1, newdata = as.data.frame(covs))

      eta <- model$coef[2]
      tau <- model$coef[1]
      # Replace tauhat with tauhat_eta0
      # When eta = eta_0, the RHS is constant, so we just need the mean.
      tau <- mean(resp - (1 + eta) * pred)

      b21.1 <- apply(-(1 + eta) * covs * scale, 2, sum)
      b21.2 <- apply( (resp - tau - 2 * (1 + eta) * pred) * covs * scale, 2, sum)

      rbind(b21.1, b21.2)
    }
    assignInNamespace("bread21", bread21, pos = "package:pbph")
  }

  ############### Make meat use null instead of hat
  changeMeatToUseNull <- function() {
    meat22.null <- function(model, eta) {
      newres <- lm(y - (1 + eta)*pred ~ 1, subset = treatment == 1, data = model$pbph$data)$res
      x <- model.matrix(model)
      crossprod(newres * x, newres * x)
    }

    corrVar <- function(eta, object,
                        breadAndMeat = pbph:::createBreadAndMeat(object, cluster = object$pbph$cluster)) {

      b21 <- pbph:::bread21(object, eta)
      m22 <- meat22.null(object, eta = eta)

      corrected <- pbph:::correctedvar(breadAndMeat$b11,
                                       b21,
                                       breadAndMeat$b22,
                                       breadAndMeat$m11,
                                       m22)

      return(corrected)
    }

    assignInNamespace("corrVar", corrVar, pos = "package:pbph")
  }

  chooseVarSave <- pbph:::makeSaveMatrix(c("form", "truth", "covered", "mean.width", "median.width"),
                                   reps = 4 * length(true_inter))

  runSim <- function(form) {
    for (j in 1:length(true_inter)) {
      ti <- true_inter[j]
      save <- vector(length = reps)
      width <- vector(length = reps)
      for (i in 1:reps) {
        covs <- data.frame(matrix(rnorm(n*p), nrow = n))
        truebeta <- rep(0, p)
        truebeta[sample(1:p, round(informative*p))] <- rnorm(round(informative*p),0,1)

        treatment <- rep(0:1, c(n*pc, n*(1 - pc)))

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
          covered <- !covered
        }
        width[i] <- diff(ci[1,])
        save[i] <- covered

      }
      # Using <<- to save to global chooseVarSave
      chooseVarSave[(form - 1)*length(true_inter) + j,] <<- c(form, ti, mean(save), mean(width[is.finite(width)]), median(width))
    }
  }


  # 4) Hybrid 2 (real)
  runSim(4)

  # 1) model based
  changeMeatToUseNull()
  runSim(1)

  # 3) Hybrid 1
  changeBreadToUseHat()
  changeMeatToUseNull()
  runSim(3)

  devtools::load_all()

  # 2) empirical
  changeBreadToUseHat()
  runSim(2)


  chooseVarSave <- as.data.frame(chooseVarSave)
  save(chooseVarSave, file = cachefile)
}

load(cachefile)
rm(cachefile)

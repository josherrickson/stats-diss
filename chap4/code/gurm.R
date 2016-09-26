cachefile <- "cache/gurm.Rdata"
if (!file.exists(cachefile)) {
  # Run a modified version of Carrie's script to clean the data as she
  # did and to merge the match from `matching.rda` in.
  wd <- getwd()
  setwd("/Volumes/data/Gurm-pci2007/")
  source("/Volumes/data/Gurm-pci2007/carriescript.R")
  setwd(wd)
  # Carrie's script saves a bunch of temp variables, all we care about
  # is `md2`, so nuke the rest to keep a clean env.
  rm(list = ls()[!grepl("md2", ls())])

  # Carrie's script made `mi_stelv` a factor inside the formulas, lets
  # just do it beforehand (mlm errors with objects like
  # `as.numeric(mi_stelv==1)` in the formula).
  md2$mi_stelv <- as.factor(md2$mi_stelv)

  # Match created pre-subproblem, needed for summary (for debugging)
  attr(md2$match, "subproblem") <- factor(rep(1, nrow(md2)))


  # drop missing data
  md2 <- md2[complete.cases(md2),]

  m1 <- glm(ih_v_c ~ # vascular complications
              gfr + # glomerular filtrationrate
              mi_stelv + # ST-segment elevation myocardial infarction
              hx_chf + # prior congestive heart failure
              hx_chf_a + # congestive heart failure on admission
              hx_ptca + # prior PCI?
              hosp_cd,
            data = md2,
            family = binomial,
            subset = pr_vcd == 0)

  ybars <- aggregate(md2$ih_v_c, by = list(md2$match), FUN = mean)
  names(ybars) <- c("match", "ybars")
  md2 <- merge(md2, ybars, by = "match")

  # 'link' gives the fitted value on the logit scale. Model 2 is
  # logit(y) - logit(y_c) = xbeta.
  p <- predict(m1, type = 'link', newdata = md2)

  m2log <- glm(ih_v_c ~ pr_vcd + offset(p),
               data = md2,
               family = binomial)

  m2lin <- glm(ih_v_c ~ pr_vcd:ybars + offset(p),
               data = md2,
               family = binomial)


  # http://stackoverflow.com/questions/35329585/how-to-get-fitted-values-from-clogit-model
  p1 <- predict(m2log, type = "response", newdata = md2)
  p2 <- predict(m2lin, type = "response", newdata = md2)


  pred <-  predict(m1, newdata = md2, type = "response")

  w <- ipr(md2$pr_vcd, md2$match)
  w <- w/sum(w)

  stage2 <- lm(ih_v_c~ 1 + offset(pred), data = md2,
               weights = w,
               subset = pr_vcd == 1)

  stage2 <- lm(ih_v_c ~ 0 + pr_vcd + offset(pred), data = md2,
               weights = w)

  b11 <- pbph:::bread11(m1) # No weights, already inverted
  b21 <- matrix(apply(w[md2$pr_vcd == 1]*model.matrix(formula(m1), data = md2)[md2$pr_vcd == 1,]*pred[md2$pr_vcd == 1]*(1 - pred[md2$pr_vcd == 1]), 2, sum), nrow = 1)
  b22 <- .5
  m11 <- pbph:::meat11(m1) # No weights
  m22 <- pbph:::meat22(stage2) # sum(w[md2$pr_vcd == 1]^2 * residuals(stage2)^2)

  # coef is -.0020288

  gurm.save <- c(mean(-md2$ih_v_c*log(p1) - (1 - md2$ih_v_c)*log(1 - p1), na.rm = TRUE),
                 mean(-md2$ih_v_c*log(p2) - (1 - md2$ih_v_c)*log(1 - p2), na.rm = TRUE),
                 stage2$coef,
                 sqrt(4*(m22 + b21 %*% b11 %*% m11 %*% t(b11) %*% t(b21)))
  )

  save(gurm.save, file = cachefile)
}

load(cachefile)
rm(cachefile)

cachefile <- "cache/gine_naive_coverage.Rdata"
if(!file.exists(cachefile)) {
  final <- read.dta("~/Dropbox/research/ben/interaction/data/AER_2010_1090_R2/data/final.dta")

  # Not sure what this means precisely, but it splits almost all of the results
  # missing treatment into NA
  final$source_oct_repay <- final$source_oct_baseline * final$source_repayment_nov

  # Create indicators to agree with how STATA handles the collinearity
  for(i in c(2,3,6:11,13))  {
    final <- cbind(final, 1*(final$focode == i))
    names(final)[length(names(final))] <- paste('focode',i,sep='')  }
  for(i in 35:39)  {
    final <- cbind(final, 1*(final$offer_week == i))
    names(final)[length(names(final))] <- paste('offerweek',i,sep='')  }
  for(i in list(c(3,39), c(4,36), c(4, 39), c(8, 38), c(11, 37), c(13, 39))) {
    final <- cbind(final, 1*(final$focode== i[1] & final$offer_week == i[2]))
    names(final)[length(names(final))] <- paste('focode',i[1],'_offerweek',i[2],sep='')  }

  # model formula
  mod_form <- frac_paid_sept30 ~ paprikaexp + married + male + default + risky +
    hungry + incomesd + late + nopreviousloan + educdum2 + educdum3 + educdum4 +
    educdum5 + educdum6 + educdum7 + educdum8 + educdum9 + educdum10 + agedum2 +
    agedum3 + agedum4 + agedum5 + agedum6 + agedum7 + agedum8 + agedum9 + agedum10 +
    agedum11 + agedum12 + agedum13 + agedum15 + focode2 + focode3 + focode6 +
    focode7 + focode8 + focode9 + focode10 + focode13 + offerweek35 + offerweek36 +
    offerweek37 + offerweek38 + offerweek39 + focode3_offerweek39 +
    focode4_offerweek36 + focode4_offerweek39 + focode8_offerweek38 +
    focode11_offerweek37 + focode13_offerweek39

  # stage 1

  # Since NA variables will be dropped anyways, we can ditch them early
  final2 <- na.omit(final[,names(final) %in% c("fingerprint", all.vars(mod_form))])


  fakefinal <- final2[final2$fingerprint == 0,]
  fakefinal$fingerprint <- c(rep(0,282), rep(1,281))

  reps <- 1000
  gine_sim_res <- matrix(nrow = reps, ncol = 2)
  for (i in 1:reps) {
    fakefinal$fingerprint <- sample(fakefinal$fingerprint)

    mod1 <- lm(mod_form, data=fakefinal[fakefinal$fingerprint == 0, ])
    mod1step <- stepAIC(mod1, trace=FALSE)

    fakefinal$pred <- predict(mod1step, fakefinal)

    mod2 <- lm(frac_paid_sept30 - pred ~ pred,
               data=fakefinal[fakefinal$fingerprint == 1,])

    mod2cor <- pbph(mod1step, fakefinal$fingerprint, fakefinal)

    gine_sim_res[i,] <- c(summary(mod2)$coef[2,4], summary(mod2cor)$coef[2,4])
  }

  save(gine_sim_res, file = cachefile)
}

load(cachefile)
rm(cachefile)

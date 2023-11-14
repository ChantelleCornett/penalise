library(mstate)
library(tidyverse)
# test commit3
data("ebmt4")
head(ebmt4)

# makes transition state matrix
tmat <- mstate::transMat(x = list(c(2, 3, 5, 6), 
                                  c(4, 5, 6), 
                                  c(4, 5, 6), 
                                  c(5, 6),
                                  c(),
                                  c()),
                         names = c("Tx", "Rec", "AE", "Rec+AE", 
                                   "Rel", "Death"))
print(tmat)

# prepares dataset for MS modelling
msebmt <- msprep(data = ebmt4, trans = tmat, 
                 time = c(NA, "rec", "ae","recae", "rel", "srv"), 
                 status = c(NA, "rec.s", "ae.s", "recae.s", "rel.s", "srv.s"), 
                 keep = c("match", "proph", "year", "agecl"))

# subset based on patient 2 to check the structure of dataset
msebmt[msebmt$id == 2, ]

# make a cox proportional hazards model with multiple outcomes
coxph(Surv(Tstart, Tstop, status) ~ strata(trans) + match*factor(trans) + 
        proph*factor(trans) + year*factor(trans) + agecl*factor(trans),
      data = msebmt,
      method = "breslow")

# expand covariates in competing risks in stacked format for those in keep statement above
test_dt <- expand.covs(msebmt, covs = c("match", 
                                        "proph",
                                        "year",
                                        "agecl"),
                       append = TRUE,
                       longnames = TRUE)

# cox proportional hazards model based only on variables we stacked
coxph(Surv(Tstart, Tstop, status) ~ strata(trans) + ., 
      data = test_dt %>%
        select(-id, -from, -to,
               -time,
               -match,
               -proph,
               -year,
               -agecl),
      method = "breslow")


# Looks at how shrinkage and penalisation can be applied
# finds maximum number of transitions from transition matrix
n_trans <- max(tmat, na.rm = TRUE)

fits_wei <- vector(mode = "list", length = n_trans)

# fits models subsetting data on the number of transitions 
for (i in 1:n_trans){
  fits_wei[[i]] <- coxph(Surv(Tstart, Tstop, status) ~ match + 
                           proph + year + agecl,
                         data = subset(msebmt, trans == i),
                         method = "breslow")
}

#transition model specific calibration slope:
test_lp <- predict(fits_wei[[2]], type = 'lp')
coxph(Surv(Tstart, Tstop, status) ~ test_lp,
                       data = subset(msebmt, trans == 2),
                       method = "breslow")

##############################################################################
##               Comparative methods study on Penalisation                  ##
##                    Methods for Multi-State Models                        ##  
##                              Penalisation                                ##
##                       Author: Chantelle Cornett                          ##
##                           Date: 19NOV2023                                ##
##############################################################################

##############################################################################
## Change Log:                                                              ##
## Person                |   Date    | Changes Made                         ##
## _________________________________________________________________________##
## Chantelle Cornett     | 19NOV2023 | File initialisation                  ##
##############################################################################

library(mstate)
library(shrink)

msdat$entry <- msdat$Tstart
msdat$exit <- msdat$Tstop
msdat$event <- msdat$status

###################################
##.     MLE NO SHRINKAGE         ##
###################################

# make a cox proportional hazards model with multiple outcomes
noShrinkSimp <- coxph(Surv(entry, exit, event) ~ strata(trans) + genders*factor(trans) + ages*factor(trans),
      data = msdat,
      method = "breslow")

###################################
##     UNIFORM SHRINKAGE         ##
###################################

uniShrinkSimp <- shrink.coxph(fit = noShrinkSimp, type = "global", method = "jackknife", join = NULL, postfit = FALSE)

###################################
##  LASSO PENALISED LIKELIHOOD   ##
###################################
fit <- coxph(Surv(time, status) ~ lasso(covariate1 + covariate2, penalty = 0.1), data = data)
lassoSimp <- penMSM(type = "lasso", msdat, X, PSM1, PSM2, lambda1, lambda2, w, betastart, nu = 0.5, tol = 1e-10,
       max.iter = 50, trace = TRUE, diagnostics = TRUE, family = "coxph", poissonresponse = NULL,
       poissonoffset = NULL, constant.approx = 1e-8)
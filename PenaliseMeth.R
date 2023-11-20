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
###################################
##.     MLE NO SHRINKAGE         ##
###################################

# make a cox proportional hazards model with multiple outcomes
noShrinkSimp <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans) + genders*factor(trans) + ages*factor(trans),
      data = msdat,
      method = "breslow")

###################################
##     UNIFORM SHRINKAGE         ##
###################################

uniShrinkSimp <- shrink.coxph(fit = noShrinkSimp, type = "global", method = "jackknife", join = NULL, postfit = FALSE)

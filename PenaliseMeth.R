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
library(survival)
library(glmnet)
library(rrpack)
library(brms)
source("redrank.R")
msdat <- temp.data.cohort


msdat$entry <- msdat$Tstart
msdat$exit <- msdat$Tstop
msdat$event <- msdat$status

###################################
##.     MLE NO SHRINKAGE         ##
###################################

# make a cox proportional hazards model with multiple outcomes
noShrinkSimp <- coxph(Surv(entry, exit, event) ~ strata(trans) + gender*factor(trans) + age*factor(trans)+ BMI*factor(trans),
      data = msdat,
      method = "breslow")

###################################
##     UNIFORM SHRINKAGE         ##
###################################

# using parallel processing
cl <- makeCluster(5)
registerDoParallel(5)
start_time_fit_parallel <- Sys.time()

uniSimp_list<-(foreach(input=1:5, .combine=list, .multicombine=TRUE, 
                             .packages=c("shrink", "mstate", "survival"))) %dopar%{
uniShrinkSimp <- shrink.coxph(fit = noShrinkSimp, type = "global", method = "jackknife", join = NULL, postfit = FALSE)
                             }

diff_fit_parallel <- start_time_fit_parallel - end_time_fit_parallel
stopCluster(cl)
print("FINISH DATA GEN")
diff_fit_parallel
###################################
##  LASSO PENALISED LIKELIHOOD   ##
###################################

# PRE-PROCESSING DATA
d <- msdat[, c("entry", "exit","trans","event")]
X <- msdat[, c("trans","gender","age","BMI")]

# 10-FOLD CROSS VALIDATION FOR OPTIMAL VALUE OF LAMBDA
xmat <- as.matrix(X)
ysurv <- msdat[, c("Tstop","status")]
cvfit <- cv.glmnet(xmat, ysurv, family = "cox", type.measure = "C")
lassoSimp <- penMSM(type = "lasso", d, X, PSM1, PSM2, lambda1, lambda2, w, betastart, nu = 0.5, tol = 1e-10,
       max.iter = 50, trace = TRUE, diagnostics = TRUE, family = "coxph", constant.approx = 1e-8)

##########################################
##   FUSED LASSO PENALISED LIKELIHOOD   ##
##########################################

fusedSimp <- penMSM(type = "fused", d, X, PSM1, PSM2, lambda1, lambda2, w, betastart, nu = 0.5, tol = 1e-10,
                    max.iter = 50, trace = TRUE, diagnostics = TRUE, family = "coxph", poissonresponse = NULL,
                    poissonoffset = NULL, constant.approx = 1e-8)

###################################
##      REDUCED RANK METHOD      ##
###################################

# The reduced rank 2 solution
rr2 <- myredrank(Surv(Tstart, Tstop, status) ~ strata(trans) + gender*factor(trans) + age*factor(trans)+ BMI*factor(trans),
               data=msdat, R=2)
rr2$Alpha; rr2$Gamma; rr2$Beta; rr2$loglik

# The full rank solution
rr3 <- myredrank(Surv(entry, exit, event) ~ strata(trans) + genders*factor(trans) + ages*factor(trans),
               data=msdat, R=3)
rr3$Alpha; rr3$Gamma; rr3$Beta; rr3$loglik

###################################
##      BAYESIAN APPROACH        ##
###################################

fit <- brm(Surv(entry, exit, event) ~ strata(trans) + gender*factor(trans) + age*factor(trans)+BMI*factor(trans),
           data = msdat,
           family = brms::cox(),
           chains = 2,
           warmup = 1000,
           seed=123,
           prior = set_prior("horseshoe(0.5)", class = "b"))

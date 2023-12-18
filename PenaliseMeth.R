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

url <-
  "https://cran.r-project.org/src/contrib/Archive/apricom/apricom_1.0.0.tar.gz"
install.packages(url, type = "source", repos = NULL)
library(devtools)
devtools::load_all("/Users/m13477cc/Downloads/apricom/R")
library(apricom)
library(mstate)
library(shrink)
library(survival)
library(glmnet)
library(rrpack)
library(brms)
library(dplyr)
library(hdnom)
library(parallel)
library(penalized)
library(doParallel)
msdat <- temp.data.cohort


msdat$entry <- msdat$Tstart
msdat$exit <- msdat$Tstop
msdat$event <- msdat$status

###################################
##.     MLE NO SHRINKAGE         ##
###################################

# make a cox proportional hazards model with multiple outcomes
# assumping ph

noShrinkSimp <-
  coxph(
    Surv(entry, exit, event) ~  gender * factor(trans) + age *
      factor(trans) + BMI * factor(trans),
    data = msdat,
    method = "breslow",
    x = TRUE
  )

# not assuming PH - need to do!

###################################
##     UNIFORM SHRINKAGE         ##
###################################

# finds maximum number of transitions from transition matrix

tmat <- trans.illdeath()

n_trans <- max(tmat, na.rm = TRUE)

fits_wei <- vector(mode = "list", length = n_trans)

s <- 1 - (length(noShrinkSimp$coefficients) / 2267749)

# fits models subsetting data on the number of transitions
for (i in 1:n_trans) {
  fits_wei[[i]] <-
    coxph(
      Surv(entry, exit, event) ~ factor(gender) + factor(age) + BMI,
      method = "breslow",
      data = subset(msdat, trans == i)
    )
}

for (i in 1:n_trans) {
  for (j in 1:4) {
    fits_wei[[i]]$coefficients[j] <- fits_wei[[i]]$coefficients[j] * s
  }
}


###################################
##  LASSO PENALISED LIKELIHOOD   ##
###################################

X <- msdat[, c("trans", "gender", "age", "BMI")]
y <- cbind(time = msdat$Tstop, status = msdat$status)
cvfit <-
  cv.glmnet(
    x = X,
    y = y,
    family = "cox",
    type.measure = "C",
    data = msdat,
    parallel = TRUE,
    alpha = 1
  )

coef(cvfit, s = 'lambda.min')
final_model <-
  glmnet(
    x = X,
    y = y,
    family = "cox",
    data = msdat,
    lambda = 0.007408425,
    alpha = 1
  )
lassoSimp <- final_model


###################################
##      REDUCED RANK METHOD      ##
###################################

# The reduced rank 2 solution
rr2 <-
  redrank(
    Surv(Tstart, Tstop, status) ~  as.factor(gender) + as.factor(age) + BMI,
    data = na.exclude(msdat),
    R = 2
  )

rr2$Alpha
rr2$Gamma
rr2$Beta
rr2$loglik


###################################
##      BAYESIAN APPROACH        ##
###################################
model_covshrink <- bf(
  exit | cens(1 - status) ~ base + cand,
  base ~ 1,
  cand ~ strata(trans) + gender * factor(trans) + age *
    factor(trans),
  nl = TRUE
)
get_prior(model_covshrink, data = msdat)

prior_covshrink <-
  prior(normal(0, 2), coef = "Intercept", nlpar = "base") +
  prior(
    horseshoe(
      df = 3,
      scale_global = 2,
      df_global = 6,
      scale_slab = 1,
      df_slab = 6
    ),
    class = "b",
    nlpar = "cand"
  )

fit <-
  brm(
    model_covshrink,
    data = msdat,
    family = brms::cox(),
    chains = 2,
    warmup = 1000,
    seed = 123,
    prior = prior_covshrink
  )

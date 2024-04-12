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
install.packages(
  "/Users/m13477cc/R/Penalise/penalise/brms-2.19.0.tar",
  repos = NULL,
  type = "source",
  dependencies = TRUE
)
library(brms)
library(mstate)
library(shrink)
library(survival)
library(glmnet)
library(rrpack)
library(dplyr)
library(hdnom)
library(parallel)
library(penalized)
library(doParallel)
library(flexsurv)
library(ggplot2)
library(SurvMetrics)

###################################
##.     MLE NO SHRINKAGE         ##
###################################

# make a cox proportional hazards model with multiple outcomes
# assumping ph

patient_histories1$time2 <- 0
for (i in 1: length(patient_histories1$time)){
patient_histories1$time2[i] <- min(patient_histories1$time[i], patient_histories1$cens_time[i])
}


patient_histories1$status <- 0
for (i in 1: length(patient_histories1$time)){
  if (patient_histories1$cens_time[i] == patient_histories1$time2[i]){
    patient_histories1$status[i] <- 1
  }
}
noShrink_wei <- vector(mode = "list", length = n_trans)

noShrink_wei[[1]] <-
  coxph(
    Surv(time2,status) ~  gender + age + BMI ,
    data = patient_histories1,
    method = "breslow",
    x = TRUE
  )

noShrink_wei[[2]] <-
  coxph(
    Surv(time) ~  gender + age + BMI ,
    data = patient_histories2,
    method = "breslow",
    x = TRUE
  )

noShrink_wei[[3]] <-
  coxph(
    Surv(time) ~  gender + age + BMI ,
    data = patient_histories3,
    method = "breslow",
    x = TRUE
  )

###################################
##     UNIFORM SHRINKAGE         ##
###################################

# finds maximum number of transitions from transition matrix

tmat <- trans.illdeath()
n_trans <- max(tmat, na.rm = TRUE)

fits_wei <- vector(mode = "list", length = n_trans)

fits_wei[[1]] <-
  coxph(
    Surv(time) ~  gender + age + BMI ,
    data = patient_histories1,
    method = "breslow",
    x = TRUE
  )

fits_wei[[2]] <-
  coxph(
    Surv(time) ~  gender + age + BMI ,
    data = patient_histories2,
    method = "breslow",
    x = TRUE
  )

fits_wei[[3]] <-
  coxph(
    Surv(time) ~  gender + age + BMI ,
    data = patient_histories3,
    method = "breslow",
    x = TRUE
  )

# Calculating the heuristic shrinkage factor
s1 <- 1 - (length(fits_wei[[1]]$coefficients) / 2557)
s2 <- 1 - (length(fits_wei[[2]]$coefficients) / 3.78)
s3 <- 1 - (length(fits_wei[[3]]$coefficients) / 1.2)

# scales the coefficients according to the shrinkage factor

for (j in 1:3) {
  fits_wei[[1]]$coefficients[j] <- fits_wei[[1]]$coefficients[j] * s1
  fits_wei[[2]]$coefficients[j] <- fits_wei[[2]]$coefficients[j] * s2
  fits_wei[[3]]$coefficients[j] <- fits_wei[[3]]$coefficients[j] * s3
}

###################################
##  LASSO PENALISED LIKELIHOOD   ##
###################################

model_covshrinklas <- bf(time  ~ gender  + age + BMI, family = cox())

prior_covshrinklas <-
  prior(lasso(),
        class = "b")

fitlas <-
  brms::brm(
    model_covshrinklas,
    data = patient_histories1,
    family = brms::cox,
    chains = 2,
    warmup = 2000,
    iter = 2500,
    seed = 123,
    prior = prior_covshrinklas,
    cores = 2,
    threads = 2,
    control = list(max_treedepth = 30)
  )

fitlas2 <-
  brms::brm(
    model_covshrinklas,
    data = patient_histories2,
    family = brms::cox,
    chains = 2,
    warmup = 2000,
    iter = 2500,
    seed = 123,
    prior = prior_covshrinklas,
    cores = 2,
    threads = 2,
    control = list(max_treedepth = 30)
  )

fitlas3 <-
  brms::brm(
    model_covshrinklas,
    data = patient_histories3,
    family = brms::cox,
    chains = 1,
    warmup = 2000,
    iter = 2500,
    seed = 123,
    prior = prior_covshrinklas,
    cores = 2,
    threads = 2,
    control = list(max_treedepth = 30)
  )

###################################
##      BAYESIAN APPROACH        ##
###################################
model_covshrink <- bf(time ~ cand,
                      cand ~ gender  + age + BMI,
                      nl = TRUE)

prior_covshrink <-
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

fithist1 <-
  brm(
    model_covshrink,
    data = patient_histories1,
    family = brms::weibull,
    chains = 1,
    warmup = 1500,
    seed = 123,
    prior = prior_covshrink,
    control = list(max_treedepth = 30)
  )

model_covshrink2 <- bf(time ~ cand,
                       cand ~ gender  + age + BMI,
                       nl = TRUE)

prior_covshrink2 <-
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

fithist2 <-
  brm(
    model_covshrink2,
    data = patient_histories2,
    family = brms::weibull,
    chains = 1,
    warmup = 1500,
    seed = 123,
    prior = prior_covshrink2,
    control = list(max_treedepth = 30)
  )

model_covshrink3 <- bf(time ~ cand,
                       cand ~ gender  + age + BMI,
                       nl = TRUE)

prior_covshrink3 <-
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

fithist3 <-
  brm(
    model_covshrink3,
    data = patient_histories3,
    family = brms::weibull,
    chains = 1,
    warmup = 1500,
    seed = 123,
    prior = prior_covshrink3,
    control = list(max_treedepth = 30)
  )

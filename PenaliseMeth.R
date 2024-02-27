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
library(flexsurv)
library(ggplot2)
msdat <- state3

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
mssub1 <- subset(msdat, trans == 1)
mssub2 <- subset(msdat, trans == 2)
mssub3 <- subset(msdat, trans == 3)
noShrink_wei <- vector(mode = "list", length = n_trans)

  noShrink_wei[[1]]<-
    coxph(
      Surv(time) ~  gender + age + BMI ,
      data = msdat1,
      method = "breslow",
      x = TRUE
    )
    

  noShrink_wei[[2]] <-
    coxph(
      Surv(time) ~  gender + age + BMI ,
      data = msdat2,
      method = "breslow",
      x = TRUE
    )
  
  noShrink_wei[[3]] <-
    coxph(
      Surv(time) ~  gender + age + BMI ,
      data = msdat3,
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
fits_wei[[1]]<-
  coxph(
    Surv(time) ~  gender + age + BMI ,
    data = msdat1,
    method = "breslow",
    x = TRUE
  )

fits_wei[[2]] <-
  coxph(
    Surv(time) ~  gender + age + BMI ,
    data = msdat2,
    method = "breslow",
    x = TRUE
  )

fits_wei[[3]] <-
  coxph(
    Surv(time) ~  gender + age + BMI ,
    data = msdat3,
    method = "breslow",
    x = TRUE
  )

s1 <- 1 - (length(fits_wei[[1]]$coefficients) / 709.3)
s2 <- 1 - (length(fits_wei[[2]]$coefficients) / 709.3)
s3 <- 1 - (length(fits_wei[[3]]$coefficients) / 42.71)
# fits models subsetting data on the number of transitions

  for (j in 1:3) {
    fits_wei[[1]]$coefficients[j] <- fits_wei[[1]]$coefficients[j] * s1
  }

  for (j in 1:3) {
    fits_wei[[2]]$coefficients[j] <- fits_wei[[2]]$coefficients[j] * s2
  }
  for (j in 1:3) {
    fits_wei[[3]]$coefficients[j] <- fits_wei[[3]]$coefficients[j] * s3
  }

###################################
##  LASSO PENALISED LIKELIHOOD   ##
###################################

install.packages("/Users/m13477cc/R/Penalise/penalise/brms-2.19.0.tar", repos=NULL, type="source", dependencies = TRUE)
library(brms)
model_covshrinklas <- bf(
  time  ~ gender  + age + BMI, family = cox())


msdat1 <- subset(msdat, (from == 1 & to == 2))
msdat2 <- subset(msdat, (from == 1 & to == 3))
msdat3 <- subset(msdat, from == 2 & to == 3)
msdat1$exit <- msdat1$Tstop
get_prior(model_covshrinklas, data = msdat1)

prior_covshrinklas <-
  prior(
    lasso(
    ),
    class = "b"
  )

fitlas <-
  brms::brm(
    model_covshrinklas,
    data = msdat1,
    family = brms::cox,
    chains = 2,
    warmup = 2000,
    iter = 2500,
    seed = 123,
    prior = prior_covshrinklas,
    cores=2,
    threads=2,
    control = list(max_treedepth = 30) 
  )
fitlas2 <-
  brms::brm(
    model_covshrinklas,
    data = msdat2,
    family = brms::cox,
    chains = 2,
    warmup = 2000,
    iter = 2500,
    seed = 123,
    prior = prior_covshrinklas,
    cores=2,
    threads=2,
    control = list(max_treedepth = 30) 
  )
fitlas3 <-
  brms::brm(
    model_covshrinklas,
    data = msdat3,
    family = brms::cox,
    chains = 2,
    warmup = 2000,
    iter = 2500,
    seed = 123,
    prior = prior_covshrinklas,
    cores=2,
    threads=2,
    control = list(max_treedepth = 30) 
  )

library(tidybayes)
msdatlas <- msdat1 %>%
  add_residual_draws(fitlas) %>%
  median_qi()

ggplot(data=msdatlas, aes(sample = .residual)) +
  geom_qq() +
  geom_qq_line()+
  geom_abline(slope=1, intercept = 0, color = "red")

###################################
##      REDUCED RANK METHOD      ##
###################################

# The reduced rank 2 solution
attr(msdat, "trans") <- tmat

rr2 <-
  redrank(
    Surv(Tstart, Tstop, status) ~  as.factor(gender)+ as.factor(age) + BMI,
    data = msdat,
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
  time | cens(1 - status) ~ base + cand,
  base ~ 1,
  cand ~ gender  + age + BMI,
  nl = TRUE
)

msdat1 <- subset(msdat, trans == 1)
msdat1$exit <- msdat1$Tstop
get_prior(model_covshrink, data = msdat1)

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

fit1 <-
  brm(
    model_covshrink,
    data = msdat1,
    family = brms::weibull,
    chains = 2,
    warmup = 1500,
    seed = 123,
    prior = prior_covshrink,
    control = list(max_treedepth = 20) 
  )

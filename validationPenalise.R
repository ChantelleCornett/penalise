##############################################################################
##               Comparative methods study on Penalisation                  ##
##                    Methods for Multi-State Models                        ##
##                              Penalisation                                ##
##                       Author: Chantelle Cornett                          ##
##                           Date: 04JAN2024                                ##
##############################################################################

##############################################################################
## Change Log:                                                              ##
## Person                |   Date    | Changes Made                         ##
## _________________________________________________________________________##
## Chantelle Cornett     | 04JAN2024 | File initialisation                  ##
##############################################################################

library(predRupdate)

# First need a validation dataset
library("flexsurv")
library("hesim")
library("data.table")

# Simulate multi-state dataset -------------------------------------------------
sim_onc3_data <- function(n = 100000, seed = NULL){
  
  if (!is.null(seed)) set.seed(seed)
  
  # Data 
  age_mu <- 60
  data <- data.table(
    intercept = 1,
    strategy_id = 1,
    strategy_name = sample(c("SOC", "New 1", "New 2"), n, replace = TRUE,
                           prob = c(1/3, 1/3, 1/3)),
    patient_id = 1:n,
    gender = rbinom(n, 1, .5),
    age = sample(c(1,2,3), n, replace=TRUE),
    BMI = rnorm(n, mean = 24, sd = 3)
  )
  data[, `:=` (new1 = ifelse(strategy_name == "New 1", 1, 0),
               new2 = ifelse(strategy_name == "New 2", 1, 0))]
  attr(data, "id_vars") <- c("strategy_id", "patient_id")
  
  # Transition matrix
  tmat <- rbind(
    c(NA, 1,  2),
    c(NA,  NA, 3),
    c(NA, NA, NA)
  )
  trans_dt <- create_trans_dt(tmat)
  
  # Parameters for each transition
  get_scale <- function(shape, mean) {
    scale <- mean/(gamma(1 + 1/shape))
    scale_ph <- scale^{-shape}
    return(scale_ph)
  }
  
  matrixv <- function(v) {
    x <- matrix(v); colnames(x) <- "intercept"
    return(x)
  }
  
  params_wei <- function(shape, mean,
                         beta_new1 = log(1), 
                         beta_new2 = log(1),
                         beta_age, beta_female, beta_bmi){
    log_shape <- matrixv(log(shape))
    scale = get_scale(shape, mean)
    beta_intercept <- log(scale) - mean(data$age) * beta_age - mean(data$gender)* beta_female - mean(data$BMI)* beta_bmi
    scale_coefs <-  matrix(c(beta_intercept, beta_new1, beta_new2, 
                             beta_age, beta_female, beta_bmi), 
                           ncol = 6)
    colnames(scale_coefs) <- c("intercept", "new1", "new2", "age", "gender", "BMI")
    params_surv(coefs = list(shape = log_shape,
                             scale = scale_coefs),
                dist = "weibullPH")
  }
  
  mstate_params <- params_surv_list(
    
    # 1. S -> P
    params_wei(shape = 2, mean = 50, 
               beta_new1 = log(.7), beta_new2 = log(.6),
               beta_female = log(1.4), beta_age = log(1.03), beta_bmi = log(1.5)),
    
    # 2. S -> D
    params_wei(shape = 25, mean = 1000,
               beta_new1 = log(.85), beta_new2 = log(.75),
               beta_female = log(1.2), beta_age = log(1.02), beta_bmi = log(1.5)),
    
    # 3. P -> D
    params_wei(shape = 3.5, mean = 50, beta_new1 = log(1),
               beta_female = log(1.3), beta_age = log(1.02),  beta_bmi = log(1.5))
  )
  
  # Create multi-state model
  mstatemod <- create_IndivCtstmTrans(mstate_params, 
                                      input_data = data,
                                      trans_mat = tmat,
                                      clock = "reset",
                                      start_age = data$age)
  
  # Simulate data
  ## Observed "latent" transitions
  sim <- mstatemod$sim_disease(max_age = 100)
  sim[, c("sample", "grp_id", "strategy_id") := NULL]
  sim <- cbind(
    data[match(sim$patient_id, data$patient_id)][, patient_id := NULL],
    sim
  )
  sim[, ":=" (intercept = NULL, strategy_id = NULL, status = 1, added = 0)]
  
  ## Add all possible states for each transition
  ### Observed 1->2 add 1->3
  sim_13 <- sim[from == 1 & to == 2]
  sim_13[, ":=" (to = 3, status = 0, final = 0,  added = 1)]
  sim <- rbind(sim, sim_13)
  
  ### Observed 1->3 add 1->2
  sim_12 <- sim[from == 1 & to == 3 & added == 0]
  sim_12[, ":=" (to = 2, status = 0, final = 0, added = 1)]
  sim <- rbind(sim, sim_12)
  
  ### Sort and clean
  sim <- merge(sim, trans_dt, by = c("from", "to")) # Add transition ID
  setorderv(sim, c("patient_id", "from", "to"))
  sim[, added := NULL]
  
  ## Add right censoring
  rc <- data.table(patient_id = 1:n,
                   time = stats::rexp(n, rate = 1/15))
  sim[, time_rc := rc[match(sim$patient_id, rc$patient_id)]$time]
  sim[, status := ifelse(time_stop < 15 & time_stop < time_rc, status, 0)]
  sim[, time_stop := pmin(time_stop, 15, time_rc)]
  sim <- sim[time_start <= pmin(time_rc, 15)]
  
  ## Final data cleaning
  sim[, strategy_id := fcase(
    strategy_name == "SOC", 1L,
    strategy_name == "New 1", 2L,
    strategy_name == "New 2", 3L
  )]
  sim[, strategy_name := factor(strategy_id, 
                                levels = c(1, 2, 3),
                                labels = c("SOC", "New 1", "New 2"))]
  label_states <- function (x) {
    fcase(
      x == 1, "1",
      x == 2, "2",
      x == 3, "3"
    )
  }
  sim[, from := label_states(from)]
  sim[, to := label_states(to)]
  sim[, c("new1", "new2", "final", "time_rc") := NULL]
  
  # Return
  sim[, time := time_stop - time_start]
  return(sim[, ])
}
onc3 <- sim_onc3_data(n = 100000, seed = 220)

# Check that coefficient estimates are consistent with "truth"
fit_weibull <- function(i) {
  flexsurvreg(Surv(time, status) ~ gender + age + BMI,
              data = onc3, subset = (transition_id == i), dist = "weibullPH")
}
fit_weibull(1)
fit_weibull(2)
fit_weibull(3)

# Panel data version -----------------------------------------------------------
onc3p <- copy(onc3)
onc3p[, n := 1:.N, by = c("patient_id", "time_start")]
onc3p[, c("transition_id", "time") := NULL]

# Time 0
onc3p_t0 <- onc3p[time_start == 0 & n == 1]
onc3p_t0[, c("time_stop", "n", "to", "status") := NULL]
setnames(onc3p_t0, c("time_start", "from"), c("time", "state"))

# Time > 0
onc3p[, mstatus := mean(status), by = c("patient_id", "time_start")]
onc3p <- onc3p[status == 1 | (mstatus == 0 & n == 1)]
onc3p[, state := ifelse(mstatus == 0, from, to)]
onc3p[, c("time_start", "n", "from",
          "to", "status", "mstatus") := NULL]
setnames(onc3p, "time_stop", "time")

# Full panel
onc3p <- rbind(onc3p_t0, onc3p)
setorderv(onc3p, c("patient_id", "time"))
onc3p[, state_id := factor(
  state, 
  levels = c("1", "2", "3"),
  labels = 1:3)]

onc3 <- onc3[,-c(3,11)]
colnames(onc3) <- c("from","to","gender","age","BMI","id","Tstart","Tstop","status","trans","time")

nrow(subset(onc3, from ==1 & to == 3))
nrow(subset(onc3, from ==1 & to ==2))
nrow(subset(onc3, from ==2 & to ==3))

state3valid <- onc3

state3valid$age <- as.factor(state3valid$age)
state3valid$gender <- as.factor(state3valid$gender)
state3valid <-dummy_vars(state3valid)
state3valid <- cbind(state3valid, onc3$age)
colnames(state3valid) <- c("from","to","BMI","id","Tstart","Tstop","status","trans","time","gender1","gender","lessTwenty","twentyToForty","greaterForty", "age")

#### UNIFORM SHRINKAGE ####

library(survival)
validWei1 <- subset(state3valid, (from == 1& to == 2))
library(flexsurv)

########### NO PEN 2 #################
validNoPen2 <- coxph(Surv(time) ~ offset(0.01772*age + 0.10141 * BMI -0.79066*gender), data=patient2_histories2)

# Calibration plots
library(survival)
library(predtools)
library(probably)
library(timeROC)
validWei1$prednp1 <- predict(validNoPen1, data = validWei1, type = "survival")
real <- coxph(Surv(time,status) ~ 1, data = validWei1)
validWei1$realprobs <- predict(real, data= validWei1)
validWei1$predtime <- predict(real, data= validWei1, type="survival")
cal_plot_regression(validWei1, truth = "realprobs", estimate = "prednp1", smooth = TRUE)
rocnp1 <- timeROC(T=validWei1$time[1:200],cause = 1, delta = validWei1$status[1:200], marker=validWei1$realprobs[1:200], times = seq(0,15, by =0.5),iid=TRUE)
plot(rocnp1, FP = 2, add = FALSE,conf.int = FALSE,conf.band = FALSE, col = "black")
plot(rocnp1,time=15,add = FALSE, col = "black",title=FALSE, lwd=2)
########### NO PEN 1 #################
validNoPen1 <- coxph(Surv(time) ~ offset(0.0047995*age + 0.0064177 * BMI -0.1096818*gender), data=patient2_histories1)
validWei2$prednp2 <- predict(validNoPen2, data = validWei2, type = "survival")
real2 <- coxph(Surv(time) ~ 1, data = validWei2)
validWei2$realprobs <- predict(real2, data= validWei2, type="survival")
cal_plot_regression(validWei2, truth = "realprobs", estimate = "prednp2", smooth = TRUE)
########### NO PEN 3 #################
validNoPen3 <- coxph(Surv(time) ~ offset(-0.00008014*age -0.0002384* BMI + 0.002305*gender), data=patient2_histories3)
validWei3$prednp3 <- predict(validNoPen3, data = validWei3, type = "survival")
real3 <- coxph(Surv(time) ~ 1, data = validWei3)
validWei3$realprobs <- predict(real3, data= validWei3, type="survival")
cal_plot_regression(validWei3, truth = "realprobs", estimate = "prednp3", smooth = TRUE)
############# G SHRINK 1 ###############
validGShrink1 <- coxph(Surv(time) ~ offset(0.0047939*age + 0.0064102 * BMI -0.1095531*gender), data=patient2_histories1)
validWei1$predgs1 <- predict(validGShrink1, data = validWei1, type = "survival")
realg1<- coxph(Surv(time) ~ 1, data = validWei1)
validWei1$gprobs <- predict(realg1, data= validWei1, type="survival")
cal_plot_regression(validWei1, truth = "gprobs", estimate = "predgs1", smooth = TRUE)
############# G SHRINK 2 ###############
validGShrink2 <- coxph(Surv(time) ~ offset(0.003655*age + 0.020925 * BMI -0.163152*gender), data=patient2_histories2)
validWei2$predgs2 <- predict(validGShrink2, data = validWei2, type = "survival")
realg2<- coxph(Surv(time) ~ 1, data = validWei2)
validWei2$gprobs <- predict(realg2, data= validWei2, type="survival")
cal_plot_regression(validWei2, truth = "gprobs", estimate = "predgs2", smooth = TRUE)
############# G SHRINK 3 ###############
validGShrink3 <- coxph(Surv(time) ~ offset(-0.006149*age - 0.018992 * BMI -0.06170*gender), data=validWei3)
validWei3$predgs3 <- predict(validGShrink3, data = validWei3, type = "survival")
realg3<- coxph(Surv(time) ~ 1, data = validWei3)
validWei3$gprobs <- predict(realg3, data= validWei3, type="survival")
cal_plot_regression(validWei3, truth = "gprobs", estimate = "predgs3", smooth = TRUE)
############# LASSO 1 ###############
validLAS1 <- coxph(Surv(time) ~ offset(0.005216672*age + 0.018911185 * BMI + 0.019396405*gender), data=validWei1)

############# LASSO 2 ###############
validLAS2 <- coxph(Surv(time) ~ offset(0.005216672*age + 0.018911185 * BMI + 0.019306405*gender), data=validWei2)

############# LASSO 3 ###############
validLAS3 <- coxph(Surv(time) ~ offset(-0.004146397*age -0.021360847* BMI -0.014194088*gender), data=validWei3)


validWei2$from <- as.numeric(validWei2$from)
validWei2$to <- as.numeric(validWei2$to)
cumHaz<- basehaz(fits_wei[[2]], centered = FALSE)
cumHaz <- cumHaz[,c(2,1)]

fitsWei2_mi <- as.data.frame(t(fits_wei[[2]][["coefficients"]]))
library(predRupdate)
x <- pred_input_info(
  model_type = "survival",
  model_info = fitsWei2_mi,
  cum_hazard = cumHaz
) 

pred_validate(
  x=x,
  new_data = validWei2,
  binary_outcome =NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = 15,
  cal_plot = TRUE
)

# validation for transition from ill to dead
validWei3 <- subset(state3valid, from == 2 & to==3)

validWei3$from <- as.numeric(validWei3$from)
validWei3$to <- as.numeric(validWei3$to)
cumHaz<- basehaz(validNoPen1, centered = FALSE)
cumHaz <- cumHaz[,c(2,1)]

validNoPen1_mi <- as.data.frame(t(c(0.313246, 0.031646, 0.401673)))
colnames(validNoPen1_mi) <- c("gender", "age","BMI")
library(predRupdate)
x <- pred_input_info(
  model_type = "survival",
  model_info = validNoPen1_mi,
  cum_hazard = cumHaz
) 

pred_validate(
  x=x,
  new_data = validWei1,
  binary_outcome =NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = 15,
  cal_plot = TRUE
)

cumHaz<- basehaz(validNoPen3, centered = FALSE)
cumHaz <- cumHaz[,c(2,1)]

validNoPen3_mi <- as.data.frame(t(c(0.71741, 0.02449, 0.46039)))
colnames(validNoPen3_mi) <- c("gender", "age","BMI")

library(predRupdate)
x <- pred_input_info(
  model_type = "survival",
  model_info = validNoPen3_mi,
  cum_hazard = cumHaz
) 

pred_validate(
  x=x,
  new_data = validWei3,
  binary_outcome =NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = cumHaz$time[1029],
  cal_plot = TRUE
)

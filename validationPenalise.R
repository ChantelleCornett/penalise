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

# Need to put into the correct form for model_info
validWei1 <- subset(state3valid, trans == 1)

validWei1$from <- as.numeric(validWei1$from)
validWei1$to <- as.numeric(validWei1$to)
fit_cox <- coxph(Surv(time, status) ~ 1 ,data = msdat, subset = trans==1, method = "breslow")
cumHaz<- basehaz(fit_cox)
cumHaz <- cumHaz[,c(2,1)]

fitsWei1_mi <- as.data.frame(t(fits_wei[[1]][["coefficients"]]))

colnames(fitsWei1_mi) <- c("gender","twentyToForty","greaterForty","BMI")

# validation for transition from healthy to ill
x <- pred_input_info(
  model_type = "survival",
  model_info = fitsWei1_mi,
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

# validation for transition from healthy to dead
validWei2 <- subset(state3valid, trans == 2)

validWei2$from <- as.numeric(validWei2$from)
validWei2$to <- as.numeric(validWei2$to)
fit_cox <- coxph(Surv(time, status) ~ 1 ,data = msdat, subset = trans==2, method = "breslow")
cumHaz<- basehaz(fit_cox)
cumHaz <- cumHaz[,c(2,1)]

fitsWei2_mi <- as.data.frame(t(fits_wei[[2]][["coefficients"]]))

colnames(fitsWei2_mi) <- c("gender","twentyToForty","greaterForty","BMI")

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
validWei3 <- subset(state3valid, trans == 3)

validWei3$from <- as.numeric(validWei3$from)
validWei3$to <- as.numeric(validWei3$to)
fit_cox <- coxph(Surv(time, status) ~ 1 ,data = state3valid1, method = "breslow")
cumHaz<- basehaz(fit_cox)
cumHaz <- cumHaz[,c(2,1)]

fitsWei3_mi <- as.data.frame(t(fits_wei[[3]][["coefficients"]]))

colnames(fitsWei3_mi) <- c("gender","twentyToForty","greaterForty","BMI")

x <- pred_input_info(
  model_type = "survival",
  model_info = fitsWei3_mi,
  cum_hazard = cumHaz
) 

pred_validate(
  x=x,
  new_data = validWei3,
  binary_outcome =NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = 15,
  cal_plot = TRUE
)

#### LASSO PENALISED LIKELIHOOD ####

state3valid$from <- as.numeric(state3valid$from)
state3valid$to <- as.numeric(state3valid$to)
fit_cox <- coxph(Surv(time, status) ~ 1 ,data = validWei1, method = "breslow")
cumHaz<- basehaz(fit_cox)
cumHaz <- cumHaz[,c(2,1)]
##########################################
lassoList <- c()
for (i in 1:length(lassoSimp$beta)){
  lassoList <- c(lassoList, lassoSimp$beta[i])
}
lasso_mi <- as.data.frame(t(lassoList))

colnames(lasso_mi) <- c("gender","age","BMI")

x <- pred_input_info(
  model_type = "survival",
  model_info = lasso_mi,
  cum_hazard = cumHaz
) 

pred_validate(
  x=x,
  new_data = state3valid,
  binary_outcome =NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = 15,
  cal_plot = TRUE
)

#### MLE NO SHRINKAGE ####

noShrinkSimp_mi <- as.data.frame(t(noShrinkSimp[["coefficients"]]))

colnames(noShrinkSimp_mi) <- c("gender","twentyToForty","greaterForty","age","BMI","gender:twentyToForty","gender:greaterForty","age:twentyToForty","age:greaterForty", "BMI:twentyToForty","BMI:greaterForty")

x <- pred_input_info(
  model_type = "survival",
  model_info = noShrinkSimp_mi,
  cum_hazard = NULL
) 

pred_validate(
  x=x,
  new_data = state3valid,
  binary_outcome = NULL,
  survival_time = time,
  event_indicator = status,
  time_horizon = NULL,
  cal_plot = TRUE,
  ...
)

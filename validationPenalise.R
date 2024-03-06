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

library("predRupdate")
library("flexsurv")
library("hesim")
library("data.table")
library("survival")
library("predtools")
library("probably")
library("timeROC")

######################################
############ UNPENALIZED #############
######################################

########### NO PEN 1 #################

validNoPen1 <-
  coxph(Surv(time) ~ offset(0.0047995 * age + 0.0064177 * BMI - 0.1096818 *
                              gender),
        data = patient2_histories1)

cumHaz <- basehaz(validNoPen1, centered = FALSE)
cumHaz <- cumHaz[, c(2, 1)]

validNoPen1_mi <- as.data.frame(t(c(-0.1096818, 0.0047995, 0.0064177)))
colnames(validNoPen1_mi) <- c("gender", "age", "BMI")

x <- pred_input_info(model_type = "survival",
                     model_info = validNoPen1_mi,
                     cum_hazard = cumHaz)

pred_validate(
  x = x,
  new_data = patient2_histories1,
  binary_outcome = NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = 15,
  cal_plot = TRUE
)

# ROC plots

patient2_histories1$prednp1 <-
  predict(validNoPen1, data = patient2_histories1, type = "survival")

real1 <- coxph(Surv(time) ~ 1, data = patients2_histories1)

patient2_histories1$realprobs <-
  predict(real1, data = patient2_histories1, type = "survival")


rocnp1 <-
  timeROC(
    T = patient2_histories1$time[1:200],
    cause = 1,
    delta = patient2_histories1$status[1:200],
    marker = patient2_histories1$realprobs[1:200],
    times = seq(0, 15, by = 0.5),
    iid = TRUE
  )

########### NO PEN 2 #################

validNoPen2 <-
  coxph(Surv(time) ~ offset(0.01772 * age + 0.10141 * BMI - 0.79066 * gender),
        data = patient2_histories2)

cumHaz <- basehaz(validNoPen2, centered = FALSE)
cumHaz <- cumHaz[, c(2, 1)]

validNoPen2_mi <- as.data.frame(t(c(-0.79066, 0.01772, 0.10141)))
colnames(validNoPen2_mi) <- c("gender", "age", "BMI")

x <- pred_input_info(model_type = "survival",
                     model_info = validNoPen2_mi,
                     cum_hazard = cumHaz)

pred_validate(
  x = x,
  new_data = patient2_histories2,
  binary_outcome = NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = 15,
  cal_plot = TRUE
)

# ROC plots

patient2_histories2$prednp1 <-
  predict(validNoPen2, data = patient2_histories2, type = "survival")

real2 <- coxph(Surv(time) ~ 1, data = patients2_histories2)

patient2_histories2$realprobs <-
  predict(real2, data = patient2_histories2, type = "survival")


rocnp2 <-
  timeROC(
    T = patient2_histories2$time[1:200],
    cause = 1,
    delta = patient2_histories2$status[1:200],
    marker = patient2_histories2$realprobs[1:200],
    times = seq(0, 15, by = 0.5),
    iid = TRUE
  )

########### NO PEN 3 #################

validNoPen3 <-
  coxph(Surv(time) ~ offset(-0.00008014 * age - 0.0002384 * BMI + 0.002305 *
                              gender),
        data = patient2_histories3)

cumHaz <- basehaz(validNoPen3, centered = FALSE)
cumHaz <- cumHaz[, c(2, 1)]

validNoPen3_mi <- as.data.frame(t(c(0.002305, -0.00008014, -0.0002384)))
colnames(validNoPen3_mi) <- c("gender", "age", "BMI")

x <- pred_input_info(model_type = "survival",
                     model_info = validNoPen3_mi,
                     cum_hazard = cumHaz)

pred_validate(
  x = x,
  new_data = patient2_histories3,
  binary_outcome = NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = 15,
  cal_plot = TRUE
)

# ROC plots

patient2_histories3$prednp3 <-
  predict(validNoPen3, data = patient2_histories3, type = "survival")

real3 <- coxph(Surv(time) ~ 1, data = patients2_histories3)

patient2_histories3$realprobs <-
  predict(real3, data = patient2_histories3, type = "survival")

rocnp3 <-
  timeROC(
    T = patient2_histories3$time[1:200],
    cause = 1,
    delta = patient2_histories3$status[1:200],
    marker = patient2_histories3$realprobs[1:200],
    times = seq(0, 15, by = 0.5),
    iid = TRUE
  )

########################################
######### GLOBAL SHRINKAGE #############
########################################

############# G SHRINK 1 ###############

validGShrink1 <-
  coxph(Surv(time) ~ offset(0.0047939 * age + 0.0064102 * BMI - 0.1095531 *
                              gender),
        data = patient2_histories1)

cumHaz <- basehaz(validGShrink1, centered = FALSE)
cumHaz <- cumHaz[, c(2, 1)]

validGShrink1_mi <- as.data.frame(t(c(-0.1095531, 0.0047939, 0.0064102)))
colnames(validGShrink1_mi) <- c("gender", "age", "BMI")

x <- pred_input_info(model_type = "survival",
                     model_info = validGShrink1_mi,
                     cum_hazard = cumHaz)

pred_validate(
  x = x,
  new_data = patient2_histories1,
  binary_outcome = NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = 15,
  cal_plot = TRUE
)

# ROC plots

patient2_histories1$predgs1 <-
  predict(validGShrink1, data = patient2_histories1, type = "survival")

real1 <- coxph(Surv(time) ~ 1, data = patients2_histories1)

patient2_histories1$realprobs <-
  predict(real1, data = patient2_histories1, type = "survival")

rocgs1 <-
  timeROC(
    T = patient2_histories1$time[1:200],
    cause = 1,
    delta = patient2_histories1$status[1:200],
    marker = patient2_histories1$realprobs[1:200],
    times = seq(0, 15, by = 0.5),
    iid = TRUE
  )


############# G SHRINK 2 ###############

validGShrink2 <-
  coxph(Surv(time) ~ offset(0.003655 * age + 0.020925 * BMI - 0.163152 *
                              gender),
        data = patient2_histories2)

cumHaz <- basehaz(validGShrink2, centered = FALSE)
cumHaz <- cumHaz[, c(2, 1)]

validGShrink2_mi <- as.data.frame(t(c(-0.163152, 0.003655, 0.020925)))
colnames(validGShrink2_mi) <- c("gender", "age", "BMI")

x <- pred_input_info(model_type = "survival",
                     model_info = validGShrink2_mi,
                     cum_hazard = cumHaz)

pred_validate(
  x = x,
  new_data = patient2_histories2,
  binary_outcome = NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = 15,
  cal_plot = TRUE
)

# ROC plots

patient2_histories2$predgs2 <-
  predict(validGShrink2, data = patient2_histories2, type = "survival")

real2 <- coxph(Surv(time) ~ 1, data = patients2_histories2)

patient2_histories2$realprobs <-
  predict(real2, data = patient2_histories2, type = "survival")

rocgs2 <-
  timeROC(
    T = patient2_histories2$time[1:200],
    cause = 1,
    delta = patient2_histories2$status[1:200],
    marker = patient2_histories2$realprobs[1:200],
    times = seq(0, 15, by = 0.5),
    iid = TRUE
  )


############# G SHRINK 3 ###############

validGShrink3 <-
  coxph(Surv(time) ~ offset(0.0001202 * age + 0.0003577 * BMI - 0.0034578 *
                              gender),
        data = patient2_histories3)

cumHaz <- basehaz(validGShrink3, centered = FALSE)
cumHaz <- cumHaz[, c(2, 1)]

validGShrink3_mi <- as.data.frame(t(c(-0.0034578, 0.0001202, 0.0003577)))
colnames(validGShrink3_mi) <- c("gender", "age", "BMI")

x <- pred_input_info(model_type = "survival",
                     model_info = validGShrink3_mi,
                     cum_hazard = cumHaz)

pred_validate(
  x = x,
  new_data = patient2_histories3,
  binary_outcome = NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = 15,
  cal_plot = TRUE
)

# ROC plots

patient2_histories3$predgs3 <-
  predict(validGShrink3, data = patient2_histories3, type = "survival")

real3 <- coxph(Surv(time) ~ 1, data = patients2_histories3)

patient2_histories3$realprobs <-
  predict(real3, data = patient2_histories3, type = "survival")

rocgs3 <-
  timeROC(
    T = patient2_histories3$time[1:200],
    cause = 1,
    delta = patient2_histories3$status[1:200],
    marker = patient2_histories3$realprobs[1:200],
    times = seq(0, 15, by = 0.5),
    iid = TRUE
  )

#####################################
########## LASSO PEN ################
#####################################

############# LASSO 1 ###############

validLShrink <-
  coxph(Surv(time) ~ offset(0.0047990 * age + 0.006402035 * BMI - 0.109169861 *
                              gender),
        data = patient2_histories1)

cumHaz <- basehaz(validLShrink, centered = FALSE)
cumHaz <- cumHaz[, c(2, 1)]

validLShrink1_mi <- as.data.frame(t(c(-0.109169861, 0.0047990, 0.006402035)))
colnames(validLShrink1_mi) <- c("gender", "age", "BMI")

x <- pred_input_info(model_type = "survival",
                     model_info = validLShrink1_mi,
                     cum_hazard = cumHaz)

pred_validate(
  x = x,
  new_data = patient2_histories1,
  binary_outcome = NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = 15,
  cal_plot = TRUE
)

# ROC plots

patient2_histories1$predls1 <-
  predict(validLShrink1, data = patient2_histories1, type = "survival")

real1 <- coxph(Surv(time) ~ 1, data = patients2_histories1)

patient2_histories1$realprobs <-
  predict(real1, data = patient2_histories1, type = "survival")

rocls1 <-
  timeROC(
    T = patient2_histories1$time[1:200],
    cause = 1,
    delta = patient2_histories1$status[1:200],
    marker = patient2_histories1$realprobs[1:200],
    times = seq(0, 15, by = 0.5),
    iid = TRUE
  )

############# LASSO 2 ###############

validLShrink2 <-
  coxph(Surv(time) ~ offset(0.0072141454 * age + 0.01690258 * BMI - 0.0172223 *
                              gender),
        data = patient2_histories2)


cumHaz <- basehaz(validLShrink2, centered = FALSE)
cumHaz <- cumHaz[, c(2, 1)]

validLShrink2_mi <- as.data.frame(t(c(-0.0172223, 0.0072141454, 0.01690258)))
colnames(validLShrink2_mi) <- c("gender", "age", "BMI")

x <- pred_input_info(model_type = "survival",
                     model_info = validLShrink2_mi,
                     cum_hazard = cumHaz)

pred_validate(
  x = x,
  new_data = patient2_histories2,
  binary_outcome = NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = 15,
  cal_plot = TRUE
)

# ROC plots

patient2_histories2$predls2 <-
  predict(validLShrink2, data = patient2_histories2, type = "survival")

real2 <- coxph(Surv(time) ~ 1, data = patients2_histories2)

patient2_histories2$realprobs <-
  predict(real2, data = patient2_histories2, type = "survival")

rocls2 <-
  timeROC(
    T = patient2_histories2$time[1:200],
    cause = 1,
    delta = patient2_histories2$status[1:200],
    marker = patient2_histories2$realprobs[1:200],
    times = seq(0, 15, by = 0.5),
    iid = TRUE
  )

############# LASSO 3 ###############

validLAS3 <-
  coxph(
    Surv(time) ~ offset(-0.0002839332 * age - 0.0006460146 * BMI + 0.0005969886 *
                          gender),
    data = patient2_histories3
  )


cumHaz <- basehaz(validLShrink3, centered = FALSE)
cumHaz <- cumHaz[, c(2, 1)]

validLShrink3_mi <- as.data.frame(t(c(0.0005969886, -0.0002839332,-0.0006460146)))
colnames(validLShrink3_mi) <- c("gender", "age", "BMI")

x <- pred_input_info(model_type = "survival",
                     model_info = validLShrink3_mi,
                     cum_hazard = cumHaz)

pred_validate(
  x = x,
  new_data = patient2_histories3,
  binary_outcome = NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = 15,
  cal_plot = TRUE
)

# ROC plots

patient2_histories3$predls3 <-
  predict(validLShrink3, data = patient2_histories3, type = "survival")

real3 <- coxph(Surv(time) ~ 1, data = patients2_histories3)

patient2_histories3$realprobs <-
  predict(real3, data = patient2_histories3, type = "survival")

rocls3 <-
  timeROC(
    T = patient2_histories3$time[1:200],
    cause = 1,
    delta = patient2_histories3$status[1:200],
    marker = patient2_histories3$realprobs[1:200],
    times = seq(0, 15, by = 0.5),
    iid = TRUE
  )


###########################################
################# HORSESHOE ###############
###########################################

############# H 1 ###############

validHShrink <-
  coxph(
    Surv(time) ~ offset(-0.0000372326 * age + 0.0005008169 * BMI + 0.0038252140 *
                          gender),
    data = patient2_histories1
  )

cumHaz <- basehaz(validHShrink, centered = FALSE)
cumHaz <- cumHaz[, c(2, 1)]

validHShrink1_mi <- as.data.frame(t(c(0.0038252140, -0.0000372326, 0.0005008169)))
colnames(validHShrink1_mi) <- c("gender", "age", "BMI")

x <- pred_input_info(model_type = "survival",
                     model_info = validHShrink1_mi,
                     cum_hazard = cumHaz)

pred_validate(
  x = x,
  new_data = patient2_histories1,
  binary_outcome = NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = 15,
  cal_plot = TRUE
)

# ROC plots

patient2_histories1$predhs1 <-
  predict(validHShrink, data = patient2_histories1, type = "survival")

real1 <- coxph(Surv(time) ~ 1, data = patients2_histories1)

patient2_histories1$realprobs <-
  predict(real1, data = patient2_histories1, type = "survival")

rochs1 <-
  timeROC(
    T = patient2_histories1$time[1:200],
    cause = 1,
    delta = patient2_histories1$status[1:200],
    marker = patient2_histories1$realprobs[1:200],
    times = seq(0, 15, by = 0.5),
    iid = TRUE
  )

############# H 2 ###############
validHShrink2 <-
  coxph(Surv(time) ~ offset(0.0072141454 * age + 0.01690258 * BMI - 0.0172223 *
                              gender),
        data = patient2_histories2)

cumHaz <- basehaz(validHShrink2, centered = FALSE)
cumHaz <- cumHaz[, c(2, 1)]

validHShrink2_mi <- as.data.frame(t(c(-0.0172223, 0.0072141454, 0.01690258)))
colnames(validHShrink2_mi) <- c("gender", "age", "BMI")

x <- pred_input_info(model_type = "survival",
                     model_info = validHShrink2_mi,
                     cum_hazard = cumHaz)

pred_validate(
  x = x,
  new_data = patient2_histories2,
  binary_outcome = NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = 15,
  cal_plot = TRUE
)

# ROC plots

patient2_histories2$predhs2 <-
  predict(validHShrink2, data = patient2_histories2, type = "survival")

real2 <- coxph(Surv(time) ~ 1, data = patients2_histories2)

patient2_histories2$realprobs <-
  predict(real2, data = patient2_histories2, type = "survival")

rochs2 <-
  timeROC(
    T = patient2_histories2$time[1:200],
    cause = 1,
    delta = patient2_histories2$status[1:200],
    marker = patient2_histories2$realprobs[1:200],
    times = seq(0, 15, by = 0.5),
    iid = TRUE
  )


############# H 3 ###############

validHShrink3 <-
  coxph(Surv(time) ~ offset(0.0003452605 * age + 0.0010034381 * BMI - 0.0030383568 *
                              gender),
        data = patient2_histories3)


cumHaz <- basehaz(validHShrink3, centered = FALSE)
cumHaz <- cumHaz[, c(2, 1)]

validHShrink3_mi <- as.data.frame(t(c(-0.0030383568, 0.0003452605, 0.0010034381)))
colnames(validHShrink3_mi) <- c("gender", "age", "BMI")

x <- pred_input_info(model_type = "survival",
                     model_info = validHShrink3_mi,
                     cum_hazard = cumHaz)

pred_validate(
  x = x,
  new_data = patient2_histories3,
  binary_outcome = NULL,
  survival_time = "time",
  event_indicator = "status",
  time_horizon = 15,
  cal_plot = TRUE
)

# ROC plots

patient2_histories3$predhs3 <-
  predict(validHShrink3, data = patient2_histories3, type = "survival")

real3 <- coxph(Surv(time) ~ 1, data = patients2_histories3)

patient2_histories3$realprobs <-
  predict(real3, data = patient2_histories3, type = "survival")

rochs3 <-
  timeROC(
    T = patient2_histories3$time[1:200],
    cause = 1,
    delta = patient2_histories3$status[1:200],
    marker = patient2_histories3$realprobs[1:200],
    times = seq(0, 15, by = 0.5),
    iid = TRUE
  )

##############################################################################
##               Comparative methods study on Penalisation                  ##
##                    Methods for Multi-State Models                        ##
##                      Data Simulation for 5 State                         ##
##                       Author: Chantelle Cornett                          ##
##                           Date: 24NOV2023                                ##
##############################################################################

##############################################################################
## Change Log:                                                              ##
## Person                |   Date    | Changes Made                         ##
## _________________________________________________________________________##
## Chantelle Cornett     | 24NOV2023 | File initialisation                  ##
##############################################################################

######################################
##          DATA SIMULATION         ##
######################################
library("msm")
library("tidyr")

# Will do this via forking
library("parallel")

detectCores() # 10 cores

n.cohort <- 200000
n <- n.cohort
max.follow <- 365

### Generate the data, but paralellise the process to improve speed (note I could move to CSF at a later date if this process is highly time consuming,
### but at the moment I think other things require the CSF more (fitting the actual multistate models))
print("START DATA GEN")
cl <- makeCluster(5)
registerDoParallel(5)
start_time_fit_parallel <- Sys.time()

data_parallel_list<-(foreach(input=1:5, .combine=list, .multicombine=TRUE, 
                             .packages=c("gems", "dplyr")) %dopar%{
                               
                               x.baseline <- data.frame("BMI" = rnorm(n.cohort, 24, 4),
                                                        "age" = sample(c(1,2,3),n.cohort, replace = TRUE, prob = c(0.332,0.531,0.137 )),
                                                        "gender" = sample(c(0,1),n.cohort, replace = TRUE, prob = c(0.5,0.5)))
                               bl <- x.baseline
                               
                               ### Baseline hazards
                               shape12 <- 1
                               scale12 <- 0.295
                               
                               shape23 <- 1
                               scale23 <- 0.310
                               
                               shape34 <- 1
                               scale34 <- 0.333
                               
                               shape45 <- 1
                               scale45 <- 0.375
                               
                               shape15 <- 1
                               scale15 <- 0.248
                               
                               shape25 <- 1
                               scale25 <- 0.258
                               
                               shape35 <- 1
                               scale35 <- 0.271
                               
                               ## Covariate effects
                               
                               beta12.x1 <- 0.5
                               beta12.x2 <- 0.5
                               beta12.x3 <- 0.5
                               
                               beta23.x1 <- 0.5
                               beta23.x2 <- 0.5
                               beta23.x3 <- 0.5
                               
                               beta34.x1 <- 0.5
                               beta34.x2 <- 0.5
                               beta34.x3 <- 0.5
                               
                               beta45.x1 <- 0.5
                               beta45.x2 <- 0.5
                               beta45.x3 <- 0.5
                               
                               beta15.x1 <- 0.5
                               beta15.x2 <- 0.5
                               beta15.x3 <- 0.5
                               
                               beta25.x1 <- 0.5
                               beta25.x2 <- 0.5
                               beta25.x3 <- 0.5
                               
                               beta35.x1 <- 0.5
                               beta35.x2 <- 0.5
                               beta35.x3 <- 0.5
                               
                               x.in <- x.baseline
                               numsteps <- max.follow
                               ## Generate an empty hazard matrix
                               hf <- generateHazardMatrix(5)
                               #hf
                               
                               ## Change the entries of the transitions we want to allow
                               ## Define the transitions as weibull
                               hf[[1, 2]] <- function(t, shape, scale, beta.x1, beta.x2, beta.x3) {
                                 exp(bl["age"]*beta.x1 - bl["gender"]*beta.x2 + bl["BMI"]*beta.x3)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
                               
                               hf[[2, 3]] <- function(t, shape, scale, beta.x1, beta.x2, beta.x3) {
                                 exp(bl["age"]*beta.x1 - bl["gender"]*beta.x2 + bl["BMI"]*beta.x3)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
                               
                               hf[[3, 4]] <- function(t, shape, scale, beta.x1, beta.x2, beta.x3) {
                                 exp(bl["age"]*beta.x1 - bl["gender"]*beta.x2 + bl["BMI"]*beta.x3)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
                               
                               hf[[4, 5]] <- function(t, shape, scale, beta.x1, beta.x2, beta.x3) {
                                 exp(bl["age"]*beta.x1 - bl["gender"]*beta.x2 + bl["BMI"]*beta.x3)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
                               
                               hf[[1, 5]] <- function(t, shape, scale, beta.x1, beta.x2, beta.x3) {
                                 exp(bl["age"]*beta.x1 - bl["gender"]*beta.x2 + bl["BMI"]*beta.x3)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
                               
                               hf[[2, 5]] <- function(t, shape, scale, beta.x1, beta.x2, beta.x3) {
                                 exp(bl["age"]*beta.x1 - bl["gender"]*beta.x2 + bl["BMI"]*beta.x3)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
                               
                               hf[[3, 5]] <- function(t, shape, scale, beta.x1, beta.x2, beta.x3) {
                                 exp(bl["age"]*beta.x1 - bl["gender"]*beta.x2 + bl["BMI"]*beta.x3)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}
                               
                               ## Generate an empty parameter matrix
                               par <- generateParameterMatrix(hf)
                               
                               ## Use the vector of scales in each transition hazard
                               par[[1, 2]] <- list(shape = shape12, scale = scale12, 
                                                   beta.x1 = beta12.x1, beta.x2 = beta12.x2, beta.x3 = beta12.x3)
                               par[[2, 3]] <- list(shape = shape23, scale = scale23, 
                                                   beta.x1 = beta23.x1, beta.x2 = beta23.x2, beta.x3 = beta23.x3)
                               par[[3, 4]] <- list(shape = shape34, scale = scale34, 
                                                   beta.x1 = beta34.x1, beta.x2 = beta34.x2, beta.x3 = beta34.x3)
                               par[[4, 5]] <- list(shape = shape45, scale = scale45, 
                                                   beta.x1 = beta45.x1, beta.x2 = beta45.x2, beta.x3 = beta45.x3)
                               par[[1, 5]] <- list(shape = shape15, scale = scale15, 
                                                   beta.x1 = beta15.x1, beta.x2 = beta15.x2, beta.x3 = beta15.x3)
                               par[[2, 5]] <- list(shape = shape25, scale = scale25, 
                                                   beta.x1 = beta25.x1, beta.x2 = beta25.x2, beta.x3 = beta25.x3)
                               par[[3, 5]] <- list(shape = shape35, scale = scale35, 
                                                   beta.x1 = beta35.x1, beta.x2 = beta35.x2, beta.x3 = beta35.x3)
                               
                               ## Generate the cohort
                               
                               cohort <- simulateCohort(transitionFunctions = hf, parameters = par,
                                                        cohortSize = n, baseline = bl, to = max.follow, sampler.steps = numsteps)
                               cohort.out <- data.frame(cohort@time.to.state, cohort@baseline, patid = 1:nrow(cohort@time.to.state))
                               
                               ## Turn event times into a dataframe and make the colnames not have any spaces in them
                               dat.mstate.temp <- select(cohort.out, paste("State.", 1:5, sep = ""))
                               colnames(dat.mstate.temp) <- paste0("state", 1:5)
                               
                               ## Now set any transitions that didn't happen to the maximum value of follow up
                               ## Therefore any event that happens, will happen before this. If a transition never happens, an individual will be censored
                               ## at this point in time (when follow up stops)
                               dat.mstate.temp.noNA <- dat.mstate.temp
                               dat.mstate.temp.noNA <- data.frame(t(apply(dat.mstate.temp.noNA, 1, function(x) {ifelse(is.na(x), n.cohort, x)})))
                               
                               
                               ## Add censoring variables
                               dat.mstate.temp.noNA[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2)] <- matrix(0, ncol = ncol(dat.mstate.temp), nrow = nrow(dat.mstate.temp))
                               
                               colnames(dat.mstate.temp.noNA)[6:10] <- 
                                 paste0("state", 1:5, ".s")
                               ###################################################################
                               ##                  UPDATED TO HERE                              ##
                               ###################################################################
                               
                               ## If it is not an NA value (from original dataset), set the censoring indicator to 1
                               dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,2]),(2+ncol(dat.mstate.temp))] <- 1
                               dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,3]),(3+ncol(dat.mstate.temp))] <- 1
                               dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,4]),(4+ncol(dat.mstate.temp))] <- 1
                               dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,5]),(5+ncol(dat.mstate.temp))] <- 1
                               
                               ## Rename dataset to what it was before, and remove excess dataset
                               dat.mstate.temp <- dat.mstate.temp.noNA
                               
                               ## Now need to add baseline data
                               dat.mstate.temp$age <- cohort.out$age
                               dat.mstate.temp$gender <- cohort.out$gender
                               dat.mstate.temp$BMI <- cohort.out$BMI
                               dat.mstate.temp$patid <- cohort.out$patid
                               
                               ### Now we can use msprep from the mstate package to turn into wide format
                               ## First create a transition matrix corresponding to the columns
                               tmat <- matrix(c(NA,NA,NA,NA,NA,1,NA,NA,NA,NA,2,5,NA,NA,NA,3,6,8,NA,NA,4,7,9,10,NA),
                                              nrow=5, ncol = 5)
                               
                               ## Now can prepare the data into wide format
                               dat.mstate.temp.wide <- msprep(dat.mstate.temp, trans = tmat, 
                                                              time = c(NA, paste0("state", 2:5)),
                                                              status = c(NA, paste0("state", 2:5, ".s")), 
                                                              keep = c("age","gender", "BMI","patid"))
                               
                               ## Want to expand the covariates to allow different covariate effects per transition
                               covs <- c("age", "gender", "BMI")
                               dat.mstate.temp.wide <- expand.covs(dat.mstate.temp.wide, covs, longnames = FALSE)
                               dat.mstate.temp.wide
                               
                             })

end_time_fit_parallel <- Sys.time()
diff_fit_parallel <- start_time_fit_parallel - end_time_fit_parallel
stopCluster(cl)
print("FINISH DATA GEN")
diff_fit_parallel

### Combine data into one data frame
temp.data.cohort <- rbind(data_parallel_list[[1]],data_parallel_list[[2]],data_parallel_list[[3]],data_parallel_list[[4]],data_parallel_list[[5]])
### Assign patient ID's and remove rownames
temp.data.cohort$patid <- 1:nrow(temp.data.cohort)
rownames(temp.data.cohort) <- NULL
save(temp.data.cohort, file = "5state.csv")

##############################################################################
##               Comparative methods study on Penalisation                  ##
##                    Methods for Multi-State Models                        ##
##                            Data Simulation                               ##
##                       Author: Chantelle Cornett                          ##
##                           Date: 12NOV2023                                ##
##############################################################################

##############################################################################
## Change Log:                                                              ##
## Person                |   Date    | Changes Made                         ##
## _________________________________________________________________________##
## Chantelle Cornett     | 12NOV2023 | File initialisation                  ##
##                       | 14NOV2023 | First simulation                     ##
##                       | 16NOV2023 | Using Alex's code                    ##
##                       | 19NOV2023 | Using MSM to simulate data.          ##
##                       | 24NOV2023 | Add Weibull and parallelise.         ##
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
scale12 <- 0.219

shape13 <- 1
scale13 <- 0.219

shape23 <- 1
scale23 <- 0.233

## Covariate effects
beta12.x1 <- 0.5
beta12.x2 <- 0.5
beta12.x3 <- 0.5
beta13.x1 <- 0.5
beta13.x2 <- 0.5
beta13.x3 <- 0.5
beta23.x1 <- 0.5
beta23.x2 <- 0.5
beta23.x3 <- 0.5
x.in <- x.baseline
numsteps <- max.follow
## Generate an empty hazard matrix
hf <- generateHazardMatrix(3)
#hf

## Change the entries of the transitions we want to allow
## Define the transitions as weibull
hf[[1, 2]] <- function(t, shape, scale, beta.x1, beta.x2, beta.x3) {
  exp(bl["age"]*beta.x1 - bl["gender"]*beta.x2 + bl["BMI"]*beta.x3)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}

hf[[1, 3]] <- function(t, shape, scale, beta.x1, beta.x2, beta.x3) {
  exp(bl["age"]*beta.x1 - bl["gender"]*beta.x2 + bl["BMI"]*beta.x3)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}

hf[[2, 3]] <- function(t, shape, scale, beta.x1, beta.x2, beta.x3) {
  exp(bl["age"]*beta.x1 - bl["gender"]*beta.x2 + bl["BMI"]*beta.x3)*(shape/scale)*((t + sum(history))/scale)^(shape - 1)}

## Generate an empty parameter matrix
par <- generateParameterMatrix(hf)

## Use the vector of scales in each transition hazard
par[[1, 2]] <- list(shape = shape12, scale = scale12, 
                    beta.x1 = beta12.x1, beta.x2 = beta12.x2, beta.x3 = beta12.x3)
par[[1, 3]] <- list(shape = shape13, scale = scale13, 
                    beta.x1 = beta13.x1, beta.x2 = beta13.x2, beta.x3 = beta13.x3)
par[[2, 3]] <- list(shape = shape23, scale = scale23, 
                    beta.x1 = beta23.x1, beta.x2 = beta23.x2, beta.x3 = beta23.x3)

## Generate the cohort

cohort <- simulateCohort(transitionFunctions = hf, parameters = par,
                         cohortSize = n, baseline = bl, to = max.follow, sampler.steps = numsteps)
cohort.out <- data.frame(cohort@time.to.state, cohort@baseline, patid = 1:nrow(cohort@time.to.state))
  
  ## Turn event times into a dataframe and make the colnames not have any spaces in them
  dat.mstate.temp <- select(cohort.out, paste("State.", 1:3, sep = ""))
  colnames(dat.mstate.temp) <- paste0("state", 1:3)
  
  ## Now set any transitions that didn't happen to the maximum value of follow up
  ## Therefore any event that happens, will happen before this. If a transition never happens, an individual will be censored
  ## at this point in time (when follow up stops)
  dat.mstate.temp.noNA <- dat.mstate.temp
  dat.mstate.temp.noNA <- data.frame(t(apply(dat.mstate.temp.noNA, 1, function(x) {ifelse(is.na(x), n.cohort, x)})))

  
  ## Add censoring variables
  dat.mstate.temp.noNA[(ncol(dat.mstate.temp)+1):(ncol(dat.mstate.temp)*2)] <- matrix(0, ncol = ncol(dat.mstate.temp), nrow = nrow(dat.mstate.temp))

  colnames(dat.mstate.temp.noNA)[4:6] <- 
    paste0("state", 1:3, ".s")
  
  
  ## If it is not an NA value (from original dataset), set the censoring indicator to 1
    dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,2]),(2+ncol(dat.mstate.temp))] <- 1
    dat.mstate.temp.noNA[!is.na(dat.mstate.temp[,3]),(3+ncol(dat.mstate.temp))] <- 1
  
  ## Rename dataset to what it was before, and remove excess dataset
  dat.mstate.temp <- dat.mstate.temp.noNA
  
  ## Now need to add baseline data
  dat.mstate.temp$age <- cohort.out$age
  dat.mstate.temp$gender <- cohort.out$gender
  dat.mstate.temp$BMI <- cohort.out$BMI
  dat.mstate.temp$patid <- cohort.out$patid
  
  ### Now we can use msprep from the mstate package to turn into wide format
  ## First create a transition matrix corresponding to the columns
  tmat <- trans.illdeath()

  ## Now can prepare the data into wide format
  dat.mstate.temp.wide <- msprep(dat.mstate.temp, trans = tmat, 
                                 time = c(NA, paste0("state", 2:3)),
                                 status = c(NA, paste0("state", 2:3, ".s")), 
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

save(temp.data.cohort, file="3state.csv")

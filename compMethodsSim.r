##############################################################################
##               Comparative methods study on Penalisation                  ##
##                    Methods for Multi-State Models                        ##
##                       Author: Chantelle Cornett                          ##
##                           Date: 12NOV2023                                ##
##############################################################################

##############################################################################
## Change Log:                                                              ##
## Person                |   Date    | Changes Made                         ##
## _________________________________________________________________________##
## Chantelle Cornett     | 12NOV2023 | File initialisation                  ##
##                       | 14NOV2023 | First simulation                     ##
##############################################################################

##############################################################################
## To Do:                                                                   ##
## * Simulate patient information                                           ##
##     * Patient ID (DONE)                                                  ##
##.    * Gender (DONE)                                                      ##
##.    * Year of birth (DONE)                                               ##
##.    * Socioeconomic status (DONE)                                        ##
##.    * Region                                                             ##
##.    * Event info                                                         ##
##############################################################################

######################################
##          DATA SIMULATION         ##
######################################

# Install and load the mstate package
install.packages("mstate")
install.packages("gems")
library("mstate")

########## Patient ID ################

patid <- seq(1,2000,1)

############# Gender #################

gender <- list()
for (i in 1:length(patid)){
  u <- runif(1)
  if (u < 0.5){
    gender[i] <- 0
  } 
  else{
    gender[i] <- 1
  }
}

######## Year of Birth ###############

yob <- list()
for (i in 1:length(patid)){
  u <- 0
  while(u < 1 | u > 100){
    u <- rnorm(1, mean = 46, sd = 14)
  } 
  yob[i] <- 2023 - round(u, digits=0)
}

############ SES #####################

ses <- as.list(sample(c(1,2,3,4,5), size = length(patid), replace = TRUE, 
              prob = c(0.2,0.2,0.2,0.2,0.2)))

############ Region ##################
# 10 SHA and then Wales, NI, Scot

region <- as.list(sample(seq(1,13,1), size = length(patid), replace = TRUE))





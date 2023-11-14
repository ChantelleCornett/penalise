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
##     * Patient ID                                                         ##
##.    * Gender                                                             ##
##.    * Year of birth                                                      ##
##.    * Socioeconomic status                                               ##
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



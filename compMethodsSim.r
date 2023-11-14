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
##############################################################################

######################################
##          DATA SIMULATION         ##
######################################

# Install and load the mstate package
install.packages("mstate")
install.packages("gems")
library("mstate")
library("mstate")

set.seed(123)
tmat <- mstate::transMat(x = list(c(1,2,3), 
                         c(2,3), 
                         c()),
                       names = c("Healthy", "Unwell", "Death"))
print(tmat)

id <- seq(1,100,1)
srt_ages <- rnorm(10, mean = 40, sd = 15)

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
##.    * Region (DONE)                                                      ##
##.    * Event info (date, type)                                            ##
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


########### Events ##################

# Set the seed for reproducibility
set.seed(123)

# Function to simulate a multi-state illness-death model
simulate_multi_state <- function(n, time_max) {
  time <- rep(NA, n)
  state <- rep(NA, n)
  
  for (i in 1:n) {
    t <- 0
    current_state <- 1  # Start in state 1 (e.g., Healthy)
    
    while (t < time_max) {
      # Generate time until the next event
      delta_t <- rexp(1, rate = -Q[current_state, current_state])
      
      # Update time and state
      t <- t + delta_t
      time[i] <- t
      state[i] <- current_state
      
      # Check for death event
      if (current_state == n_states) {
        break
      }
      
      # Generate next state
      u <- runif(1)
      prob_transitions <- cumsum(Q[current_state, ] / -Q[current_state,
                                                         current_state])
      next_state <- sum(u > prob_transitions) + 1
      
      current_state <- next_state
    }
  }
  
  return(data.frame(time = time, state = state))
}

# Define the number of individuals and number of states
n_individuals <- 2000
n_states <- 3

# Define the transition intensity matrix
Q <- matrix(c(-0.1, 0.1, 0, 0, -0.2, 0.2, 0, 0, -0.05), nrow = n_states,
            byrow = TRUE)

# Simulate data
sim_data <- simulate_multi_state(n = n_individuals, time_max = 10)

# Plot the simulated data
plot(sim_data$time, sim_data$state, type = 's', col = 'blue', xlab = 'Time',
     ylab = 'State',
     main = 'Simulated Multi-State Illness-Death Model')


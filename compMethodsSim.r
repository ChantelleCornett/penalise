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
##############################################################################



######################################
##          DATA SIMULATION         ##
######################################
library("msm")
library("doBy")
library("tidyr")

# Transition matrix

tmat <- trans.illdeath()

# need to simulate event times for 100 people

sim.df <-
  data.frame(subject = rep(1:100, rep(13, 100)), time = rexp(100, 0.2)) %>% arrange(subject, time)

# Transition rates
lambda_h_i = 0.02  # Transition rate from healthy to illness
mu_i_d = 0.01      # Transition rate from illness to death

# Q-matrix
qmatrix <- matrix(
  c(-lambda_h_i, lambda_h_i, 0,
    0,-mu_i_d, mu_i_d,
    0, 0, 0),
  nrow = 3,
  byrow = TRUE
)

df <- simmulti.msm(sim.df, qmatrix, death = c(3))

df2 <- summaryBy(time ~ subject + state, FUN = c(sum), data = df)

# Getting data into same format as ebmt3

df3 <-
  pivot_wider(df2, names_from = state, values_from = c(time.sum)) # transpose
df3 <-
  df3[, !(names(df3) %in% c('3'))] # remove 'time' spent in state 3
df3$'2' <- df3$'2'+df3$'1' # cumulative time to state 3
colnames(df3) <- c('subject', 'timeToIll', 'timeToDeath')
df3$timeToDeath <-
  replace_na(df3$timeToDeath, replace = 200) # set missing vals to a large number

# create indicator variables for transitions to state 2 and 3
counts <- as.data.frame(table(df2$subject))
colnames(counts) <- c("subject", "count")
df4 <- merge(df3, counts, by = c("subject"))

for (i in 1:length(df4$subject)) {
  if (df4$count[i] == 1) {
    df4$illstat[i] <- 0
    df4$deathstat[i] <- 0
  } else if (df4$count[i] == 2) {
    df4$illstat[i] <- 1
    df4$deathstat[i] <- 0
  } else if (df4$count[i] == 3) {
    df4$illstat[i] <- 1
    df4$deathstat[i] <- 1
  }
}

# Get data in format suitable for use by mstate

msdat <-
  msprep(
    data = df4,
    id = "subject",
    time = c(NA, "timeToIll", "timeToDeath"),
    trans = tmat,
    status = c(NA, "illstat", "deathstat")
  )


save(msdat, file = "msdat.Rdata")
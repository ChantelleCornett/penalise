library(msm)

# Transition rates matrix
# Rows and columns represent states A, B, C
# A -> B is common, B -> C is common, A -> C is rare
Q <- matrix(c(-0.1, 0.1, 0.00001,  # From A to B (common), A to C (rare)
              0, -0.2, 0.2,      # From B to C (common)
              0, 0, 0),          # Absorbing state C has no transitions out
            byrow = TRUE, nrow = 3)

rownames(Q) <- colnames(Q) <- c("A", "B", "C")

adjust_transition_matrix <- function(Q, age, gender, BMI) {
  # This is a placeholder function. You should replace the logic
  # below with your specific rules for how age, gender, and BMI
  # affect the transition rates.
  
  # Example adjustments:
  # - Increase transition rate from state A to B by 0.01 for each year above 50.
  # - Adjust transition rates for gender (e.g., if gender == "female", decrease certain rates).
  # - Adjust transition rates based on BMI categories.
  
  # Clone the original matrix to avoid altering it directly
  Q_adjusted <- Q
  
  # Adjust rate for age (assuming age affects the first transition rate as an example)
  if (age > 50) {
    Q_adjusted[1,2] <- Q[1,2] * (1 + 0.01 * (age - 50))
  }
  
  # Example: Adjust rate based on gender (placeholder logic)
  if (gender == 1) {
    Q_adjusted[1,2] <- Q_adjusted[1,2] * 0.9 # Example adjustment
  }
  
  # Example: Adjust rates based on BMI (simplified logic)
  if (BMI > 30) { # Classifying BMI > 30 as obesity
    Q_adjusted[1,2] <- Q_adjusted[1,2] * 1.1 # Increase rate by 10% for obesity
  }
  
  return(Q_adjusted)
}
patient_data <- data.frame(
  id = 1:200000,  # Adjust for the actual number of patients
  age = runif(200000, 18, 85),  # Random ages for demonstration
  gender = sample(c(0, 1), 200000, replace = TRUE),
  BMI = runif(200000, 18, 40)  # Random BMI for demonstration
)

# Initialize an empty list for storing detailed patient trajectories
patient_list <- vector("list", length = nrow(patient_data))

set.seed(123)
for(i in 1:nrow(patient_data)) {
  patient <- patient_data[i,]
  
  # Adjust the transition matrix Q based on patient covariates
  Q_adjusted <- adjust_transition_matrix(Q, patient$age, patient$gender, patient$BMI)
  
  # Simulate the trajectory for patient i with the adjusted transition matrix
  simulated <- sim.msm(qmatrix = Q_adjusted, start = 1, maxtime= 200)
  
  # Store the results as before
  patient_list[[i]] <- data.frame(
    patient_id = patient$id,
    from = head(simulated$state, -1),
    to = tail(simulated$state, -1),
    time = diff(simulated$times),
    cumulative_time = head(simulated$times, -1),
    age = patient$age,
    gender = patient$gender,
    BMI = patient$BMI
  )
}

# Combine all patient data frames into one
patient_histories <- do.call(rbind, patient_list)

nrow(subset(patient_histories, from == 1 & to == 3))

patient_histories1 <- subset(patient_histories, from == 1 & to == 2)
patient_histories2 <- subset(patient_histories, from == 1 & to == 3)
patient_histories3 <- subset(patient_histories, from == 2 & to == 3)


set.seed(1)

patient_data2 <- data.frame(
  id = 1:200000,  # Adjust for the actual number of patients
  age = runif(200000, 18, 85),  # Random ages for demonstration
  gender = sample(c(0, 1), 200000, replace = TRUE),
  BMI = runif(200000, 18, 40)  # Random BMI for demonstration
)

# Initialize an empty list for storing detailed patient trajectories
patient_list2 <- vector("list", length = nrow(patient_data))

for(i in 1:nrow(patient_data2)) {
  patient <- patient_data2[i,]
  
  # Adjust the transition matrix Q based on patient covariates
  Q_adjusted <- adjust_transition_matrix(Q, patient$age, patient$gender, patient$BMI)
  
  # Simulate the trajectory for patient i with the adjusted transition matrix
  simulated <- sim.msm(qmatrix = Q_adjusted, start = 1, maxtime= 200)
  
  # Store the results as before
  patient_list2[[i]] <- data.frame(
    patient_id = patient$id,
    from = head(simulated$state, -1),
    to = tail(simulated$state, -1),
    time = diff(simulated$times),
    cumulative_time = head(simulated$times, -1),
    age = patient$age,
    gender = patient$gender,
    BMI = patient$BMI
  )
}

patient2_histories <- do.call(rbind, patient_list2)

nrow(subset(patient2_histories, from == 1 & to == 3))

patient2_histories1 <- subset(patient2_histories, from == 1 & to == 2)
patient2_histories2 <- subset(patient2_histories, from == 1 & to == 3)
patient2_histories3 <- subset(patient2_histories, from == 2 & to == 3)


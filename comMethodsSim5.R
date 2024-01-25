library("mstate")
devtools::install_github("hputter/mstate")
library("markovMSM")
tmat <- mstate::transMat(x = list(c(2, 5), 
                                  c(3, 5), 
                                  c(4, 5), 
                                  c(5),
                                  c()),
                         names = c("Healthy", "I", "II", "III", 
                                    "D"))
print(tmat)

data(ebmt4)

id <- seq(1,200000,1)
H <- rep(times =200000, 0)
H.s <- rep(times = 200000,1)
I <- rweibull(200000,shape = 2, scale = 56.42)
I.s <- sample(size=200000, replace = TRUE, prob = c(0.95,0.05), x=c(1,0))

II <- c()
II.s <- c()
for (i in 1:200000){
  if(I.s[i] == 1){
    II[i] <- I[i] + rweibull(1,shape = 2, scale = 22.6)
    II.s[i]<- sample(size=1, replace = TRUE, prob = c(0.95,0.05), x=c(1,0))
  } else{
    II[i] <- I[i]
    II.s[i] <- 0
  }
}

III <- c()
III.s <- c()
for (i in 1:200000){
  if(II.s[i] == 1){
    III[i] <- II[i] + rweibull(1,shape = 2, scale = 22.6)
    III.s[i]<- sample(size=1, replace = TRUE, prob = c(0.5,0.5), x=c(1,0))
  } else{
    III[i] <- II[i]
    III.s[i] <- 0
  }
}

D <- c()
D.s <- c()
for (i in 1:200000){
  if(III.s[i] == 1){
    D[i] <- III[i] + rweibull(1,shape = 2, scale = 22.6)
   D.s[i]<- sample(size=1, replace = TRUE, prob = c(0.95,0.05), x=c(1,0))
  } else{
    D[i] <- III[i]
    D.s[i] <- 0
  }
}

test <- cbind(H, H.s,I,I.s,II,II.s, III, III.s, D, D.s)

gender <- rbinom(200000, 1, .5)
age <- sample(c(1,2,3), 200000, replace=TRUE)
BMI <- rnorm(200000, mean = 24, sd = 3)

pre5state <- cbind(id,test, gender, age, BMI)
new <- as.data.frame(pre5state)
state5 <- msprep(data=as.data.frame(new), trans = tmat, time = c("H","I","II","III","D"), status = c("H.s","I.s","II.s","III.s","D.s"),id="id")

# merge characteristics back on

allstate5 <- merge(pre5state[,c(1,12,13,14)], state5,by=c("id"))

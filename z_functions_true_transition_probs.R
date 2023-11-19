### Storing these functions in a seperate workspace, so they don't clog up the other function file,
### as these are very long, and have lots of functions within them

###
###
### Function to calculate true transition probabilities for DGM1
###
###
calc.true.transition.probs.DGM1 <- function(u.eval, t.eval, x1.eval, x2.eval,
                                       shape12, scale12, #shape and scale for weibull baseline hazard for transition 1 -> 2
                                       shape13, scale13, #shape and scale for weibull baseline hazard for transition 1 -> 3
                                       shape23, scale23, #shape and scale for weibull baseline hazard for transition 2 -> 3
                                       beta12.x1, beta12.x2, #covariate effects for transiion 12
                                       beta13.x1, beta13.x2, #covariate effects for transiion 13
                                       beta23.x1, beta23.x2, #covariate effects for transiion 23
                                       ){
  ###
  ### Start by defining the cause-specific hazards
  csh12 <- function(t.in, x1, x2){
    return(exp(beta12.x1*x1 + beta12.x2*x2)*(shape12/scale12)*((t.in)/scale12)^(shape12 - 1))
  }
  
  csh13 <- function(t.in, x1, x2){
    return(exp(beta13.x1*x1 + beta13.x2*x2)*(shape13/scale13)*((t.in)/scale13)^(shape13 - 1))
  }
  
  csh23 <- function(t.in, x1, x2){
    return(exp(beta23.x1*x1 + beta23.x2*x2)*(shape23/scale23)*((t.in)/scale23)^(shape23 - 1))
  }
  }
  
  ###
  ### Define survival functions out of each state
  S1 <- function(t.in, x1, x2){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta12.x1*x1 + beta12.x2*x2)*(t.in/scale12)^shape12 +
      exp(beta13.x1*x1 + beta13.x2*x2)*(t.in/scale13)^shape13
    ## Convert into survival probability
    return(exp(-(cumhaz)))
  }
  
  S2 <- function(t.in, x1, x2){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta23.x1*x1 + beta23.x2*x2)*(t.in/scale23)^shape23
    ## Convert into survival probabilty
    return(exp(-(cumhaz)))
  }
  
  
  ###
  ###
  ### Define transition probabilities
  ###
  ###
  
  ### In format P.ij.route (from state i, to state j, via route)
  
  ###
  ### Transitions out of state 2
  ###
  P.22 <- function(u.in, t.in, x1, x2){
    return(S2(t.in, x1, x2)/S2(u.in, x1, x2))
  }
  
  P.23.direct <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh23(r.in, x1, x2)*S2(r.in, x1, x2)/S2(u.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ###
  ### Transitions out of state 1
  ###
  
  ###
  ### Staying in state 1
  P.11 <- function(u.in, t.in, x1, x2){
    return(S1(t.in, x1, x2)/S1(u.in, x1, x2))
  }
  
  ###
  ### Transitions to states 2 or 3 (and staying)
  P.12 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh12(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.22(r.in, t.in, x1, x2))
    }
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  P.13 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh13(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.33(r.in, t.in, x1, x2))
    }
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  
  ###
  ### Transitions to state 3
  
  ## Direct
  P.13.direct <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh13(r.in, x1, x2)*S1(r.in, x1, x2)/S1(u.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Via 2 only
  P.13.2 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh12(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.23.direct(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
 
  
  ## Combine all for transition probability from 1 to 3
  P.13 <- function(u.in, t.in, x1, x2){
    return(P.13.direct(u.in, t.in, x1, x2) + 
             P.13.2(u.in, t.in, x1, x2)
            )
  }
  
  ### Create output vector with transition probabilities
  output <- c("P.11" = P.11(u.eval, t.eval, x1.eval, x2.eval),
              "P.12" = P.12(u.eval, t.eval, x1.eval, x2.eval),
              "P.13" = P.13(u.eval, t.eval, x1.eval, x2.eval))
  return(output)


# 
# ### Storing these seperately so they don't clog up workspace, and make z_functions.R file difficult to navigate
# 
# ### Baseline hazards
# shape12 <- 1
# scale12 <- 1588.598
# 
# shape13 <- 1
# scale13 <- 0.5*1588.598
# 
# shape15 <- 1
# scale15 <- 5*1588.598
# 
# shape24 <- 1
# scale24 <- 1588.598
# 
# shape25 <- 1
# scale25 <- 5*1588.598
# 
# shape34 <- 1
# scale34 <- 0.5*1588.598
# 
# shape35 <- 1
# scale35 <- 5*1588.598
# 
# shape45 <- 1
# scale45 <- 5*1588.598
# 
# #qweibull(0.8, 1, 1588.598)
# 
# ## Covariate effects
# beta12.x1 <- 1
# beta12.x2 <- 1
# beta13.x1 <- 0.5
# beta13.x2 <- 0.5
# beta15.x1 <- 1
# beta15.x2 <- 0.5
# beta24.x1 <- 0.5
# beta24.x2 <- 1
# beta25.x1 <- 1
# beta25.x2 <- 1
# beta34.x1 <- 0.5
# beta34.x2 <- 0.5
# beta35.x1 <- 1
# beta35.x2 <- 0.5
# beta45.x1 <- 0.5
# beta45.x2 <- 1
# 
# ### Let's test this
# true.trans.probs <- calc.true.transition.probs.DGM1(u.eval = 0, t.eval = ceiling(7*365.25), x1.eval = 0, x2.eval = 0,
#                        shape12 = 1, scale12 = 1588.598, #shape and scale for weibull baseline hazard for transition 1 -> 2
#                        shape13 = 1, scale13 = 0.5*1588.598, #shape and scale for weibull baseline hazard for transition 1 -> 3
#                        shape15 = 1, scale15 = 5*1588.598, #shape and scale for weibull baseline hazard for transition 1 -> 5
#                        shape24 = 1, scale24 = 1588.598, #shape and scale for weibull baseline hazard for transition 2 -> 4
#                        shape25 = 1, scale25 = 5*1588.598, #shape and scale for weibull baseline hazard for transition 2 -> 5
#                        shape34 = 1, scale34 = 0.5*1588.598, #shape and scale for weibull baseline hazard for transition 3 -> 4
#                        shape35 = 1, scale35 = 5*1588.598, #shape and scale for weibull baseline hazard for transition 3 -> 5
#                        shape45 = 1, scale45 = 5*1588.598, #shape and scale for weibull baseline hazard for transition 3 -> 5
#                        beta12.x1 = 1, beta12.x2 = 1, #covariate effects for transiion 12
#                        beta13.x1 = 0.5, beta13.x2 = 0.5, #covariate effects for transiion 13
#                        beta15.x1 = 1, beta15.x2 = 0.5, #covariate effects for transiion 15
#                        beta24.x1 = 0.5, beta24.x2 = 1, #covariate effects for transiion 24
#                        beta25.x1 = 1, beta25.x2 = 1, #covariate effects for transiion 25
#                        beta34.x1 = 0.5, beta34.x2 = 0.5, #covariate effects for transiion 34
#                        beta35.x1 = 1, beta35.x2 = 0.5, #covariate effects for transiion 35
#                        beta45.x1 = 0.5, beta45.x2 = 1 #covariate effects for transiion 45
#                        )
# true.trans.probs
# 


###
###
### Function to calculate true transition probabilities for DGM2
###
###
calc.true.transition.probs.DGM2 <- function(u.eval, t.eval, x1.eval, x2.eval, out.num.states,
                                            shape12, scale12, #shape and scale for weibull baseline hazard for transition 1 -> 2
                                            shape13, scale13, #shape and scale for weibull baseline hazard for transition 1 -> 3
                                            shape23, scale23, #shape and scale for weibull baseline hazard for transition 2 -> 3
                                            beta12.x1, beta12.x2, #covariate effects for transiion 1 -> 2
                                            beta13.x1, beta13.x2, #covariate effects for transiion 1 -> 3
                                            beta23.x1, beta23.x2, #covariate effects for transiion 2 -> 3
){
  ###
  ### Start by defining the cause-specific hazards
  csh12 <- function(t.in, x1, x2){
    return(exp(beta12.x1*x1 + beta12.x2*x2)*(shape12/scale12)*((t.in)/scale12)^(shape12 - 1))
  }
  
  csh13 <- function(t.in, x1, x2){
    return(exp(beta13.x1*x1 + beta13.x2*x2)*(shape13/scale13)*((t.in)/scale13)^(shape13 - 1))
  }
  
  
  csh23 <- function(t.in, x1, x2){
    return(exp(beta23.x1*x1 + beta23.x2*x2)*(shape23/scale23)*((t.in)/scale23)^(shape23 - 1))
  }
  
}
  
  ###
  ### Define survival functions out of each state
  S1 <- function(t.in, x1, x2){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta12.x1*x1 + beta12.x2*x2)*(t.in/scale12)^shape12 +
      exp(beta13.x1*x1 + beta13.x2*x2)*(t.in/scale13)^shape13
    ## Convert into survival probabilty
    return(exp(-(cumhaz)))
  }
  
  S2 <- function(t.in, x1, x2){
    ## Calculate total cumulative hazard out of state
    cumhaz <- exp(beta23.x1*x1 + beta23.x2*x2)*(t.in/scale24)^shape23
    ## Convert into survival probabilty
    return(exp(-(cumhaz)))
  }
  
  
  ###
  ###
  ### Define transition probabilities
  ###
  ###
  
  ### In format P.ij.route (from state i, to state j, via route)
  
  ###
  ### Transitions out of state 2
  ###
  P.22 <- function(u.in, t.in, x1, x2){
    return(S2(t.in, x1, x2)/S2(u.in, x1, x2))
  }
  
  P.23.direct <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh23(r.in, x1, x2)*S2(r.in, x1, x2)/S2(u.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ###
  ### Transitions out of state 1
  ###
  
  ###
  ### Staying in state 1
  P.11 <- function(u.in, t.in, x1, x2){
    return(S1(t.in, x1, x2)/S1(u.in, x1, x2))
  }
  
  ###
  ### Transitions to states 2 or 3 (and staying)
  P.12 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh12(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.22(r.in, t.in, x1, x2))
    }
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  P.13 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh13(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.33(r.in, t.in, x1, x2))
    }
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
 
  
  ###
  ### Transitions to state 3
  
  ## Direct
  P.13.direct <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh13(r.in, x1, x2)*S1(r.in, x1, x2)/S1(u.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
  ## Via 2 only
  P.13.2 <- function(u.in, t.in, x1, x2){
    
    ### First define the function we want to integrate over (we integrate it between u and t)
    func.for.int <- function(r.in){
      return(csh12(r.in, x1, x2)*(S1(r.in, x1, x2)/S1(u.in, x1, x2))*P.23.direct(r.in, t.in, x1, x2))
    }
    
    ### Do integration
    numerical.int.out <- cubintegrate(f = func.for.int, lower = u.in, upper = t.in, method = "pcubature")
    
    return(numerical.int.out$integral)
  }
  
 
  
  ## Combine all for transition probability from 1 to 5
  P.16 <- function(u.in, t.in, x1, x2){
    return(P.16.direct(u.in, t.in, x1, x2) + 
             P.16.2(u.in, t.in, x1, x2) +
             P.16.24(u.in, t.in, x1, x2) + 
             P.16.3(u.in, t.in, x1, x2) +
             P.16.35(u.in, t.in, x1, x2))
  }
  
  ### Create output vector with transition probabilities. Have two options, for whether we want the true probs of the underling states,
  ### or the probability of being in the states according to the 5 state model
  if (out.num.states == 5){
    output <- c("P.11" = P.11(u.eval, t.eval, x1.eval, x2.eval),
                "P.12" = P.12(u.eval, t.eval, x1.eval, x2.eval),
                "P.13" = P.13(u.eval, t.eval, x1.eval, x2.eval),
                "P.14" = P.14(u.eval, t.eval, x1.eval, x2.eval) + P.15(u.eval, t.eval, x1.eval, x2.eval),
                "P.15" = P.16(u.eval, t.eval, x1.eval, x2.eval))
  } else if (out.num.states == 6){
    output <- c("P.11" = P.11(u.eval, t.eval, x1.eval, x2.eval),
                "P.12" = P.12(u.eval, t.eval, x1.eval, x2.eval),
                "P.13" = P.13(u.eval, t.eval, x1.eval, x2.eval),
                "P.14" = P.14(u.eval, t.eval, x1.eval, x2.eval),
                "P.15" = P.15(u.eval, t.eval, x1.eval, x2.eval),
                "P.16" = P.16(u.eval, t.eval, x1.eval, x2.eval))
  }

  return(output)
}


# ### Storing these seperately so they don't clog up workspace, and make z_functions.R file difficult to navigate
# u.eval <- 0
# t.eval <- ceiling(7*365.25)
# x1.eval <- 0
# x2.eval <- 0
# ### Baseline hazards
# shape12 <- 1
# scale12 <- 1588.598
# 
# shape13 <- 1
# scale13 <- 0.5*1588.598
# 
# shape16 <- 1
# scale16 <- 5*1588.598
# 
# shape24 <- 1
# scale24 <- 1588.598
# 
# shape26 <- 1
# scale26 <- 5*1588.598
# 
# shape35 <- 1
# scale35 <- 0.5*1588.598
# 
# shape36 <- 1
# scale36 <- 5*1588.598
# 
# shape46 <- 1
# scale46 <- 5*1588.598
# 
# shape56 <- 1
# scale56 <- 5*1588.598
# 
# #qweibull(0.8, 1, 1588.598)
# 
# ## Covariate effects
# beta12.x1 <- 1
# beta12.x2 <- 1
# beta13.x1 <- 0.5
# beta13.x2 <- 0.5
# beta16.x1 <- 1
# beta16.x2 <- 0.5
# beta24.x1 <- 0.5
# beta24.x2 <- 1
# beta26.x1 <- 1
# beta26.x2 <- 1
# beta35.x1 <- 0.5
# beta35.x2 <- 0.5
# beta36.x1 <- 1
# beta36.x2 <- 0.5
# beta46.x1 <- 0.5
# beta46.x2 <- 1
# beta56.x1 <- 0.5
# beta56.x2 <- 1
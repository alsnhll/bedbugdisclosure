# ----------------------------------------------------------------------------
# SetODEsDDM function:
# ----------------
# This is the same as the SetODEs function with the addition of multiple compartments
# of disclosed units, which allows for a non-exponential distribution of times spent
# in the disclosed comparment. This code involves the additional parameter
#
#   Dnum = number of disclosed compartments Sv2. Mean time of being in disclosed state 
#         (if no leaving due to renting) is 1/D with standard deviation 1/sqrt(Dnum)
#
# and Sv2 and DSv2.dt are now vectors with Dnum entries

SetODEsDDM<-function(t,y,p){
  # Dnum is number of Sv2 compartments
  Sr <- y[1]
  Ir <- y[2]
  Sv <- y[3]
  Iv <- y[4]
  Sv2 <- y[5:(5+Dnum-1)] #vector
  trt <- y[5+Dnum]
  tov <- y[5+Dnum+1]
  
  # Note f(t) = b*Ir/(Sr+b*Ir)
  with(as.list(p),{
    dSr.dt <- -beta*Sr*Ir/N + gamma*Ir + n*(1-d)*(1-k*b*Ir/(Sr+b*Ir))*sum(Sv2) + n*(1-k*b*Ir/(Sr+b*Ir))*Sv - m*Sr
    dIr.dt <- beta*Sr*Ir/N + n*k*b*Ir/(Sr+b*Ir)*Sv + n*(1-d)*k*b*Ir/(Sr+b*Ir)*sum(Sv2) + n*(1-d)*Iv - gamma*Ir - b*m*Ir 
    dSv.dt <- m*Sr + Dnum/D*Sv2[Dnum] - n*Sv 
    dIv.dt <- b*m*Ir - gamma*Iv - n*(1-d)*Iv
    
    dSv2.dt <- numeric(Dnum)
    dSv2.dt[1] <- gamma*Iv - n*(1-d)*Sv2[1] - Dnum/D*Sv2[1]
    
    if(Dnum>1){
      for(j in 2:Dnum){
        dSv2.dt[j] <- Dnum/D*Sv2[j-1] - n*(1-d)*Sv2[j] - Dnum/D*Sv2[j]
      }
    }
    
    dtrt.dt <- gamma*Ir + gamma*Iv
    dtov.dt <- n*Sv + n*(1-d)*sum(Sv2) + n*(1-d)*Iv
    dvac.dt <- Sv + sum(Sv2) +Iv
    
    return(list(c(dSr.dt, dIr.dt, dSv.dt, dIv.dt, dSv2.dt, dtrt.dt, dtov.dt, dvac.dt)))
  })
}

# ----------------------------------------------------------------------------
# GetCostDDM function:
# ----------------
# Inputs: 
#   - Vector of parameter values & initial conditions
#   - Vector of bed bug-related costs
#   - Vector of years over which to run the simulation
# Output: 
#   - A data frame tracking the total and component costs, as well as 
#     prevalence for each year of the simulation
#   - This output is similar to that of Figure 3 but for the model for fixed disclosure period

GetCostDDM <- function(p.set, bbcosts, years){
  
  # Set parameter values
  gamma <- p.set[1]
  k <- p.set[2]
  b <- p.set[3]
  m <- p.set[4]
  n <- p.set[5]
  N <- p.set[6]
  D <- p.set[7]
  beta <- p.set[8]
  d <- p.set[14]
  Dnum <- p.set[15]
  
  # Set initial conditions and time interval
  Sr0 <- p.set[10]
  Ir0 <- p.set[11]
  Sv0 <- p.set[12]
  Iv0 <- p.set[13]
  Sv20 <- 0  # Sv20 is set to 0 because we assume disclosure begins at time 0
  trt0 <- 0  # set treatment counter to 0 at time 0
  tov0 <- 0  # set turnover counter to 0 at time 0
  vac0 <- 0  # set vacancy counter to 0 at time 0
  y0 <- c(Sr0, Ir0, Sv0, Iv0, Sv20*rep(0,Dnum), trt0, tov0, vac0)
  y00 <- c(Sr0, Ir0, Sv0, Iv0, Sv20, trt0, tov0, vac0)
  t <- seq(from=0, to=365*max(years)+1, by=1)
  
  # Set bed bug related costs (treatment, turnover, and vacancy costs)
  ctrt <- bbcosts[1] 
  ctov <- bbcosts[2]
  cvac <- bbcosts[3]
  
  # We model the absence of disclosure by setting the disclosure index d = 0
  p <- list(beta=beta, gamma=gamma, b=b, d=d, D=D, k=k, m=m, n=n, N=N, Dnum=Dnum)
  p0 <- list(beta=beta, gamma=gamma, b=b, d=0, D=D, k=k, m=m, n=n, N=N)
  
  # We solve the ODE's for the model in the presence and absence of disclosure 
  # (out and out0, respectively)
  out <- ode(y=y0, times=t, func=SetODEsDDM, parms=p)
  out0 <- ode(y=y00, times=t, func=SetODEs, parms=p0)
  
  # get total number in disclosed compartment
  if(Dnum>1){
    Sv2tot <- rowSums(out[,6:(6+Dnum-1)])
  } else {
    Sv2tot <- out[,6]
  }
  
  colf <- dim(out)[2] #final variable
  
  # We use a for loop to calculate total and component costs for each year of 
  # the simulation
  cost <- trt <- tov <- vac <- prev <- pvac <- as.numeric()
  
  for(jj in 1:length(years)){
    
    # Get the first and last day of year jj
    first.day <- years[jj]*365 - 364
    last.day <- years[jj]*365 + 1
    
    # Calculate the DIFFERENCE between the disclosure and no disclosure simulations 
    # in the number of treatments (n.trt), turnovers (n.tov), and days vacant (n.vac) 
    # that fell on year jj
    n.trt <- ((out[,colf-2][last.day] - out[,colf-2][first.day]) - 
                (out0[,7][last.day] - out0[,7][first.day]))
    n.tov <- ((out[,colf-1][last.day] - out[,colf-1][first.day]) - 
                (out0[,8][last.day] - out0[,8][first.day]))
    n.vac <- ((out[,colf][last.day] - out[,colf][first.day]) - 
                (out0[,9][last.day] - out0[,9][first.day]))
    
    # Total per unit treatment cost = 
    # (# of treatments) x (avg cost of bed bug treatment) / (total # units)
    trt[jj] <- n.trt*ctrt/N
    
    # Total per unit turnover cost = 
    # (# of turnover events) x (avg cost of turnover) / (total # units)
    tov[jj] <- n.tov*ctov/N
    
    # Total per unit vacancy cost = 
    # (# months vacant) x (average monthly rent) / (total # units)
    # Note: # months = # days / 30
    vac[jj] <- (n.vac/30)*cvac/N
    
    # Total cost is equal to the sum of the component costs
    cost[jj] <- trt[jj] + tov[jj] + vac[jj]
    
    # Prevalence at the end of year jj is simply the number of units in the 
    # Ir and Iv classes on the last day of the year divided by N
    prev[jj] <- (out[,3][last.day] + out[,5][last.day])/N
    
    pvac[jj] <- (out[,4][last.day] + out[,5][last.day] + Sv2tot[last.day])/N
    
  }
  
  df <- data.frame(Year = years, Total_Cost = cost, Treatment = trt, 
                   Turnover = tov, Vacancy = vac, Prevalence=prev, 
                   Prop_Vacant=pvac)
  return(df)
}
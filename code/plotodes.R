library(deSolve)

# IMPORTANT: UPDATE THE WORKING DIRECTORY BELOW
#setwd("~/Dropbox/Research/SESYNC/Code")
#setwd("UPDATE_PATH/bedbugdisclosure")

source("code/functions.R")

# 1. Set baseline prevalence (p), renter selectivity (s), and years of simulation
p <- 0.05
s <- 0.5
nyears <- 10

# 2. Set parameter values manually
SetParametersManual <- function(){
  gamma <- 1/(6*30)
  k <- 0.3
  b <- 1.3
  m <- 1/(2*365)
  n <- 1/(3*30)
  N <- 1000
  D <- 365
  return(c(gamma=gamma, k=k, b=b, m=m, n=n, N=N, D=D))
}
preparam <- SetParametersManual()

# 3. Solve for beta that will give the desired baseline prevalence
beta <- GetBeta(preparam, p)

# 4. Either solve for initial conditions "at equilibrium"...
param <- c(preparam, beta, base.prev=p)
init <- GetInit(param)
y0 <- c(init, 0, 0, 0, 0)
names(y0)[5:8] <- c("Sv2","C_treat","C_turnover","C_vacancy")

# OR

# 4. set initial conditions manually
#y0 <- c(Sr0=899, Ir0=1, Sv0=100, Iv0=0, Sv20=0, trt0=0, tov0=0)

# 5. Format parameters and time to input into ode function
pp <- list(beta = beta, gamma = preparam[1], k = preparam[2], b = preparam[3],
           m = preparam[4], n = preparam[5], N = preparam[6], D = preparam[7],
           d = s)
t <- seq(from=0, to=365*nyears+1, by=1)

# 6. Run the ode solver function
out <- ode(y=y0, times = t, func=SetODEs, parms=pp)

# 7. Plot the output

plot(t/365, out[,2], type="l", col="black", ylim=c(0,1000), 
     ylab="Number of units", xlab="Year") #Sr
lines(t/365, out[,4], col="black", lty=3) #Sv
lines(t/365, out[,6], col="black", lty=2) #Sv'
lines(t/365, out[,3], col="red") #Ir
lines(t/365, out[,5], col="red", lty=3) #Iv
legend("right", legend=c("Sr", "Sv","Sv'", "Ir", "Iv"), bty="n",
       col=c("black", "black","black", "red", "red"), lty=c(1,3,2,1,3),lwd=c(1,1,1,1,1), cex=0.8)

Itot <- out[,3]+out[,5]
Vtot <- out[,4]+out[,5]+out[,6]

plot(t/365, Itot, type="l", col="red", lty=1, lwd=3, ylim=c(0,1000), 
     ylab="Number of units", xlab="Year") #Iv+Ir
lines(t/365, Vtot, col="blue", lty=3, lwd=3) #Sv+Sv'+Iv
legend("right", legend=c("Total I","Total v"), bty="n",
       col=c("red","blue"), lty=c(1,3),lwd=c(3,3), cex=0.8)

print(out[dim(out)[1],])
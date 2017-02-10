#################################################################
#
#   Modeling Infectious Diseases in Humans & Animals
#   Chapter 1
#   Sean Wu
#   January 31, 2017
#
#################################################################

library(deSolve)


#Epidemic Burnout
epiBurnout <- function(r0,rInf,s0=1){
  1 - rInf - s0*exp(-rInf*r0)
}

rInfMesh <- seq(.Machine$double.eps,1,length.out=501)[1:500]

r0Mesh <- sapply(rInfMesh,function(x){
  uniroot(epiBurnout,interval=c(1,20),rInf=x)$root
})

plot(r0Mesh,rInfMesh,type="l",col="steelblue",lwd=1.5,xlab=expression("Basic Reproductive Ratio, R"[0]),ylab="Fraction Infected")
grid()

#SIR model with demography
sirDemo <- function(time,state,parms){
  with(as.list(c(state,parms)),{
    dS = mu - beta*S*I - mu*S
    dI = beta*S*I - gamma*I - mu*I
    dR = gamma*I - mu*R
    return(list(c(dS,dI,dR)))
  })
}

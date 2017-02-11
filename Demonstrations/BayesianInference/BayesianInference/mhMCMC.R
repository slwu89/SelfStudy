#################################################################
#
#   Metropolis Hastings Random Walk and Adaptive MCMC
#   R Shiny Helper Functions
#   Sean Wu
#   February 10, 2017
#
#################################################################


#########################################
# MCMC & Associated Functions
#########################################

#perspCol: plot persp with color according to z-axis 
perspCol <- function(x,y,z,color,...,xlg=TRUE,ylg=TRUE)
{
  colnames(z) <- y
  rownames(z) <- x
  
  nrz <- nrow(z)
  ncz <- ncol(z) 
  
  nb.col = length(color)
  
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  facetcol <- cut(zfacet, nb.col)
  par(xlog=xlg,ylog=ylg)
  persp(
    as.numeric(rownames(z)),
    as.numeric(colnames(z)),
    as.matrix(z),
    col=color[facetcol],
    ...
  )
}

#maybe have something that reacts to ix push x1, x2, zz to global then have plotsurface only respond to theta and phi
calcSurface <- function(ix){
  
  #set up lattice 
  if(ix==1){
    x1 <- seq(-15, 15, length=200)
    x2 <- seq(-15, 15, length=200)
    func <- rosenbrock
    zz <- matrix(apply(expand.grid(x1, x2), 1, func), nrow=length(x1))
    zz <- exp(zz)
  }
  if(ix==2){
    x1 <- seq(-2, 2, length=200)
    x2 <- seq(-2, 2, length=200)
    func <- goldpr
    zz <- matrix(apply(expand.grid(x1, x2), 1, func), nrow=length(x1))
  }
  if(ix==3){
    x1 <- seq(-512, 512, length=200)
    x2 <- seq(-512, 512, length=200)
    func <- egg
    zz <- matrix(apply(expand.grid(x1, x2), 1, func), nrow=length(x1))
  }
  if(ix==4){
    x1 <- seq(-10, 10, length=200)
    x2 <- seq(-10, 10, length=200)
    func <- levy13
    zz <- matrix(apply(expand.grid(x1, x2), 1, func), nrow=length(x1))
  }
 
  return(list(x1=x1,x2=x2,zz=zz,func=func)) 
}

plotSurface <- function(dat,theta,phi,r){
  with(dat,{
    perspCol(x = x1,y = x2,z = zz,color = viridis(60),theta=theta,phi=phi,border=NA,xlab="x",ylab="y",zlab="f(x)",ticktype="detailed")
  })
}


#plotMCMC: plot 2D surface and trace
#ix: optimization function
#seed: seed for MCMC
#x1: initial location of chain
#x2: initial location of chain
#adapt: adaptive or random walk 
#iter: number of iterations for MCMC 
#adapt_size_start: for adaptive chain 
#acceptance_rate_weight: for adaptive chain
#acceptance_window: for adaptive chain
#adapt_shape_start: for adaptive chain
runMCMC <- function(ix,seed,x1,x2,adapt,iter,adapt_size_start,acceptance_rate_weight,acceptance_window,adapt_shape_start){
 
  target = switch(ix,
    "1" = rosenbrock,
    "2" = goldpr,
    "3" = egg,
    "4" = levy13
  )
  
  if(adapt){
    mcmcOut = adaptMCMC(target = target,init_theta = c(x1,x2),covmat = diag(c(1,1)),n_iterations = iter,adapt_size_start = adapt_size_start,acceptance_rate_weight = acceptance_rate_weight,acceptance_window = acceptance_window,adapt_shape_start = adapt_shape_start,seedMH = seed)
  } else {
    mcmcOut = rwMCMC(target = target,init_theta = c(x1,x2),covmat = diag(c(1,1)),n_iterations = iter,seedMH = seed)
  }
  
  return(mcmcOut)
}


#########################################
# Test Functions for Optimization
#########################################

#Rosenbrock (banana) Function
#B: controls 'bananacity'
rosenbrock <- function(x,B=0.03) {
  -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
}


#Goldenstein-Price Function
goldpr <- function(xx){
  x1 <- xx[1]
  x2 <- xx[2]
  
  fact1a <- (x1 + x2 + 1)^2
  fact1b <- 19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2
  fact1 <- 1 + fact1a*fact1b
  
  fact2a <- (2*x1 - 3*x2)^2
  fact2b <- 18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2
  fact2 <- 30 + fact2a*fact2b
  
  y <- fact1*fact2
  return(-log(y))
}

#Eggholder Function
egg <- function(xx){

  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- -(x2+47) * sin(sqrt(abs(x2+x1/2+47)))
  term2 <- -x1 * sin(sqrt(abs(x1-(x2+47))))
  
  y <- term1 + term2
  return(-y)
}

#Levy Function Number 13
levy13 <- function(xx){
  
  x1 <- xx[1]
  x2 <- xx[2]
  
  term1 <- (sin(3*pi*x1))^2
  term2 <- (x1-1)^2 * (1+(sin(3*pi*x2))^2)
  term3 <- (x2-1)^2 * (1+(sin(2*pi*x2))^2)
  
  y <- term1 + term2 + term3
  return(-y)
}

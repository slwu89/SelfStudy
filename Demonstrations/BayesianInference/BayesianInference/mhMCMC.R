#################################################################
#
#   Metropolis Hastings Random Walk and Adaptive MCMC
#   R Shiny Helper Functions
#   Sean Wu
#   February 10, 2017
#
#################################################################


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

#plotSurface: plot function surface
plotSurface <- function(ix,theta,phi){
  
  #set up lattice 
  if(ix==1){
    x1 <- seq(-15, 15, length=100)
    x2 <- seq(-15, 15, length=100)
    func <- rosenbrock
    zz <- matrix(apply(expand.grid(x1, x2), 1, func), nrow=length(x1))
    zz <- exp(zz)
  }
  if(ix==2){
    x1 <- seq(-2, 2, length=100)
    x2 <- seq(-2, 2, length=100)
    func <- goldpr
    zz <- matrix(apply(expand.grid(x1, x2), 1, func), nrow=length(x1))
  }
  
  persp.withcol(x = x1,y = x2,z = zz,color = viridis(60),theta=theta,phi=phi,border=NA,xlab="x",ylab="y",zlab="f(x)")
}


#########################################
# Rosenbrock (banana) Function
#########################################

#B: controls 'bananacity'
rosenbrock <- function(x,B=0.03) {
  -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
}


#########################################
# Goldenstein-Price Function
#########################################

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



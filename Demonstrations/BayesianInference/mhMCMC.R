#################################################################
#
#   Metropolis Hastings Random Walk and Adaptive MCMC
#   R Shiny
#   Sean Wu
#   February 10, 2017
#
#################################################################

library(viridis)

#helper function to sample equally from color space
colGG <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

persp.withcol <- function(x,y,z,color,...,xlg=TRUE,ylg=TRUE)
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



#########################################
# Rosenbrock (banana) Function
#########################################

p.log <- function(x) {
  B <- 0.03 # controls 'bananacity'
  -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
}

#calculate function density
x1 <- seq(-15, 15, length=100)
x2 <- seq(-15, 15, length=100)
d.banana <- matrix(apply(expand.grid(x1, x2), 1, p.log), nrow=100)

#plot 3d function surface
funcCol = colorRampPalette(colors=c("#132B43","#56B1F7"))(100)
funcCol = viridis(60)
zfacet = d.banana[-1,-1] + d.banana[-1,-length(x2)] + d.banana[-length(x1),-1] + d.banana[-length(x1),-length(x2)]
perspCol = cut(zfacet,100)

# persp(x1,x2,d.banana,col=funcCol[perspCol],theta = 120,phi = 30,xlab=expression(theta[1]),ylab=expression(theta[2]),zlab=expression(f(x)))


persp.withcol(x = x1,y = x2,z = exp(d.banana),color = viridis(60),theta=120,phi=30,border=NA,xlab="x",ylab="y",zlab="f(x)")


#plot output from complex sampler
banana1 <- adaptMCMC(target=p.log,init_theta=c(10,10),covmat=diag(c(1,1)),n_iterations=1e3,
                     adapt_size_start=10,acceptance_rate_weight=0,acceptance_window=0,adapt_shape_start=20,seedMH = 50)

par(mfrow=c(1,2))


image(x1, x2, exp(d.banana), col= funcCol,xlab=expression(theta[1]),ylab=expression(theta[2])) #image of function
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.5)) #contours of function density
traceCol = viridis(n = nrow(banana1$theta_trace),option = "A",begin = 0.5) #colors of MCMC trace
for(i in 1:(nrow(banana1$theta_trace)-1)){ #plot MCMC trace
  lines(banana1$theta_trace[i:(i+1),],col=traceCol[i],type="l")
}

matplot(banana1$theta_trace,type="l",col=colGG(2),lty=1,lwd=1.25,xlab="Iterations",ylab="Parameter Trace") #
grid()

par(mfrow=c(1,1))

#plot output from simple sampler
banana2 <- rwMCMC(target = p.log,init_theta = c(10,10),covmat = diag(c(1,1)),n_iterations = 1e3,info = 1e2,seedMH = 42)

par(mfrow=c(1,2))

image(x1, x2, exp(d.banana), col= colorRampPalette(colors=c("#132B43","#56B1F7"))(100),xlab=expression(theta[1]),ylab=expression(theta[2])) #image of function
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.5)) #contours of function density
traceCol = viridis(n = nrow(banana2$theta_trace),option = "D",begin = 0.5) #colors of MCMC trace
for(i in 1:(nrow(banana2$theta_trace)-1)){ #plot MCMC trace
  lines(banana2$theta_trace[i:(i+1),],col=traceCol[i],type="l")
}

matplot(banana2$theta_trace,type="l",col=colGG(2),lty=1,lwd=1.25,xlab="Iterations",ylab="Parameter Trace") #
grid()

par(mfrow=c(1,1))


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

#calculate function density
x1 <- seq(-2, 2, length=100)
x2 <- seq(-2, 2, length=100)
d.goldpr <- matrix(apply(expand.grid(x1, x2), 1, goldpr), nrow=100)

#plot 3d function surface
# funcCol = colorRampPalette(colors=c("#132B43","#56B1F7"))(100)
# funcCol = viridis(60)
# zfacet = d.goldpr[-1,-1] + d.goldpr[-1,-length(x2)] + d.goldpr[-length(x1),-1] + d.goldpr[-length(x1),-length(x2)]
# perspCol = cut(zfacet,100)

# persp(x1,x2,d.banana,col=funcCol[perspCol],theta = 120,phi = 30,xlab=expression(theta[1]),ylab=expression(theta[2]),zlab=expression(f(x)))


persp.withcol(x1,x2,d.goldpr,color = viridis(60),theta=120,phi=30,border=NA,xlab="x",ylab="y",zlab="f(x)")

goldpr1 <- adaptMCMC(target=goldpr,init_theta=c(2,2),covmat=diag(c(1,1)),n_iterations=100,
                     adapt_size_start=100,acceptance_rate_weight=0,acceptance_window=0,adapt_shape_start=150,
                     info=1e2,seedMH = 50)

par(mfrow=c(1,2))

funCol = colorRampPalette(colors=c("#132B43","#56B1F7"))(100)

image(x1, x2, d.goldpr, col= funCol,xlab=expression(theta[1]),ylab=expression(theta[2])) #image of function
contour(x1, x2, d.goldpr, add=TRUE, col=gray(0.5)) #contours of function density
traceCol = viridis(n = nrow(goldpr1$theta_trace),option = "D",begin = 0.5) #colors of MCMC trace
for(i in 1:(nrow(goldpr1$theta_trace)-1)){ #plot MCMC trace
  lines(goldpr1$theta_trace[i:(i+1),],col=traceCol[i],type="l")
}

matplot(goldpr1$theta_trace,type="l",col=colGG(2),lty=1,lwd=1.25,xlab="Iterations",ylab="Parameter Trace") #
grid()

par(mfrow=c(1,1))
traceCol = viridis(n = nrow(goldpr1$theta_trace),option = "D",begin = 0.5) #colors of MCMC trace
par(bg="white")
goldprPmat = persp.withcol(x1,x2,d.goldpr,color = funCol,theta=120,phi=30,border=NA,xlab="x",ylab="y",zlab="f(x)",ticktype="detailed")
x = trans3d(x = goldpr1$theta_trace[,1],y = goldpr1$theta_trace[,2],z = apply(X = goldpr1$theta_trace,MARGIN = 1,FUN = goldpr),pmat = goldprPmat)
for(i in 1:(nrow(goldpr1$theta_trace)-1)){ #plot MCMC trace
  lines(x$x[i:(i+1)],x$y[i:(i+1)],col=traceCol[i],type="l")
}

points(x,pch=16,col=traceCol,cex=0.5)


#################################################################
#
#   Metropolis Hastings Random Walk and Adaptive MCMC
#   R Shiny
#   Sean Wu
#   February 10, 2017
#
#################################################################



#########################################
# Rosenbrock (banana) Function
#########################################

p.log <- function(x) {
  B <- 0.03 # controls 'bananacity'
  -x[1]^2/200 - 1/2*(x[2]+B*x[1]^2-100*B)^2
}

#plot output from complex sampler
banana1 <- adaptMCMC(target=p.log,init_theta=c(10,10),covmat=diag(c(1,1)),n_iterations=1e3,
                     adapt_size_start=10,acceptance_rate_weight=0,acceptance_window=0,adapt_shape_start=20,
                     info=1e2,seedMH = 50)

par(mfrow=c(1,2))

x1 <- seq(-15, 15, length=100)
x2 <- seq(-15, 15, length=100)
d.banana <- matrix(apply(expand.grid(x1, x2), 1, p.log), nrow=100)
image(x1, x2, exp(d.banana), col=cm.colors(60))
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(banana1$theta_trace, type='l')

matplot(banana1$theta_trace,type="l")

par(mfrow=c(1,1))

#plot output from simple sampler
banana2 <- rwMCMC(target = p.log,init_theta = c(10,10),covmat = diag(c(1,1)),n_iterations = 1e3,info = 1e2,seedMH = 42)

par(mfrow=c(1,2))

image(x1, x2, exp(d.banana), col=cm.colors(60))
contour(x1, x2, exp(d.banana), add=TRUE, col=gray(0.6))
lines(banana2$theta_trace, type='l')

matplot(banana2$theta_trace,type="l")

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

goldpr1 <- adaptMCMC(target=goldpr,init_theta=c(1,1),covmat=diag(c(1,1)),n_iterations=2e3,
                     adapt_size_start=100,acceptance_rate_weight=0,acceptance_window=0,adapt_shape_start=150,
                     info=1e2,seedMH = 50)

par(mfrow=c(1,2))

x1 <- seq(-2, 2, length=100)
x2 <- seq(-2, 2, length=100)
d.goldpr <- matrix(apply(expand.grid(x1, x2), 1, goldpr), nrow=100)
image(x1, x2, d.goldpr, col=cm.colors(100))
contour(x1, x2, d.goldpr, add=TRUE, col=gray(0.6))
lines(goldpr1$theta_trace, type='l')

matplot(goldpr1$theta_trace,type="l")

par(mfrow=c(1,1))

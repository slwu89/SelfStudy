################################################################################
#
#  Revisiting Rose:
#  Comparing the Benefits and Costs of Population-Wide and Targeted Interventions
#  R implementation of simulations in the above paper
#  Sean Wu
#  February 18, 2017
#
################################################################################


library(viridis)
library(compiler)
enableJIT(level = 3)

# create ggplot2 colors
colGG <- function(n,alpha=1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100,alpha = alpha)[1:n]
}

# logistic risk curve
logRisk <- function(x, b0 = 4, b1 = 0.3){
  1/(1 + exp(b0 - (b1 *x)))
}

# generate population
genPop <- function(n, mean = 10, sd = 2){
  rnorm(n=n,mean=mean,sd=sd)
}

# targeted treatment
treatTarget <- function(pop, threshold = 12, treat = 2){
  ix <- which(pop>threshold)
  pop[ix] <- pop[ix] - treat
  return(pop)
}

# population-wide treatment
treatPop <- function(pop, treat = 0.5){
  pop - treat
}

#perspCol: plot persp with color according to z-axis
perspCol <- function(x,y,z,color,...,xlg=TRUE,ylg=TRUE){
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


#simulation functions

#runSim: runs simulation with specified parameters
# ... : passes arguments to integrate function for calculating expected disease
runSim <- function(n,b0=4,b1=0.3,mean=10,sd=2,threshold=12,treatT=2,treatP=0.5,seed,...){

  set.seed(seed)
  simOut <- list(preTreat=NULL,targetTreat=NULL,popTreat=NULL) #output object

  #pretreatment population
  basePop <- genPop(n=n,mean=mean,sd=sd) #generate population
  baseRisk <- logRisk(x = basePop,b0 = b0,b1 = b1) #multiply risk curve by pretreatment population
  ixBase <- order(basePop) #integrating the distribution calculates the expected amount of disease
  baseSpline <- splinefun(x = basePop[ixBase],y = baseRisk[ixBase],method = "natural")
  baseDisease <- integrate(baseSpline,lower=min(basePop),upper=max(basePop),...)
  simOut$preTreat <- list(pop=basePop,risk=baseRisk,ix=ixBase,
                          spline=baseSpline,disease=baseDisease) #push to output

  #targeted treatment
  targetPop <- treatTarget(pop = basePop,threshold = threshold,treat = treatT) #apply targeted treatment
  targetRisk <- logRisk(targetPop) #multiply risk curve by targeted treatment population
  ixTarget <- order(targetPop) #integrating the distribution calculates the expected amount of disease, given targeted treatment
  targetSpline <- splinefun(x = targetPop[ixTarget],y = targetRisk[ixTarget],method = "natural")
  targetDisease <- integrate(targetSpline,lower=min(targetPop),upper=max(targetPop),...)
  simOut$targetTreat <- list(pop=targetPop,risk=targetRisk,ix=ixTarget,
                             spline=targetSpline,disease=targetDisease) #push to output

  #population-wide treatment
  popPop <- treatPop(pop = basePop,treat = treatP) #apply population-wide treatment
  popRisk <- logRisk(popPop) #multiply risk curve by population-wide treatment population
  ixPop <- order(popPop) #integrating the distribution calculates the expected amount of disease, given targeted treatment
  popSpline <- splinefun(x = popPop[ixPop],y = popRisk[ixPop],method = "natural")
  popDisease <- integrate(popSpline,lower=min(popPop),upper=max(popPop),...)
  simOut$popTreat <- list(pop=popPop,risk=popRisk,ix=ixPop,
                          spline=popSpline,disease=popDisease) #push to output

  return(simOut)
}

# plotSim2d: plot risk curves and distributions of risk factor for 3 scenarios
plotSim2d <- function(simOut){

  # plotCol <- viridis(n = 3)
  plotCol <- colGG(3,alpha = 0.75)

  #set up maximum limits for plotting
  xlim <- min(c(simOut$preTreat$pop,simOut$targetTreat$pop,simOut$popTreat$pop)) #risk factor distribution (x-axis)
  xlim[2] <- max(c(simOut$preTreat$pop,simOut$targetTreat$pop,simOut$popTreat$pop))

  # ylim <- min(c(simOut$preTreat$risk,simOut$targetTreat$risk,simOut$popTreat$risk)) #risk curve distribution (y-axis)
  # ylim[2] <- max(c(simOut$preTreat$risk,simOut$targetTreat$risk,simOut$popTreat$risk))

  par(mfrow=c(2,2))

  #plot pretreatment
  plot(simOut$preTreat$pop[simOut$preTreat$ix],simOut$preTreat$risk[simOut$preTreat$ix],
       type="l",col = "black",lwd = 1.75,xlab = paste0("Simulated Risk Factor \n Expected Disease: ",signif(simOut$preTreat$disease$value),6),
       ylab = "Risk Curve", main = "No Treatment",ylim=c(0,1),xlim=xlim) #plot the risk curve
  baseDensity <- density(simOut$preTreat$pop) #density of pretreatment population risk factor
  polygon(baseDensity$x,baseDensity$y,col=plotCol[1])
  grid()

  #plot targeted treatment
  plot(simOut$targetTreat$pop[simOut$targetTreat$ix],simOut$targetTreat$risk[simOut$targetTreat$ix],
       type="l",col = "black",lwd = 1.75,xlab = paste0("Simulated Risk Factor \n Expected Disease: ",signif(simOut$targetTreat$disease$value),6),
       ylab = "Risk Curve",main = "Targeted Treatment",ylim = c(0,1),xlim=xlim) #plot the risk curve
  targetDensity <- density(simOut$targetTreat$pop)
  polygon(targetDensity$x,targetDensity$y,col=plotCol[2])
  grid()

  #plot population-wide treatment
  plot(simOut$popTreat$pop[simOut$popTreat$ix],simOut$popTreat$risk[simOut$popTreat$ix],
       type="l",col = "black",lwd = 1.75,xlab = paste0("Simulated Risk Factor \n Expected Disease: ",signif(simOut$popTreat$disease$value),6),
       ylab = "Risk Curve", main = "Population-wide Treatment",ylim=c(0,1),xlim=xlim) #plot the risk curve
  popDensity <- density(simOut$popTreat$pop)
  polygon(popDensity$x,popDensity$y,col=plotCol[3])
  grid()


  #plot shift in distribution
  plotCol1 <- colGG(3,alpha = 0.5) #colors for densities (1: targeted, 2: population, 3: baseline/pretreatment)
  plotCol2 <- colGG(3) #colors of lines of mean
  ylim <- 0 #get y-limits (maximum of estimated densities)
  ylim[2] <- max(c(targetDensity$y,popDensity$y,baseDensity$y)) + 0.05
  plot(1, type="n", xlab="Simulated Risk Factor", ylab="Density", main="Shift in Distribution of Risk", ylim=ylim,xlim=xlim) #plot the densities
  grid()

  polygon(baseDensity$x,baseDensity$y,col=plotCol1[1]) #plot pretreatment
  baseMean <- mean(simOut$preTreat$pop)
  abline(v=baseMean,col=plotCol2[1])

  polygon(targetDensity$x,targetDensity$y,col=plotCol1[2]) #plot targeted treatment
  targetMean <- mean(simOut$targetTreat$pop)
  abline(v=targetMean,col=plotCol2[2])

  polygon(popDensity$x,popDensity$y,col=plotCol1[3]) #plot population-wide treatment
  popMean <- mean(simOut$popTreat$pop)
  abline(v=popMean,col=plotCol2[3])

  legend(x = xlim[1],y = ylim[2]-0.01,legend = c(paste0("Pretreatment Mean: ",signif(baseMean,3)),
                                                 paste0("Targeted Mean: ",signif(targetMean,3)),
                                                 paste0("Population-wide Mean: ",signif(popMean,3))),
         bty="n",pch=15,col=plotCol2)

  par(mfrow=c(1,1))

}

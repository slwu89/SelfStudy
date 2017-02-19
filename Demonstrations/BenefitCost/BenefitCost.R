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
  ifelse(pop > threshold, pop - treat, pop)
}

# population-wide treatment
treatPop <- function(pop, treat = 0.5){
  pop - treat
}


# replication of results
plotCol <- viridis(n = 3)
benefit <- 30000
targetCostTreat <- 1000
#pretreatment population
basePop <- genPop(n=1e4) #generate population

baseRisk <- logRisk(basePop) #multiply risk curve by pretreatment population

ixBase <- order(basePop) #integrating the distribution calculates the expected amount of disease
baseSpline <- splinefun(x = basePop[ixBase],y = baseRisk[ixBase],method = "natural")
baseDisease <- integrate(baseSpline,lower=min(basePop),upper=max(basePop))

plot(basePop[ixBase],baseRisk[ixBase],type="l",col=plotCol[1]) #plot the risk curve
grid()

#targeted treatment population
numTarget <- sum(ifelse(basePop > 12, T, F))
targetPop <- treatTarget(pop = basePop,threshold = 12,treat = 2) #apply targeted treatment

targetRisk <- logRisk(targetPop) #multiply risk curve by targeted treatment population

ixTarget <- order(targetPop) #integrating the distribution calculates the expected amount of disease, given targeted treatment
targetSpline <- splinefun(x = targetPop[ixTarget],y = targetRisk[ixTarget],method = "natural")
targetDisease <- integrate(targetSpline,lower=min(targetPop),upper=max(targetPop))

plot(targetPop[ixTarget],targetRisk[ixTarget],type="l",col=plotCol[2]) #plot the risk curve
grid()

#benefit (benefit of case avoided X how many cases were avoided)
targetBenefit <- benefit * (baseDisease$value - targetDisease$value)

#cost (number treated X cost of treatment for each)
targetCost <- numTarget * targetCostTreat

targetBenefit/targetCost

#population-wide treatment
popPop <- treatPop(pop = basePop) #apply population-wide treatment

popRisk <- logRisk(popPop) #multiply risk curve by population-wide treatment population

ixPop <- order(popPop) #integrating the distribution calculates the expected amount of disease, given targeted treatment
popSpline <- splinefun(x = popPop[ixPop],y = popRisk[ixPop],method = "natural")
integrate(popSpline,lower=min(popPop),upper=max(popPop))

plot(popPop[ixPop],popRisk[ixPop],type="l",col=plotCol[3]) #plot the risk curve
grid()



#simulation functions

#runSim: runs simulation with specified parameters
# ... : passes arguments to integrate function for calculating expected disease
runSim <- function(n,b0=4,b1=0.3,mean=10,sd=2,threshold=12,treatT=2,treatP=0.5,...){

  simOut <- list(preTreat=NULL,targetTreat=NULL,popTreat=NULL) #output object

  #pretreatment population
  basePop <- genPop(n=n) #generate population
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

#costBenefitSim: run cost benefit analysis with specified parameters
#n: size of population
#b0, b1: parameters of logistic curve
#mean, sd: parameters of simulated risk factor
#threshold: treatment threshold for targeted treatment
#treatT: targeted treatment reduction in risk factor
#treatP: population-wide treatment reduction in risk factor
#benefit: monetary benefit of avoiding one case of disease
#targetCostTreat: cost for each targeted treatment
# costBenefitSim <- function(n,b0=4,b1=0.3,mean=10,sd=2,threshold=12,treatT=2,treatP=0.5,
#                            benefit=30000,targetCostTreat=1000){
#
#   simOut <- list(preTreat=NULL,targetTreat=NULL,popTreat=NULL) #output object
#
#   #pretreatment population
#   basePop <- genPop(n=n) #generate population
#   baseRisk <- logRisk(x = basePop,b0 = b0,b1 = b1) #multiply risk curve by pretreatment population
#   ixBase <- order(basePop) #integrating the distribution calculates the expected amount of disease
#   baseSpline <- splinefun(x = basePop[ixBase],y = baseRisk[ixBase],method = "natural")
#   baseDisease <- integrate(baseSpline,lower=min(basePop),upper=max(basePop))
#   simOut$preTreat <- list(pop=basePop,risk=baseRisk,ix=ixBase,
#                           spline=baseSpline,disease=baseDisease) #push to output
#
#   #targeted treatment
#   numTarget <- sum(ifelse(basePop > threshold,T,F)) #how many patients must be treated
#   targetPop <- treatTarget(pop = basePop,threshold = threshold,treat = treatT) #apply targeted treatment
#   targetRisk <- logRisk(targetPop) #multiply risk curve by targeted treatment population
#   ixTarget <- order(targetPop) #integrating the distribution calculates the expected amount of disease, given targeted treatment
#   targetSpline <- splinefun(x = targetPop[ixTarget],y = targetRisk[ixTarget],method = "natural")
#   targetDisease <- integrate(targetSpline,lower=min(targetPop),upper=max(targetPop))
#   simOut$targetTreat <- list(pop=targetPop,risk=targetRisk,ix=ixTarget,
#                              spline=targetSpline,disease=targetDisease) #push to output
#
#   #benefit (benefit of case avoided X how many cases were avoided)
#   targetBenefit <- benefit * (((baseDisease$value/100) * n) - ((targetDisease$value/100) * n))
#   # targetBenefit <- benefit * (baseDisease$value - targetDisease$value)
#   #cost (number treated X cost of treatment for each)
#   targetCost <- numTarget * targetCostTreat
#
#   #population-wide treatment
#   popPop <- treatPop(pop = basePop,treat = treatP) #apply population-wide treatment
#   popRisk <- logRisk(popPop) #multiply risk curve by population-wide treatment population
#   ixPop <- order(popPop) #integrating the distribution calculates the expected amount of disease, given targeted treatment
#   popSpline <- splinefun(x = popPop[ixPop],y = popRisk[ixPop],method = "natural")
#   popDisease <- integrate(popSpline,lower=min(popPop),upper=max(popPop))
#   simOut$popTreat <- list(pop=popPop,risk=popRisk,ix=ixPop,
#                           spline=popSpline,disease=popDisease) #push to output
#
#   return(simOut)
# }

#given a baseline population calc the benefit/cost ratio
#basePop is "preTreat" element from simOut
calcTargetRatio <- function(threshold,targetCostTreat,basePop,benefit){

  numTarget <- sum(ifelse(basePop$pop > threshold,T,F)) #how many patients must be treated
  targetPop <- treatTarget(pop = basePop$pop,threshold = threshold,treat = 2) #apply targeted treatment
  targetRisk <- logRisk(targetPop) #multiply risk curve by targeted treatment population
  ixTarget <- order(targetPop) #integrating the distribution calculates the expected amount of disease, given targeted treatment
  targetSpline <- splinefun(x = targetPop[ixTarget],y = targetRisk[ixTarget],method = "natural")
  targetDisease <- integrate(targetSpline,lower=min(targetPop),upper=max(targetPop))

  #benefit (benefit of case avoided X how many cases were avoided)
  targetBenefit <- benefit * (basePop$disease$value - targetDisease$value)

  #cost (number treated X cost of treatment for each)
  targetCost <- numTarget * targetCostTreat

  return(targetBenefit/targetCost)
}

calcTargetRatio <- Vectorize(calcTargetRatio,vectorize.args = c("threshold","targetCostTreat"))

thresholdVector <- seq(from=8,to=16,by=0.05)
costVector <- seq(from=600,to=1400,by=10)

targetSurface <- outer(X = thresholdVector,Y = costVector,FUN = calcTargetRatio, basePop=basePop, benefit=30000)

perspCol(x = thresholdVector,y = costVector,z = targetSurface,color = viridis(60),theta=330,phi=30,border=NA,xlab="x1",ylab="x2",zlab="f(x)",ticktype="detailed")

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




#plotSim2d: plot risk curves and distributions of risk factor for 3 scenarios
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
       ylab = "Risk Curve", main = "No Treatment",ylim=c(0,1)) #plot the risk curve
  baseDensity <- density(simOut$preTreat$pop) #density of pretreatment population risk factor
  polygon(baseDensity$x,baseDensity$y,col=plotCol[1])
  grid()

  #plot targeted treatment
  plot(simOut$targetTreat$pop[simOut$targetTreat$ix],simOut$targetTreat$risk[simOut$targetTreat$ix],
       type="l",col = "black",lwd = 1.75,xlab = paste0("Simulated Risk Factor \n Expected Disease: ",signif(simOut$targetTreat$disease$value),6),
       ylab = "Risk Curve",main = "Targeted Treatment",ylim = c(0,1)) #plot the risk curve
  targetDensity <- density(simOut$targetTreat$pop)
  polygon(targetDensity$x,targetDensity$y,col=plotCol[2])
  grid()

  #plot population-wide treatment
  plot(simOut$popTreat$pop[simOut$popTreat$ix],simOut$popTreat$risk[simOut$popTreat$ix],
       type="l",col = "black",lwd = 1.75,xlab = paste0("Simulated Risk Factor \n Expected Disease: ",signif(simOut$popTreat$disease$value),6),
       ylab = "Risk Curve", main = "Population-wide Treatment",ylim=c(0,1)) #plot the risk curve
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

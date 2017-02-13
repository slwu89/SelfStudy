#################################################################
#
#   Metropolis Hastings Random Walk and Adaptive MCMC
#   R Shiny Server
#   Sean Wu
#   February 10, 2017
#
#################################################################


###########################################
# load libraries and source cpp
###########################################

library(shiny)

#Rcpp libraries
library(Rcpp)
library(RcppArmadillo)
library(RcppProgress)

#graphical libraries
library(viridis)

#compile Metropolis-Hastings MCMC code
Rcpp::sourceCpp("mhMCMC.cpp")

#source R functions
source("mhMCMC.R")


###########################################
# code for server
###########################################

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  ##########################################################
  #  Metropolis-Hastings MCMC
  ########################################################## 
  
  valuesMH <- reactiveValues() 
  
  #theta and phi adjust the angle of view for persp plots
  observe({
    valuesMH$theta <- input$mhTheta
    valuesMH$phi <- input$mhPhi
  })
  
  #only evaluate function surface when user changes input function for test optimization
  observe({
    valuesMH$funcSurface <- calcSurface(input$mhFunction)
  })
  
  #plot function surface
  output$mhFunction <- renderPlot({
    plotSurface(dat = valuesMH$funcSurface,theta = valuesMH$theta,phi = valuesMH$phi)
  },width = 600,height = 600)
  
  #run MCMC
  observeEvent(input$mhRun,{
    #run MCMC
    valuesMH$mcmc <- runMCMC(ix = input$mhFunction,seed = input$mhSeed,x1 = input$mhX1,x2 = input$mhX2,mcmcType = input$mhMCMC,diag = input$mhDiag,
                    iter = input$mhIter,adapt_size_start = input$mhAdaptSizeStart,acceptance_rate_weight = input$mhAcceptanceRateWeight,acceptance_window = input$mhAcceptanceWindow,adapt_shape_start = input$mhAdaptShapeStart)
    #2d plots
    output$mhMCMC2d <- renderPlot({
      plotMCMC(funcDat = valuesMH$funcSurface,mcmcDat = valuesMH$mcmc)
    },width = 1000,height = 600)
    
    #plot MCMC perspective
    output$mhMCMC3d <- renderPlot({
      plotMCMCpersp(funcDat = valuesMH$funcSurface,mcmcDat = valuesMH$mcmc,theta = valuesMH$theta,phi = valuesMH$phi)
    },width = 800,height = 800)
  })
  
})

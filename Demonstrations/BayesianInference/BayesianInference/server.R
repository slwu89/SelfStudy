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
  
  observe({
    
    valuesMH$mhFunction <- input$mhFunction
    valuesMH$theta <- input$mhTheta
    valuesMH$phi <- input$mhPhi
    
  })
  
  #plot function surface
  output$mhFunction <- renderPlot({
    plotSurface(ix = valuesMH$mhFunction,theta = valuesMH$theta,phi = valuesMH$phi)
  },width = 600,height = 600)
  
})

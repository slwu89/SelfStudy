#################################################################
#
#   Metropolis Hastings Random Walk and Adaptive MCMC
#   R Shiny UI
#   Sean Wu
#   February 10, 2017
#
#################################################################


library(shiny)
library(shinythemes)
# library(shinyjs)
# library(shinydashboard)
library(markdown)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(theme=shinytheme("flatly"),withMathJax(),
            
            navbarPage("Bayesian Inference",
                       
                       navbarMenu("Markov Chain Monte Carlo",
                                  
                                  tabPanel("Metropolis-Hastings",
                                           
                                           fluidRow(
                                             
                                             column(width = 4, offset = 0, style="padding:0px;",
                                                    includeMarkdown("mhMCMC.md"),
                                                    hr(),
                                                    h3("MCMC Parameters"),
                                                    selectInput(inputId="mhFunction",label="Optimization Function", 
                                                                choices = list("Rosenbrock" = 1, "Goldenstein-Price" = 2), 
                                                                selected = 1),
                                                    sliderInput(inputId="mhTheta",
                                                                label=HTML("&Theta; (Azimuth)"),
                                                                min = 0, max = 360, value = 120, step = 1),
                                                    sliderInput(inputId="mhPhi",
                                                                label=HTML("&Phi; (Colatitude)"),
                                                                min = 0, max = 360, value = 30, step = 1)
                                                    ),
                                             column(width = 8, offset = 0, style="padding:10px;",
                                                    plotOutput("mhFunction",width = "auto",height = "auto")
                                                    )
                                             
                                                  ) #end Metropolis-Hastings fluidRow
                                           
                                           ) #end Metropolis-Hastings tabPanel
                                  
                                  ) #end Markov Chain Monte Carlo navbarMenu
                       
                       ) #end Bayesian Inference navbarPage
            
  ) #end fluidPage
) #end shinyUI

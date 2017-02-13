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
library(shinydashboard)
library(markdown)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(theme=shinytheme("flatly"),withMathJax(),
            
            navbarPage("Bayesian Inference",
                       
                       navbarMenu("Markov Chain Monte Carlo",
                                  
                                  tabPanel("Metropolis-Hastings",
                                           
                                           fluidRow(
                                             
                                             column(width = 2, offset = 0, style="padding:0px;",
                                                    h3("MCMC Parameters"),
                                                    hr(),
                                                    selectInput(inputId="mhMCMC",label="Algorithm", 
                                                                choices = list("Random Walk Metropolis-Hastings" = 1, "Robust Adaptive Metropolis" = 2), 
                                                                selected = 1),
                                                    numericInput(inputId="mhIter",label="Iterations",
                                                                 value=1e3,min=1,max=1e4
                                                                 ),
                                                    numericInput(inputId="mhSeed",label="Seed",
                                                                 value=42,min=1
                                                                 ),
                                                    numericInput(inputId="mhAdaptSizeStart",label="Adapt Size Start",
                                                                 value=0,min=0
                                                                 ),
                                                    numericInput(inputId="mhAdaptShapeStart",label="Adapt Shape Start",
                                                                 value=0,min=0),
                                                    numericInput(inputId="mhAcceptanceRateWeight",label="Acceptance Rate Weight",
                                                                 value=0,min=0,max=1),
                                                    numericInput(inputId="mhAcceptanceWindow",label="Acceptance Window",
                                                                 value=0,min=0,max=1)
                                                    ),
                                              column(width = 2, offset = 0, style="padding:10px;",
                                                     sliderInput(inputId="mhTheta",
                                                                 label=HTML("&Theta; (Azimuth)"),
                                                                 min = 0, max = 360, value = 120, step = 1),
                                                     sliderInput(inputId="mhPhi",
                                                                 label=HTML("&Phi; (Colatitude)"),
                                                                 min = 0, max = 360, value = 30, step = 1),
                                                     numericInput(inputId="mhDiag",label="Sigma Diagonal",
                                                                  value=1,min=1,max=500),
                                                     numericInput(inputId="mhX1",label=HTML("X<sub>1</sub>"),
                                                                  value=10),
                                                     numericInput(inputId="mhX2",label=HTML("X<sub>2</sub>"),
                                                                  value=10),
                                                     actionButton(inputId="mhRun",label="Run MCMC",style="background-color:#428bca;")
                                                     ),
                                             tabBox(id="tabBoxMCMC",
                                                    tabPanel("Function Surface",
                                                        fluidPage(
                                                          fluidRow(
                                                            column(width = 4, offset = 0, style="padding:10px;",
                                                                   selectInput(inputId="mhFunction",label="Optimization Function", 
                                                                               choices = list("Camel" = 1, "Goldenstein-Price" = 2, "Eggholder" = 3, "Levy N.13" = 4), 
                                                                               selected = 1),
                                                                    includeMarkdown("mhFunctions.md")
                                                                   ),
                                                            column(width = 8, offset = 0, style="padding:0px;",
                                                                   plotOutput("mhFunction",width = "auto",height = "auto")
                                                            )
                                                          ),
                                                          fluidRow(
                                                            column(width = 12,offset = 0, style="padding:0px;",
                                                                   includeMarkdown("mhMCMC.md")
                                                            )
                                                          )
                                                        )       
                                                      ), #end Function Surface tabPanel
                                                      tabPanel("MCMC",
                                                        fluidPage(
                                                          fluidRow(
                                                            column(width = 8, offset = 0, style="padding:0px;",
                                                                   plotOutput("mhMCMC2d",width = "auto",height = "auto")
                                                            )
                                                          ),
                                                          fluidRow(
                                                            column(width = 12,offset = 0, style="padding:0px;",
                                                                   includeMarkdown("mhMCMC.md")
                                                            )
                                                          )
                                                        )
                                                      ), #end MCMC tabPanel
                                                    tabPanel("MCMC Surface",
                                                             fluidPage(
                                                               fluidRow(
                                                                 column(width = 8, offset = 0, style="padding:0px;",
                                                                        plotOutput("mhMCMC3d",width = "auto",height = "auto")
                                                                 )
                                                               ),
                                                               fluidRow(
                                                                 column(width = 12,offset = 0, style="padding:0px;",
                                                                        includeMarkdown("mhMCMC.md")
                                                                 )
                                                               )
                                                             )
                                                      ) #end MCMC surface tabPanel
                                                    ) #end tabBoxMCMC tabBox
                                                  ) #end Metropolis-Hastings fluidRow
                                           
                                           ) #end Metropolis-Hastings tabPanel
                                  
                                  ) #end Markov Chain Monte Carlo navbarMenu
                       
                       ) #end Bayesian Inference navbarPage
            
  ) #end fluidPage
) #end shinyUI

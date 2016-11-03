#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI 
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Optimal Decision Strategy through Dynamic Programming"),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("t_max",
                  "Maximum Time",
                  min = 1,
                  max = 500,
                  value = 30),
      sliderInput("npatch",
                  "Number of Patches",
                  min = 2,
                  max = 10,
                  value = 3),
      textInput("psurvive",
                "Survival Probabilities",
                value="0.99,0.95,0.98"),
      textInput("pfood",
                "Food Acquisition Probabilities",
                value="0.2,0.5,0.0"),
      textInput("cost",
                "Travel Cost (in Biomass)",
                value="1,1,1"),
      textInput("feedgain",
                "Biomass Gain from Feeding",
                value="2,4,0"),
      textInput("repr",
                "Fitness Gain from Reproduction",
                value="0,0,4"),
      selectInput("x_crit",
                  "Minimum Biomass for Survival",
                  choices=c(0,1,2,3,4,5)),
      selectInput("x_max",
                  "Maximum Biomass",
                  choices=c(20,30,40,50,60,70)),
      selectInput("x_rep",
                  "Biomass to Reproduce",
                  choices=c(3,4,5,6,7,8)),
      selectInput("fend",
                  "Final Maximum Fitness",
                  choices=c(50,60,70,80,90,100)),
      actionButton(inputId="run",label="Evaluate")
    ),
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("optimalPatch")
    )
  )
))

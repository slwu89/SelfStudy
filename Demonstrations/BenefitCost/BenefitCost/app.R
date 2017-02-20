#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(shinyjs)
source("functions.R")

# Define UI for application that draws a histogram
ui <- fluidPage(theme=shinytheme("flatly"),shinyjs::useShinyjs(),

    navbarPage("Revisiting Rose: Comparing the Benefits and Costs of Population-Wide and Targeted Interventions"),

   fluidRow(
     column(width=3,offset = 0, style="padding:20px;",
            numericInput(inputId="n",label="Population Size",
                         value=1e3,min=10,max=1e4
            ),
            numericInput(inputId="seed",label="Seed",
                         value=42,min=1,max=2e16
            ),
            sliderInput(inputId="mean",
                        label="Mean of Risk Factor",
                        min = 8, max = 16, value = 10, step = 1),
            sliderInput(inputId="sd",label="SD of Risk Factor",
                        min = 0.1, max = 4, value = 2, step = 0.1),
            sliderInput(inputId="threshold",label="Targeted Treatment Threshold",
                        min = 6, max = 16, value = 12, step = 0.1),
            sliderInput(inputId="treatT",label="Targeted Treatment",
                        min = 0.1, max = 4, value = 2, step = 0.1),
            sliderInput(inputId="treatP",label="Population-wide Treatment",
                        min = 0.1, max = 4, value = 0.5, step = 0.1),
            actionButton(inputId="resetP",label="Reset Parameters",
                        style='text-align:center; font-size:75%; width:55%; height:80%;
                        background-color:#428bca; position:center; margin:5 auto;
                         display:block')
            ),

     column(width=9,offset = 0, style="padding:10px;",
            plotOutput("plot2d",width = "auto",height = "auto")
            )
   )

)

# Define server logic required to draw a histogram
server <- function(input, output) {

  parms <- reactiveValues()

  #theta and phi adjust the angle of view for persp plots
  observe({
    parms$n <- input$n
    parms$mean <- input$mean
    parms$sd <- input$sd
    parms$threshold <- input$threshold
    parms$treatT <- input$treatT
    parms$treatP <- input$treatP
    parms$seed <- input$seed
  })

  observeEvent(input$resetP,{
    shinyjs::reset(id="n")
    shinyjs::reset(id="seed")
    shinyjs::reset(id="mean")
    shinyjs::reset(id="sd")
    shinyjs::reset(id="threshold")
    shinyjs::reset(id="treatT")
    shinyjs::reset(id="treatP")
  })

  output$plot2d <- renderPlot({
    plotSim2d(runSim(n = parms$n,mean = parms$mean,sd = parms$sd,threshold = parms$threshold,treatT = parms$treatT,treatP = parms$treatP,seed=parms$seed))
  },height = 950,width = 950)
}

# Run the application
shinyApp(ui = ui, server = server)


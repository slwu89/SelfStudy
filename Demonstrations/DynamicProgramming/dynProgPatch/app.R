###############################################################################
######Demonstrations of Dynamic Programming for Optimal Decision Strategy######
######Sean Wu 11/2/2016########################################################
###############################################################################

library(RColorBrewer)
library(shiny)


##################################
###Simple Patch-Selection Model###
##################################

#based on code in "A Practical Guide to Ecological Modeling"

#fitness:
#x: biomass
#t: current time
#x_max: maximal biomass allowed
#x_crit: minimum biomass to survive
#f: matrix of optimal fitness values
fitness <- function(x,t,x_max,x_crit,f){
  xx <- pmin(x,x_max)
  xx <- pmax(xx,x_crit)
  out <- f[t,(xx-x_crit+1)]
  return(out)
}

#dynProgPatch: run the dynamic programming patch model
#npatch: number of patches
#psurvive: numeric vector of survival probability upon travel to that patch
#pfood: numeric vector of probability to find food at that patch
#cost: numeric vector of cost in biomass to travel to that patch
#feedgain: numeric vector of biomass gain from food at that patch
#repr: numeric vector of maximal reproduction gain at that patch
#x_crit: minimum biomass to survive
#x_max: maximal biomass allowed
#x_rep: minimum biomass to reproduce
#t_max: number of time steps
#fend: final maximum possible fitness at t_max
dynProgPatch <- function(npatch=3,psurvive=c(0.99,0.95,0.98),pfood=c(0.2,0.5,0.0),cost=c(1,1,1),
                         feedgain=c(2,4,0),repr=c(0,0,4),x_crit=0,x_max=30,x_rep=4,t_max=20,fend=60,
                         updateProgress=NULL){
  
  x_class <- x_crit:x_max #biomass classes
  nmass <- length(x_class) #number of mass classes
  times <- 1:(t_max-1) #number of time steps
  kx <- 0.25*x_max #half-saturation coefficient for Michaelis-Menten response

  #allocate containers
  f <- matrix(0,nrow=t_max,ncol=nmass) #optimal fitness values
  bestpatch <- matrix(0,nrow=t_max-1,ncol=nmass-1) #best patch
  V <- vector(mode="numeric",length=npatch) #current fitness
  
  f[t_max,] <- fend*(x_class-x_crit) / (x_class-x_crit+kx) #final fitness evaluated as Michaelis-Menten response
  
  #run dynamic programming optimzation loop
  for(t in rev(times)){ #iterate backwards in time
    
    if(is.function(updateProgress)){
      updateProgress()
    }
    
    for(x in x_class[-1]){ #calculate for all classes except x_crit
      
      dfit <- pmax(0,pmin(x-x_rep,repr)) #possible immediate fitness gain for class x at each patch
      #expected gain in fitness:
      #probability of survival X (probability to find food X biomass at next time if found food + probability to not find food X biomass at next time if did not find food)
      #prob of survival X ((prob of food X fitness) if fed + ((1-prob) of food X fitness if not fed))
      expectgain <- psurvive * (pfood * fitness(x=(x-cost+feedgain-dfit),t=t+1,x_max=x_max,x_crit=x_crit,f=f) +
                                  (1-pfood) * fitness(x=(x-cost-dfit),t=t+1,x_max=x_max,x_crit=x_crit,f=f)) 
      
      V <- dfit + expectgain #current fitness
      V[expectgain==0] <- 0 #kill animals below critical value
      f[t,(x-x_crit+1)] <- max(V) #optimal fitness
      bestpatch[t,(x-x_crit)] <- which.max(V) #best patch
      
    } #next biomass class X
    
  }#next time t
  
  image(x=times,y=x_class[-1],z=bestpatch,zlim=c(0,npatch),main="Optimal Patch Selection",
        col=c("Black",brewer.pal(npatch,"Spectral")),xlab="Time",ylab="Biomass",
        cex.main=2,cex.lab=1.5)
  contour(x=1:t_max,y=x_class,z=f,add=TRUE,drawlabels=TRUE,labcex=2)
  legend(x="topleft",fill=c("Black",brewer.pal(npatch,"Spectral")),legend=c("Dead",paste0(1:npatch)),bty="n")
  legend(x="topright",legend="Fitness",lty=1,bty="n")
  # return(list(times=times,f=f,bestpatch=bestpatch,x_class=x_class,npatch=npatch,t_max=t_max))
}

#plotPatch: plot the output
plotPatch <- function(times,x_class,bestpatch,npatch,t_max){
  image(x=times,y=x_class[-1],z=bestpatch,zlim=c(0,npatch),main="Optimal Patch Selection",
        col=c("Black",brewer.pal(npatch,"Spectral")),xlab="Time",ylab="Biomass",
        cex.main=2,cex.lab=1.5)
  contour(x=1:t_max,y=x_class,z=f,add=TRUE,drawlabels=TRUE,labcex=2)
  legend(x="topleft",fill=c("Black",brewer.pal(npatch,"Spectral")),legend=c("Dead",paste0(1:npatch)),bty="n")
  legend(x="topright",legend="Fitness",lty=1,bty="n")
}



###################
###Shiny Web App###
###################

# Define UI
ui <- fluidPage(

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
)

# Define server
server <- function(input, output) {

  # observeEvent(input$run,{
  #   output$optimalPatch <- renderPlot({
  # 
  #     #progress bar
  #     progress <- shiny::Progress$new()
  #     progress$set(message = "Computing Optimal Patch Strategy", value = 0)
  #     on.exit(progress$close())
  # 
  #     #redefine input as numerics
  #     psurviveIn <- as.numeric(unlist(strsplit(input$psurvive,",")))
  #     pfoodIn <- as.numeric(unlist(strsplit(input$pfood,",")))
  #     costIn <- as.numeric(unlist(strsplit(input$cost,",")))
  #     feedgainIn <- as.numeric(unlist(strsplit(input$feedgain,",")))
  #     reprIn <- as.numeric(unlist(strsplit(input$repr,",")))
  #     x_crit <- as.numeric(input$x_crit)
  #     x_max <- as.numeric(input$x_max)
  #     x_rep <- as.numeric(input$x_rep)
  #     fend <- as.numeric(input$fend)
  # 
  #     updateProgress <- function() {
  #       # progress$set(value = value, detail = detail)
  #       progress$inc(1/(input$t_max-1))
  #     }
  # 
  #     dynProgPatch(npatch = input$npatch,psurvive = psurviveIn,pfood = pfoodIn,cost = costIn,feedgain = feedgainIn,
  #                  repr = reprIn,x_crit = x_crit,x_max = x_max,x_rep = x_rep,t_max = input$t_max,fend = fend,
  #                  updateProgress=updateProgress)
  #   })
  # })
  

    output$optimalPatch <- renderPlot({
      
      input$run
      
      if (input$run == 0){
        return()
      }
        
      
      #progress bar
      progress <- shiny::Progress$new()
      progress$set(message = "Computing Optimal Patch Strategy", value = 0)
      on.exit(progress$close())
      
      #redefine input as numerics
      psurviveIn <- as.numeric(unlist(strsplit(input$psurvive,",")))
      pfoodIn <- as.numeric(unlist(strsplit(input$pfood,",")))
      costIn <- as.numeric(unlist(strsplit(input$cost,",")))
      feedgainIn <- as.numeric(unlist(strsplit(input$feedgain,",")))
      reprIn <- as.numeric(unlist(strsplit(input$repr,",")))
      x_crit <- as.numeric(input$x_crit)
      x_max <- as.numeric(input$x_max)
      x_rep <- as.numeric(input$x_rep)
      fend <- as.numeric(input$fend)
      
      updateProgress <- function() {
        # progress$set(value = value, detail = detail)
        progress$inc(1/(input$t_max-1))
      }
      
      isolate(dynProgPatch(npatch = input$npatch,psurvive = psurviveIn,pfood = pfoodIn,cost = costIn,feedgain = feedgainIn,
                   repr = reprIn,x_crit = x_crit,x_max = x_max,x_rep = x_rep,t_max = input$t_max,fend = fend,
                   updateProgress=updateProgress))
    })
  

}

# Run the application
shinyApp(ui = ui, server = server)
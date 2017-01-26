#################################################################
#
#   Mushroom Hunt
#   Soc 273A
#   Sean Wu
#   January 26, 2017
#
#################################################################

using Distributions
using Plots


##########################################
# Distribute Points
##########################################

#pointsPoisson: generate homogenous Poisson point process
#n: number of points
#xLim: x limits
#yLim: y limits
function pointsPoisson(n,xLim=(0,1),yLim=(0,1))
  xDiff = xLim[2]-xLim[1]
  yDiff = yLim[2]-yLim[1]
  xy = zeros(Float64, n, 2)
  for i in 1:n
    xy[i,:] = [xLim[1]+rand()*xDiff, yLim[1]+rand()*yDiff]
  end
  return xy
end

#pointsClustered: generate a Neyman-Scott process
#n: number of parents
#c: number of children
function pointsClustered(n,c,sd=0.01,xLim=(0,1),yLim=(0,1))
  #generate parents
  xDiff = xLim[2]-xLim[1]
  yDiff = yLim[2]-yLim[1]
  xy = zeros(Float64, n+c, 2)
  for i in 1:n
    xy[i,:] = [xLim[1]+rand()*xDiff, yLim[1]+rand()*yDiff]
  end
  #generate children
  for i in 1:c
    cPar = sample(1:n)
    xy[i+n,:] = [xy[cPar,1]+rand(Normal(0,sd)) ,xy[cPar,2]+rand(Normal(0,sd))]
  end
  return xy
end

#test function
poisXY = pointsPoisson(20)
scatter(poisXY[:,1],poisXY[:,2],color=:steelblue,label="Patches")

clusXY = pointsClustered(5,15,0.05)
scatter(clusXY[:,1],clusXY[:,2],color=:steelblue,label="Patches")


##########################################
# Define Agents
##########################################

#make a turtle
function makeTurtle(ix,x,y)
  return Dict(:ix => ix, :x => x, :y => y, :pointHist => [])
end

#make some turtles
turtles = []
for i in 1:50
  push!(turtles,makeTurtle(i,rand(),rand()))
end

#plot the turtles
scatter!([turtle[:x] for turtle in agents],[turtle[:y] for turtle in agents],color=:red,label="Turtles")

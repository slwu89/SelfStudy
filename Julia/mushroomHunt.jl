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
scatter(poisXY[:,1],poisXY[:,2],color=:steelblue,label="Patches" bg=RGB(.2,.2,.2))

clusXY = pointsClustered(5,15,0.025)
scatter(clusXY[:,1],clusXY[:,2],color=:steelblue,label="Patches",m=(5),bg=RGB(.4,.4,.4))

#xy coords of mushrooms
mushrooms = clusXY = pointsPoisson(25)
#world limits as tuples
xLim = (0,1)
yLim = (0,1)

##########################################
# Define Agents
##########################################

#make a turtle
function makeTurtle(ix,x,y,xyDir)
  return Dict(:ix => ix, :xy => [x,y], :xyV => Float64,
  :timeSinceLastFound => 999, :pointHist => [])
end

#make some turtles
turtles = []
for i in 1:5
  push!(turtles,makeTurtle(i,0.5,0.5,rand()*360))
end

#get turtle xy
function getTurtleXY()
  global turtles
  xyArray = Array{Float64}(length(turtles),2)
  for i in 1:length(turtles)
    xyArray[i,:] = [turtles[i][:xy][1],turtles[i][:xy][2]]
  end
  return xyArray
end

#plot the turtles
turtleXY = getTurtleXY()
scatter!(turtleXY[:,1],turtleXY[:,2],color=:red,label="Turtles",m=(5,:star7))


##########################################
# Search Code
##########################################

#search routine
speed = 0.001 #speed of turtles
senseDist = 0.1 #sensing distance

#create pairwise Euclidean distance matrix between turtles and mushrooms
function turtleMushroomDistanceMatrix()
  global turtles, mushrooms
  #distance matrix n_turtles X n_mushrooms
  distMat = Array{Float64}(length(turtles),size(mushrooms)[1])
  for i in 1:length(turtles)
    for j in 1:size(mushrooms)[1]
      distMat[i,j] = sqrt((turtles[i][:xy][1] - mushrooms[j,1])^2 + (turtles[i][:xy][2] - mushrooms[j,2])^2)
    end
  end
  return distMat
end

#search routine for turtle ix
#ix: index of turtle
#tmDist: turtle-mushroom distance matrix
#if multiple mushrooms are within its sensing distance it will go for the closest one
function searchMushrooms(ix,tmDist)
  boolSense = tmDist[ix,:] .< senseDist
  #if no mushrooms in range
  if any(!boolSense)

  elseif sum(boolSense) > 1 #if multiple mushrooms in range
    target = findmin(tmDist[ix,:])[2]
  else
    target = find(boolSense)
  end
end

#move turtle ix
function moveTurtle(ix)
  global turtles

  turtles[ix][:xy] = turtles[ix][:xy] + veclocity

end

tmDist = turtleMushroomDistanceMatrix()

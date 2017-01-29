#################################################################
#
#   Wright-Fisher model of mutation, selection and random genetic drift
#   R Shiny app with Julia source code
#   Sean Wu
#   January 28, 2017
#
#################################################################

library(XRJulia)

#named list of parameters to send to Julia
paramsList <- list(
  pop_size = 50,
  seq_length = 100,
  generations = 500,
  mutation_rate = 0.0001, #per gen per individual per site
  fitness_effect = 1.1, #fitness effect if a functional mutation occurs
  fitness_chance = 0.1, #chance that a mutation has a fitness effect
  alleleFreq = c(0.5,0.5) #allele frequencies for two letter alphabet
)

# paramsVec = c(50,100,500,0.0001,1.1,0.1)

# js = juliaSend(paramsList) #currently deprecated in Julia v0.5.0; wait for XRJulia update

# jlParams = juliaEval("Dict(:popSize=>%s, :seqLen=>%s, :generations=>%s, :mutationRate=>%s, :fitnessEffect=>%s, :fitnessChance=>%s)",
#                  paramsVec[1],paramsVec[2],paramsVec[3],paramsVec[4],paramsVec[5],paramsVec[6])

jlParams = juliaEval("Dict(:popSize=>%s, :seqLen=>%s, :generations=>%s, :mutationRate=>%s, :fitnessEffect=>%s, :fitnessChance=>%s, :alleleFreq=>%s)",
                     paramsList[[1]],paramsList[[2]],paramsList[[3]],paramsList[[4]],paramsList[[5]],paramsList[[6]],paramsList[[7]])


juliaGet(jlParams)

juliaSource()
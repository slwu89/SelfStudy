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

jlParams = juliaEval("Dict(:popSize=>%s, :seqLen=>%s, :generations=>%s, :mutationRate=>%s, :fitnessEffect=>%s, :fitnessChance=>%s, :alleleFreq=>%s)",
                     paramsList[[1]],paramsList[[2]],paramsList[[3]],paramsList[[4]],paramsList[[5]],paramsList[[6]],paramsList[[7]])


juliaGet(jlParams)

generateHaplotype = juliaEval("
          function generate_haplotype(pars)
                               hap = wsample(['A','T'],pars[:alleleFreq])
                               [hap = hap * wsample(['A','T'],pars[:alleleFreq]) for i in 2:pars[:seqLen]]
                               return hap
                               end
                               ")



# jlMod = juliaSource("/Users/slwu89/Desktop/git/SelfStudy/Demonstrations/MutationDriftSelection/MutationDriftSelection.jl")
# 
# juliaAddToPath("/Users/slwu89/Desktop/git/SelfStudy/Demonstrations/MutationDriftSelection/")
# juliaImport("MutationDriftSelection.jl")

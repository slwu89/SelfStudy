#################################################################
#
#   Wright-Fisher model of mutation, selection and random genetic drift
#   Julia port of Python 2.5 code (SISMID 2016)
#   Sean Wu
#   January 11, 2017
#
#################################################################

using Distributions

#########################################################
#  Make population dynamic model
##########################################################

#Basic Parameters

pop_size = 100
seq_length = 10
alphabet = ["A","T"]
base_haplotype = "AAAAAAAAAA"
fitness_effect = 1.1 #fitness effect if a functional mutation  occurs
fitness_chance = 0.1 #chance that a mutation has a fitness effect


#Population of haplotypes maps to counts and fitnesses:
#Store this as a lightweight dictionary that maps a string to a count

pop = Dict("AAAAAAAAAA" => 40)
pop["AAATAAAAAA"] = 30
pop["AATTTAAAAA"] = 30

#Map haplotype string to fitness float
fitness = Dict("AAAAAAAAAA" => 1.0)
fitness["AAATAAAAAA"] = 1.05
fitness["AATTTAAAAA"] = 1.10

#Add mutation

mutation_rate = 0.005 # per gen per individual per site

function get_mutation_count()
    mean = mutation_rate * pop_size * seq_length
    return rand(Poisson(mean))
end

function get_random_haplotype()
    haplotypes = collect(keys(pop))
    frequencies = [x/pop_size for x in values(pop)]
    total = sum(frequencies)
    frequencies = [x/total for x in frequencies]
    return wsample(haplotypes,frequencies)
end

function get_mutant(haplotype)
  site = rand(1:seq_length)
  possible_mutations = alphabet
  possible_mutations = filter(x -> x != haplotype[site], alphabet)
  mutation = sample(possible_mutations)
  new_haplotype = string(haplotype[1:(site-1)],mutation,haplotype[(site+1):length(haplotype)])
  return new_haplotype
end

#mutations have fitness effects

function get_fitness(haplotype)
  old_fitness = fitness[haplotype]
  if rand() < fitness_chance
    return old_fitness * fitness_effect
  else
    return old_fitness
  end
end

#If a mutation event creates a new haplotype, assign it a random fitness

function mutation_event()
  haplotype = get_random_haplotype()
end




def mutation_event():
    haplotype = get_random_haplotype()
    if pop[haplotype] > 1:
        pop[haplotype] -= 1
        new_haplotype = get_mutant(haplotype)
        if new_haplotype in pop:
            pop[new_haplotype] += 1
        else:
            pop[new_haplotype] = 1
        if new_haplotype not in fitness:
            fitness[new_haplotype] = get_fitness(haplotype)

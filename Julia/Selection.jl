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
allele_frequencies = [0.5,0.5]

#Population of haplotypes maps to counts and fitnesses:
#Store this as a lightweight dictionary that maps a string to a count

pop = Dict("AAAAAAAAAA" => 40)
pop["AAATAAAAAA"] = 30
pop["AATTTAAAAA"] = 30

#Map haplotype string to fitness float
fitness = Dict("AAAAAAAAAA" => 1.0)
fitness["AAATAAAAAA"] = 1.05
fitness["AATTTAAAAA"] = 1.10

#Generate random haplotypes
function generate_haplotype()
  hap = wsample(alphabet,allele_frequencies)
  [hap = hap * wsample(alphabet,allele_frequencies) for i in 2:seq_length]
  return hap
end

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
  if pop[haplotype] > 1
    pop[haplotype] -= 1
    new_haplotype = get_mutant(haplotype)
      if new_haplotype in collect(keys(pop))
        pop[new_haplotype] += 1
      else
        pop[new_haplotype] = 1
      end
      if !(new_haplotype in collect(keys(fitness)))
        fitness[new_haplotype] = get_fitness(haplotype)
      end
    end
end

function mutation_step()
  mutation_count = get_mutation_count()
  for i in 1:mutation_count
    mutation_event()
  end
end

#Genetic drift and fitness affect which haplotypes make it to the next generation
#fitness weights the multinomial draw

function get_offspring_counts()
  haplotypes = collect(keys(pop))
  frequencies = [pop[haplotype] / pop_size for haplotype in haplotypes]
  fitnesses = [fitness[haplotype] for haplotype in haplotypes]
  weights = [x * y for (x,y) in zip(frequencies,fitnesses)]
  total = sum(weights)
  weights = [x / total for x in weights]
  return rand(Multinomial(pop_size, weights))
end

function offspring_step()
  counts = get_offspring_counts()
  for (haplotype,count) in zip(collect(keys(pop)),counts)
    if count > 0
      pop[haplotype] = count
    else
      delete!(pop,haplotype)
    end
  end
end

#Combine and iterate

function time_step()
  mutation_step()
  offspring_step()
end

generations = 5

function simulate()
  for i in 1:generations
    time_step()
  end
end

#Record
#We want to keep a record of past population frequencies to understand dynamics through time. At each step in the simulation, we append to a history object.

history = [] #init history

function simulate()
  clone_pop = pop
  append!(history,clone_pop)
  for i in 1:generations
    time_step()
    clone_pop = pop
    append!(history,clone_pop)
  end
end


#########################################################
#  Analyze Trajectories
##########################################################

#Calculate diversity
function get_distance(seq_a,seq_b)
  diffs = 0
  seqLen = length(seq_a)
  # assert length(seq_q) == length(seq_b) FIX LATER
  for (chr_a, chr_b) in zip(seq_a, seq_b)
    if chr_a != chr_b
      diffs += 1
    end
  end
  return diffs / seqLen
end

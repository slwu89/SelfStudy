#################################################################
#
#   Wright-Fisher model of mutation, selection and random genetic drift
#   Julia port of Python 2.5 code (SISMID 2016)
#   Sean Wu
#   January 11, 2017
#
#################################################################

using Distributions
using Gadfly
using DataFrames

#########################################################
#  Make population dynamic model
#########################################################

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

#Reset all global Parameters
function reset_pop()
  global pop, fitness, history
  pop = Dict("AAAAAAAAAA" => 40)
  pop["AAATAAAAAA"] = 30
  pop["AATTTAAAAA"] = 30
  fitness = Dict("AAAAAAAAAA" => 1.0)
  fitness["AAATAAAAAA"] = 1.05
  fitness["AATTTAAAAA"] = 1.10
  history = []
end

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
  local clone_pop
  clone_pop = deepcopy(pop)
  push!(history,clone_pop)
  for i in 1:generations
    time_step()
    clone_pop = deepcopy(pop)
    push!(history,clone_pop)
  end
end

simulate()

#########################################################
#  Analyze Trajectories
#########################################################

#Calculate hamming distance
function get_distance(seq_a,seq_b)
  diffs = 0
  seqLen = length(seq_a)
  assert(length(seq_a) == length(seq_b))
  for (chr_a, chr_b) in zip(seq_a, seq_b)
    if chr_a != chr_b
      diffs += 1
    end
  end
  return diffs / seqLen
end

#Calculate population level diversity
function get_diversity(population)
  haplotypes = collect(keys(population))
  haplotype_count = length(haplotypes)
  diversity = 0
  for i in 1:haplotype_count, j in 1:haplotype_count
    haplotype_a = haplotypes[i]
    haplotype_b = haplotypes[j]
    frequency_a = population[haplotype_a] / pop_size
    frequency_b = population[haplotype_b] / pop_size
    frequency_pair = frequency_a * frequency_b
    diversity += frequency_pair * get_distance(haplotype_a,haplotype_b)
  end
  return diversity
end

function get_diversity_trajectory()
  trajectory = [get_diversity(generation) for generation in history]
  return trajectory
end


#########################################################
#  Plot Diversity, Divergence, Haplotype & SNP Trajectories
#########################################################

#same plots

diversityDf = DataFrame(traj=get_diversity_trajectory())

plot(diversityDf, y=:traj, Geom.line, Geom.point,
  Theme(default_color=color("#447CCD")),Guide.XLabel("Generation"),Guide.YLabel("Diversity"))

###Analyze and plot Divergence###
function get_divergence(population)
  haplotypes = collect(keys(population))
  divergence = 0
  for haplotype in haplotypes
    frequency = population[haplotype] / pop_size
    divergence += frequency * get_distance(base_haplotype, haplotype)
  end
  return divergence
end

function get_divergence_trajectory()
  trajectory = [get_divergence(generation) for generation in history]
  return trajectory
end

#plot divergence
divergenceDf = DataFrame(traj=get_divergence_trajectory())

plot(divergenceDf, y=:traj, Geom.line, Geom.point,
  Theme(default_color=color("#447CCD")),Guide.XLabel("Generation"),Guide.YLabel("Divergence"))

  #########################################################
  #  Haplotype Trajectories
  #########################################################

function get_frequency(haplotype,generation)
  pop_at_generation = history[generation]
  if haplotype in collect(keys(pop_at_generation))
    return pop_at_generation[haplotype] / pop_size
  else
    return 0
  end
end

function get_trajectory(haplotype)
  trajectory = [get_frequency(haplotype,gen) for gen in 1:generations]
  return trajectory
end

function get_all_haplotypes()
  haplotypes = []
  for generation in history
    for haplotype in collect(keys(generation))
      push!(haplotypes,haplotype)
    end
  end
  return haplotypes
end

colors = ["#781C86", "#571EA2", "#462EB9", "#3F47C9", "#3F63CF", "#447CCD", "#4C90C0", "#56A0AE", "#63AC9A", "#72B485", "#83BA70", "#96BD60", "#AABD52", "#BDBB48", "#CEB541", "#DCAB3C", "#E49938", "#E68133", "#E4632E", "#DF4327", "#DB2122"]

colors_lighter = ["#A567AF", "#8F69C1", "#8474D1", "#7F85DB", "#7F97DF", "#82A8DD", "#88B5D5", "#8FC0C9", "#97C8BC", "#A1CDAD", "#ACD1A0", "#B9D395", "#C6D38C", "#D3D285", "#DECE81", "#E8C77D", "#EDBB7A", "#EEAB77", "#ED9773", "#EA816F", "#E76B6B"]

function stacked_trajectory_plot(xlabel="generation")
  haplotypes = get_all_haplotypes()
  trajectories = [get_trajectory(haplotype) for haplotype in haplotypes]
  #need to do the plotting here
end


# def stacked_trajectory_plot(xlabel="generation"):
#     mpl.rcParams['font.size']=18
#     haplotypes = get_all_haplotypes()
#     trajectories = [get_trajectory(haplotype) for haplotype in haplotypes]
#     plt.stackplot(range(generations), trajectories, colors=colors_lighter)
#     plt.ylim(0, 1)
#     plt.ylabel("frequency")
#     plt.xlabel(xlabel)

#########################################################
#  SNP Trajectories
#########################################################

function get_snp_frequency(site,generation)
  minor_allele_frequency = 0.0
  pop_at_generation = history[generation]
  for haplotype in collect(keys(pop_at_generation))
    allele = haplotype[site]
    frequency = pop_at_generation[haplotype] / pop_size
    if allele != 'A' #note that allele is char not string
      minor_allele_frequency += frequency
    end
  end
  return minor_allele_frequency
end

function get_snp_trajectory(site)
  trajectory = [get_snp_frequency(site,gen) for gen in 1:generations]
  return trajectory
end

function get_all_snps()
  snps = Set()
  for generation in history
    for haplotype in collect(keys(generation))
      for site in 1:seq_length
        if haplotype[site] != 'A'
          push!(snps,site)
        end
      end
    end
  end
  return snps
end

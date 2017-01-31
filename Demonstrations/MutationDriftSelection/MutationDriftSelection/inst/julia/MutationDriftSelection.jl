#################################################################
#
#   Wright-Fisher model of mutation, selection and random genetic drift
#   Julia port of Python 2.5 code (SISMID 2016)
#   Julia source for R Shiny app
#   Sean Wu
#   January 28, 2017
#
#################################################################


#########################################################
#  MutationDriftSelection module:
#
#  define module to be loaded into R session
#  through XR, XRJulia interfaces
#  note: pars is dictionary of parameters
#  passed from the R session.
#########################################################

module MutationDriftSelection

using Distributions

alphabet = ['A','T']

#Generate random haplotypes
function generate_haplotype(pars)
  hap = wsample(alphabet,pars[:alleleFreq])
  [hap = hap * wsample(alphabet,pars[:alleleFreq]) for i in 2:pars[:seqLen]]
  return hap
end

export generate_haplotype

end

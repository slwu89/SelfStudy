#################################################################
#
#   Pop Gen for Insect Vectors
#   Sean Wu
#   February 9, 2017
#
#################################################################

# #qPrime: recurrence for q (gametic frequency of the HEG allele)
# #e: homing rate (probability of successful gene conversion)
# #s: fitness cost of disrupting gene function for homozygote
# #sh: fitness cost of disrupting gene function for heterozygote
# #q: gametic frequency of HEG allele
# #p: gametic frequency of wild-type allele
# qPrime <- function(e,s,h,q,p){
#   (((1-s)*q^2) + ((1-(s*h))*p*q*(1+e))) / (1 - (s*q^2) - (2*s*h*p*q))
# }
# 
# qPrime(e = 0.99,s = 0.1,h = 1,q = 0.05,p = 0.95)


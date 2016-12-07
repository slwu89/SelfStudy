#################################################################
#
#   Power laws, Pareto distributions and Zipf's law
#   MEJ Newman, 2004
#   Sean Wu
#   December 7, 2016
#
#################################################################

library(poweRlaw)

##########################################
# Main Functions
##########################################

data = rplcon(n=1e4,xmin=1,alpha=2.5)
hist(log(data),freq=F)
plot(ecdf(data))


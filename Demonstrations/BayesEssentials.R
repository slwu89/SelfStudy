#################################################################
#
#   Bayesian Essentials with R
#   Chapter 2
#   Sean Wu
#   February 4, 2017
#
#################################################################

library(bayess)

qqnorm((shift-mean(shift))/sd(shift),pch=19,col="gold2")
abline(a=0,b=1,lty=2,col="indianred",lwd=2)

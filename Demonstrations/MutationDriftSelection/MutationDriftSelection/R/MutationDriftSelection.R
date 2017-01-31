####################################################
#
# Mutation, Drift, Selection
# Sean Wu
#
#
#
####################################################

test <- function(){
  x = rnorm(50,0,1)
  print(summary(x))
}

testJL <- XRJulia::juliaEval("reverse(%s)", 1:5)

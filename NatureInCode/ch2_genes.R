############################
######Chapter 2. Genes######
######Sean Wu 5/27/2016#####
############################


####################
###Hardy-Weinberg###
####################
a1a1 <- 0.15
a2a2 <- 0.35
a1a2 <- 1 - (a1a1 + a2a2)
p <- a1a1 + (a1a2/2)
q <- 1-p

calculate_allele_frequencies <- function(a1a1,a2a2,a1a2){
  p <- a1a1 + (a1a2/2)
  q <- 1-p
  return(c(p=p,q=q))
}

create_next_generation <- function(a1a1,a2a2,a1a2){
  allele_freq <- calculate_allele_frequencies(a1a1=a1a1,a2a2=a2a2,a1a2=a1a2)
  p <- allele_freq[["p"]]
  q <- allele_freq[["q"]]
  a1a1 = p^2
  a1a2 = 2 * p * q
  a2a2 = q^2
  return(c(a1a1=a1a1,a1a2=a1a2,a2a2=a2a2))
}

for(i in 1:10){
  if(i==1){
    print(paste("Generation i: 0","A1A1:",.GlobalEnv$a1a1,"A2A2:",.GlobalEnv$a2a2,"A1A2",.GlobalEnv$a1a2))
    gen_i <- create_next_generation(.GlobalEnv$a1a1,.GlobalEnv$a2a2,.GlobalEnv$a1a2)
  } else {
    gen_i <- create_next_generation(gen_i[["a1a1"]],gen_i[["a2a2"]],gen_i[["a1a2"]])
  }
  print(paste("Generation i:",i,"A1A1:",gen_i[["a1a1"]],"A2A2:",gen_i[["a2a2"]],"A1A2",gen_i[["a1a2"]]))
}
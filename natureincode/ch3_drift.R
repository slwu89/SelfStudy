####################################
######Chapter 3. Genetic Drift######
######Sean Wu 6/15/2016#############
####################################

#The Randomness of Finite Populations

p <- 0.5
N <- 1e3
generations <- 1e3

next_generation <- function(){
  draws <- 2*N
  a1 <- 0
  a2 <- 0
  for(i in 1:draws){
    if(runif(1) <= .GlobalEnv$p){
      a1 <- a1 + 1
    } else {
      a2 <- a2 + 1
    }
  }
  .GlobalEnv$p <- a1 / draws
}

for(i in 1:generations){
  next_generation()
  print(paste("generation",i,"p:",round(p,3),"q:",round(1-p,3)))
}

#Visualizing Drift

library(ggplot2)
library(reshape2)

p <- 0.5
N <- 1e3
generations <- 1e3
data <- NULL

next_generation <- function(){
  draws <- 2*N
  a1 <- 0
  a2 <- 0
  for(i in 1:draws){
    if(runif(1) <= .GlobalEnv$p){
      a1 <- a1 + 1
    } else {
      a2 <- a2 + 1
    }
  }
  .GlobalEnv$p <- a1 / draws
  .GlobalEnv$data <- c(.GlobalEnv$data,.GlobalEnv$p)
}

for(i in 1:generations){
  next_generation()
}

data <- data.frame(gen=1:1e3,p=data)

ggplot() + 
  geom_line(data=data,aes(x=gen,y=p),colour="steelblue") +
  labs(x="Generation",y="p") +
  theme_bw()

#Multiple Simulation Runs

one_run <- function(p,N,generations){
  out <- NULL
  for(i in 1:generations){
    draws <- 2*N
    a1 <- 0
    a2 <- 0
    for(j in 1:draws){
      if(runif(1) <= p){
        a1 <- a1 + 1
      } else {
        a2 <- a2 + 1
      }
    }
    p <- a1 / draws
    out[i] <- p
  }
  return(out)
}

library(parallel)
library(doSNOW)
library(foreach)

cl <- makeCluster(spec=detectCores())
registerDoSNOW(cl)

all_runs <- function(n_run,p,N,generations){
  
  out <- foreach(i=1:n_run,.verbose=TRUE,.export=c("one_run")) %dopar% {
    one_run(p=p,N=N,generations=generations)
  }
  return(out)
}

all_data <- all_runs(n_run=10,p=0.5,N=2e3,generations=1e3)
all_data <- melt(all_data)
all_data$gen <- rep(1:1e3,10)

ggplot() +
  geom_line(data=all_data,aes(x=gen,y=value,group=L1,colour=as.factor(L1))) +
  guides(colour=FALSE) +
  labs(x="Generation",y="p") +
  theme_bw()

#Effective Population Size

one_run_Ne <- function(p,N,generations){
  out <- list(p=NULL,N=NULL)
  for(i in 1:generations){
    if(i %% 10 == 0){
      N_i <- 10
    } else {
      N_i <- N
    }
    draws <- 2*N_i
    a1 <- 0
    a2 <- 0
    for(j in 1:draws){
      if(runif(1) <= p){
        a1 <- a1 + 1
      } else {
        a2 <- a2 + 1
      }
    }
    p <- a1 / draws
    out$p[i] <- p
    out$N[i] <- N_i
  }
  return(out)
}

all_runs_Ne <- function(n_run,p,N,generations){
  
  out <- foreach(i=1:n_run,.verbose=TRUE,.export=c("one_run_Ne")) %dopar% {
    one_run_Ne(p=p,N=N,generations=generations)
  }
  return(out)
  
}

all_data_Ne <- all_runs_Ne(n_run=10,p=0.5,N=2e3,generations=1e3)
all_data_Ne <- melt(all_data_Ne)
all_data_Ne$gen <- rep(1:1e3,10)

ggplot() +
  geom_line(data=all_data_Ne[all_data_Ne$L2=="p",],aes(x=gen,y=value,group=L1,colour=as.factor(L1))) +
  guides(colour=FALSE) +
  labs(x="Generation",y="p") +
  theme_bw()

#Unequal Sex Ratio

Ne_N_ratio <- function(Nm_Nf){
  ans <- 4 * (Nm_Nf) * (1-Nm_Nf)
  return(ans)
}

Ne_N_data <- sapply(seq(from=0,to=1,by=1e-3),Ne_N_ratio)
Ne_N_data <- data.frame(Ne_N=Ne_N_data,PropotionM=seq(from=0,to=1,by=1e-3))

ggplot() +
  geom_line(data=Ne_N_data,aes(y=Ne_N,x=PropotionM),colour="steelblue") +
  theme_bw()

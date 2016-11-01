################################################
######Chapter 5. Migration: Spatial Models######
######Sean Wu 5/27/2016#########################
################################################

library(Rcpp)

####################
###Create lattice###
####################

#A1A1=1, A2A2=2, A1A2=3
cell_init <- function(i,j,p){
  rand <- runif(1)
  if(rand < p^2){
    return(1)
  } else if(rand > 1-((1-p)^2)) {
    return(2)
  } else {
    return(3)
  }
}

cell_initV <- Vectorize(cell_init,c("i","j"))

lattice_init <- function(r,c,p){
  outer(1:r,1:c,cell_initV,p=p)
}


#########################################
###Helper functions to run simulations###
#########################################

#functions to return an index in the viable mating range, with bounding at borders to account for toroidal lattice structure
get_index <- function(i,d){
  min <- i-d
  max <- i+d
  return(floor(runif(1) * (max - min + 1)) + (min))
}

#rewrite in C++
cppFunction("int get_indexC(int i, int d){
  int ans;
  double random;
  random = ((double)rand()/(double)RAND_MAX);
  int min;
  int max;
  min = i - d;
  max = i + d;
  ans = floor(random * (max-min+1)) + min;
  return(ans);
}")

bounded_index <- function(i,d,dim){
  index <- get_index(i=i,d=d)
  if(index <= 0){
    index <- index + dim
  }
  if(index > dim){
    index <- index - dim
  }
  return(index)
}

#rewrite in C++
cppFunction("int bounded_indexC(int i, int d, int dim){
  // create index
  int index;
  double random;
  random = ((double)rand()/(double)RAND_MAX);
  int min;
  int max;
  min = i - d;
  max = i + d;
  index = floor(random * (max-min+1)) + min;
  // bounded index
  int ans;
  if(index <= 0){
    ans = index + dim;
  } else if(index > dim){
      ans = index - dim;
  } else {
    ans = index;
  }
  return(ans);
}")

#function to return a mate for the cell at position i,j
get_mate <- function(i,j,d,grid){
  i_prime <- bounded_index(i=i,d=d,dim=nrow(grid))
  j_prime <- bounded_index(i=j,d=d,dim=ncol(grid))
  return(grid[i_prime,j_prime])
}

#use C++ functions
get_mateC <- function(i,j,d,grid){
  i_prime <- bounded_indexC(i=i,d=d,dim=nrow(grid))
  j_prime <- bounded_indexC(i=j,d=d,dim=ncol(grid))
  return(grid[i_prime,j_prime])
}

#function to return offspring for two parents (A1A1=1, A2A2=2, A1A2=3)
get_offspring <- function(par1,par2){
  if(par1==1 & par2==1){ #both A1A1
    return(1)
  } else if((par1==1 & par2==3) | (par1==3 & par2==1)){ #one A1A1 one A1A2
    if(runif(1) < 0.5){
      return(1)
    } else {
      return(3)
    }
  } else if((par1==1 & par2==2) | (par1==2 & par2==1)){ #one A1A1 one A2A2
    return(3)
  } else if(par1==3 & par2==3){ #both A1A2
    rand <- runif(1)
    if(rand < 0.25){
      return(1)
    } else if(rand > 0.75){
      return(2)
    } else {
      return(3)
    }
  } else if((par1==2 & par2==3) | (par1==3 & par2==2)){ #one A1A2 one A2A2
    if(runif(1) < 0.5){
      return(2)
    } else {
      return(3)
    }
  } else { #both A2A2
    return(2)
  }
}

#rewrite in C++
cppFunction("int get_offspringC(int par1, int par2){
  int offspring;
  double random;
  random = ((double)rand()/(double)RAND_MAX);
  if(par1 == 1 && par2 == 1) { //both A1A1
    offspring = 1;
  } else if((par1 == 1 && par2 == 3) || (par1 == 3 && par2 == 1)) { //one A1A1 one A1A2
    if(random < 0.5) {
      offspring = 1;
    } else {
      offspring = 3;
    }
  } else if((par1 == 1 && par2 == 2) || (par1 == 2 & par2 == 1)) { //one A1A1 one A2A2
    offspring = 3;
  } else if(par1 == 3 && par2 == 3) { //both A1A2
    if(random < 0.25){
      offspring = 1;
    } else if(random > 0.75){
      offspring = 2;
    } else {
      offspring = 3;
    }
  } else if((par1 == 2 && par2 == 3) || (par1 == 3 & par2 == 2)) { //one A1A2 one A2A2
    if(random < 0.5){
      offspring = 2;
    } else {
      offspring = 3;
    }
  } else { //both A2A2
    offspring = 2;
  }
  return(offspring);
}")

#function that selects a mate and produces an offspring for a single cell
cell_generation <- function(i,j,grid,d){
  mate <- get_mate(i=i,j=j,d=d,grid=grid)
  offspring <- get_offspring(par1=grid[i,j],par2=mate)
  return(offspring)
}

#use C++ functions
cell_generationC <- function(i,j,grid,d){
  mate <- get_mateC(i=i,j=j,d=d,grid=grid)
  offspring <- get_offspringC(par1=grid[i,j],par2=mate)
  return(offspring)
}

#function to run a single generation
run_generation <- function(grid,d=1){
  dim <- nrow(grid)
  iterator <- expand.grid(1:dim,1:dim) #set up iterator for mapply 
  new_grid <- mapply(cell_generation,i=iterator$Var1,j=iterator$Var2,MoreArgs=list(grid=grid,d=d),SIMPLIFY=TRUE) #use mapply to create offspring grid based on original grid
  new_grid <- matrix(data=new_grid,nrow=dim,ncol=dim) #need to turn vector into matrix; R lists matricies in column-major order
  return(new_grid)
}

#use C++ functions
run_generationC <- function(grid,d=1){
  dim <- nrow(grid)
  iterator <- expand.grid(1:dim,1:dim) #set up iterator for mapply 
  new_grid <- mapply(cell_generationC,i=iterator$Var1,j=iterator$Var2,MoreArgs=list(grid=grid,d=d),SIMPLIFY=TRUE) #use mapply to create offspring grid based on original grid
  new_grid <- matrix(data=new_grid,nrow=dim,ncol=dim) #need to turn vector into matrix; R lists matricies in column-major order
  return(new_grid)
}

# run_generation <- function(grid){
#   new_grid <- matrix(NA,nrow=nrow(grid),ncol=ncol(grid))
#   for(i in 1:nrow(grid)){
#     for(j in 1:ncol(grid)){
#       new_grid[i,j] <- cell_generation(i=i,j=j,grid=grid)
#     }
#   }
#   return(new_grid)
# }


#####################
###Run Simulations###
#####################


init_grid <- lattice_init(r=100,c=100,p=0.5)
table(init_grid) #check to make sure allele frequencies make sense
sim_out <- list()

for(i in 1:200){
  if(i==1){
    grid_i <- run_generationC(init_grid)
    sim_out[[i]] <- grid_i
  }
  grid_i <- run_generationC(grid=grid_i)
  sim_out[[i]] <- grid_i
  print(paste("Generation i:",i))
}

image(init_grid,col=c("darkblue","steelblue","white"),axes=FALSE)
image(sim_out[[100]],col=c("darkblue","steelblue","white"),axes=FALSE)

IMAGE.dir <- "C:/Users/WuS/Desktop/migration_out/"
for(i in 1:length(sim_out)){
  write.table(x=sim_out[[i]],file=paste0(IMAGE.dir,"time",i,".csv"),col.names=FALSE,sep=",")
}

saveGIF(expr={
  for(i in 1:length(sim_out)){
    print(paste0("plotting iteration: ",i))
    image(sim_out[[i]],col=c("darkblue","steelblue","white"),axes=FALSE)
  }
},
movie.name = "/Users/slwu89/Desktop/migration.gif",
ani.height=640,ani.width=640,interval=0.1)

#plot model output with ggplot2
# library(ggplot2)
# library(reshape2)
#   
# gg_image <- function(grid){
#   melt_dat <- melt(grid)
#   plot <- ggplot() +
#     geom_raster(data=melt_dat,aes(x=Var1,y=Var2,fill=as.factor(value))) +
#     scale_fill_manual(values=c("1"="darkblue","2"="steelblue","3"="white")) +
#     guides(fill=FALSE) +
#     theme(axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.title=element_blank(),title=element_blank(),panel.border=element_rect(fill=NA,colour="black",size=2))
#   print(plot)
# }
#   
# gg_image(init_grid)
# gg_image(sim_out[[100]])


################
###Inbreeding###
################


get_inbreed <- function(grid){
  A1A1 <- sum(grid==1)
  A2A2 <- sum(grid==2)
  A1A2 <- sum(grid==3)
  N <- sum(A1A1,A2A2,A1A2)
  p <- ((2*A1A1)+A1A2) / (2*N)
  h_0 <- A1A2/N #observed heterozygosity
  h_e <- 2 * p * (1-p) #expected heterozygosity
  inbreed_coef <- (h_e - h_0) / h_e
  return(inbreed_coef)
}

sapply(sim_out,get_inbreed)


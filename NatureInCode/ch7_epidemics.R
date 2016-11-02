################################
######Chapter 7. Epidemics######
######Sean Wu 5/27/2016#########
################################

library(Rcpp)
library(animation)

####################
###Create lattice###
####################

#S=1, I=2, R=3
init_grid <- matrix(data=1,nrow=100,ncol=100)
init_grid[runif(1,1,100),runif(1,1,100)] <- 2

#bounded index
cppFunction("int bound_index(int i, int dim){
  int index;
  if(i <= 0){
    index = i + dim;
  } else if(i > dim){
    index = i - dim;
  } else {
    index = i;
  }
  return(index);
}")

#new_grid will be updated at each i,j ; grid does not change for each iteration. one iteration goes through all i,j (entire lattice)
#note that this function only needs to be called if original grid[i,j] == 2
expose_neighbors <- function(i,j,beta,grid,new_grid){
  
  #set up iterators for neighboring cells
  i_range <- (i-1):(i+1)
  j_range <- (j-1):(j+1)
  i_range <- sapply(i_range,bound_index,dim=nrow(grid))
  j_range <- sapply(j_range,bound_index,dim=ncol(grid))
  
  #iterate through neighboring cells
  for(ii in i_range){
    for(jj in j_range){
      if(ii==i & jj==j){ #do not need to iterate over focal cell
        next()
      }
      if(grid[ii,jj]==1){ #only infect current susceptibles
        if(runif(1) < beta){
          new_grid[ii,jj] <- 2
        }
      }
    }
  }
  
  return(new_grid)
}

expose_neighbors_smallWorld <- function(i,j,beta,grid,new_grid){
  
  #set up iterators for neighboring cells
  i_range <- (i-1):(i+1)
  j_range <- (j-1):(j+1)
  i_range <- sapply(i_range,bound_index,dim=nrow(grid))
  j_range <- sapply(j_range,bound_index,dim=ncol(grid))
  
  #iterate through neighboring cells
  for(ii in i_range){
    for(jj in j_range){
      if(ii==i & jj==j){ #do not need to iterate over focal cell
        next()
      }
      if(runif(1) < 0.01){ #small world contact
        rand_i <- bound_index(i=runif(1,min=1,max=nrow(grid)),dim=nrow(grid))
        rand_j <- bound_index(i=runif(1,min=1,max=ncol(grid)),dim=ncol(grid))
        if(grid[rand_i,rand_j]==1){ #only infect current susceptibles
          if(runif(1) < beta){
            new_grid[rand_i,rand_j] <- 2
          }
        }
      } else {
        if(grid[ii,jj]==1){ #only infect current susceptibles
          if(runif(1) < beta){
            new_grid[ii,jj] <- 2
          }
        }
      }
    }
  }
  
  return(new_grid)
}

#run a single iteration
run_time_step <- function(beta,gamma,grid){
  new_grid <- grid
  for(i in 1:nrow(grid)){
    for(j in 1:ncol(grid)){
      if(grid[i,j]==2){ #only run neighbor infection & recovery routine if focal cell is currently infected
        new_grid <- expose_neighbors_smallWorld(i,j,beta,grid,new_grid) #infection
        if(runif(1) < gamma){
          new_grid[i,j] <- 3
        }
      } 
    }
  }
  return(new_grid)
}


#run simulation
sim_out <- list()

for(i in 1:200){
  if(i==1){
    grid_i <- run_time_step(beta=0.05,gamma=0.1,init_grid)
    sim_out[[i]] <- grid_i
    print(paste("Day i:",i))
  }
  grid_i <- run_time_step(beta=0.05,gamma=0.1,grid_i)
  sim_out[[i]] <- grid_i
  print(paste("Day i:",i))
}

saveGIF(expr={
    for(i in 1:length(sim_out)){
      image(sim_out[[i]],col=c("#dcdcdc","#c82605","#6fc041"),axes=FALSE)
    }
  },
  movie.name = "/Users/slwu89/Desktop/epidemic.gif",
  ani.height=640,ani.width=640)

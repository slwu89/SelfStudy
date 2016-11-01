#include <Rcpp.h>
using namespace Rcpp;

/*
 * Chapter 7. Epidemics
 * Sean Wu 10/7/2016
 */

//container for output
typedef std::vector<Rcpp::CharacterMatrix> output_obj;

// [[Rcpp::export]]
int get_bounded_index(int index, int dim){
  int bounded_index = index;
  if(index < 0){
    bounded_index = index + dim;
  }
  if(index >= dim){
    bounded_index = index - dim;
  }
  return(bounded_index);
}

// [[Rcpp::export]]
CharacterMatrix run_time_step(double beta, double gamma, CharacterMatrix grid){
  
  //copy initial grid; use clone to make true copy rather than pass pointer
  CharacterMatrix temp_grid = Rcpp::clone(grid); 
  
  for(int i=0; i<grid.nrow(); i++){
    for(int ii=0; ii<grid.ncol(); ii++){
      // Rcout << grid(i,ii) << std::endl;
      if(grid(i,ii) == "I"){
        // Rcout << "i: " << i << ", ii: " << ii << ", is infected!!!" << std::endl;
        //expose neighbors
        for(int n_i=i-1; n_i<=i+1; n_i++){
          for(int n_ii=ii-1; n_ii<=ii+1; n_ii++){
            if(n_i == i && n_ii == ii){
              continue;
            } else {
              //try infection
              int bound_n_i = get_bounded_index(n_i,grid.nrow());
              int bound_n_ii = get_bounded_index(n_ii,grid.ncol());
              if(grid(bound_n_i,bound_n_ii) == "S"){
                if(R::runif(0.0,1.0) < beta){
                  // Rcout << "infection!!" << std::endl;
                  temp_grid(bound_n_i,bound_n_ii) = "I";
                }
              }
            }
          }
        }
        //try recovery
        if(R::runif(0.0,1.0) < gamma){
          temp_grid(i,ii) = "R";
        }
      } else {
        continue;
      }
    }
  }
  
  return(temp_grid);
}

// [[Rcpp::export]]
CharacterMatrix run_time_step_small(double beta, double gamma, CharacterMatrix grid){
  
  //copy initial grid; use clone to make true copy rather than pass pointer
  CharacterMatrix temp_grid = Rcpp::clone(grid); 
  
  for(int i=0; i<grid.nrow(); i++){
    for(int ii=0; ii<grid.ncol(); ii++){
      if(grid(i,ii) == "I"){
        //expose neighbors
        for(int n_i=i-1; n_i<=i+1; n_i++){
          for(int n_ii=ii-1; n_ii<=ii+1; n_ii++){
            if(n_i == i && n_ii == ii){
              continue;
            } else {
              //random contact
              if(R::runif(0.0,1.0) < 0.01){
                int random_i = get_bounded_index(floor(R::runif(0.0,grid.nrow())),grid.nrow());
                int random_j = get_bounded_index(floor(R::runif(0.0,grid.ncol())),grid.ncol());
                if(grid(random_i,random_j) == "S"){
                  if(R::runif(0.0,1.0) < beta){
                    temp_grid(random_i,random_j) = "I";
                  }
                }
              } else {
                //neighbor contact
                int bound_n_i = get_bounded_index(n_i,grid.nrow());
                int bound_n_ii = get_bounded_index(n_ii,grid.ncol());
                if(grid(bound_n_i,bound_n_ii) == "S"){
                  if(R::runif(0.0,1.0) < beta){
                    temp_grid(bound_n_i,bound_n_ii) = "I";
                  }
                }
              }
            }
          }
        }
        //try recovery
        if(R::runif(0.0,1.0) < gamma){
          temp_grid(i,ii) = "R";
        }
      } else {
        continue;
      }
    }
  }
  
  return(temp_grid);
}

// [[Rcpp::export]]
List run_simulation(int t_end,double beta, double gamma, CharacterMatrix init_grid, bool small=true){
  
  //container for output
  output_obj out;
  // out.reserve(t_end);
  out.push_back(init_grid);
  
  CharacterMatrix grid_i;
  
  for(int i=0; i<t_end; i++){
    
    Rcout << "on iteration: " << (i+1) << ", of: " << t_end << std::endl;
    
    if(i==0){
      if(small){
        grid_i = run_time_step_small(beta,gamma,init_grid);
      } else {
        grid_i = run_time_step(beta,gamma,init_grid);
      }
      out.push_back(grid_i);
    } else {
      if(small){
        grid_i = run_time_step_small(beta,gamma,grid_i);
      } else {
        grid_i = run_time_step(beta,gamma,grid_i);
      }
      out.push_back(grid_i);
    }

  }
  
  return(List::create(Named("output")=out));
}



/***R
lattice <- matrix("S",nrow=300,ncol=300)
lattice[ceiling(runif(1,0,10)),ceiling(runif(1,0,10))] <- "I"

beta <- 0.055
gamma <- 0.125

sim_outC <- run_simulation(t_end = 300,beta = 0.05,gamma = 0.1,init_grid = lattice)


# sim_out <- list()
# sim_out[[1]] <- lattice
# for(i in 1:10){
#   print(paste0("i: ",i))
#   if(i==1){
#     grid_i <- run_time_step(beta,gamma,lattice)
#     plot_lattice(grid_i)
#   } else {
#     grid_i <- run_time_step(beta,gamma,grid_i)
#     plot_lattice(grid_i)
#   }
#   sim_out[[(i+1)]] <- grid_i
# 
# }


#save default plotting parameters
default <- par()

#function to plot character matrix lattice
plot_lattice <- function(lattice){
  gg_col <- colorRampPalette(colors=c("#132B43","#56B1F7"))(3)
  lattice_plot <- apply(lattice,c(1,2),function(x){
    if(x=="S"){
      return(1)
    } else if(x=="I"){
      return(2)
    } else {
      return(3)
    }
  })
  par(mar=c(0,0,0,0),mgp=c(0,0,0))
  image(lattice_plot,col=gg_col,axes=FALSE)
}

#animate the output
library(animation)

saveGIF(expr={
  for(i in 1:length(sim_outC$output)){
    print(paste0("plotting i: ",i))
    plot_lattice(sim_outC$output[[i]])
    text(x=.9,y=.95,labels=paste0("time: ",i),col="white")
  }
},movie.name="/Users/slwu89/Desktop/epidemic.gif",ani.height=768,ani.width=768,
interval=0.1)

#restore default plotting parameters
par(default)
*/
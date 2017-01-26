#################################################################
#
#   Mushroom Hunt
#   Soc 273A
#   Sean Wu
#   January 26, 2017
#
#################################################################

using Distributions

function poisson_point_process(rectangular_bounds, cnt)
    left, right, down, up=rectangular_bounds
    width_scale=right-left
    length_scale=up-down
    x=zeros(Float64, 2, cnt)
    for draw_idx in 1:cnt
    	x[:,draw_idx]=[left+rand()*width_scale, down+rand()*length_scale]
    end
    x
end

##########################################
# Distribute Points
##########################################

# function pointsPoisson(n,)

abstract agent

type turtle

end

module DynamicalBilliardsVisualizations

using DynamicalBilliards, Makie

import DynamicalBilliards: bdplot, bdplot!, bdplot_animstep!, bdplot_interactive, bdplot_video, bdplot_boundarymap

include("defs_plotting.jl")
include("defs_animating.jl")
include("premade_anim_functions.jl")

end
"""
A Julia package for dynamical billiard systems in two dimensions.

The goals of the package is to provide a flexible, easy-to-use and intuitive framework for
fast implementation of billiard systems of arbitrary construction.
"""
module DynamicalBilliards

include("ParticlesObstacles.jl")
include("Propagation.jl")
include("StandardBilliards.jl")
include("PlotBilliards.jl")
include("RaySplitting.jl")

end

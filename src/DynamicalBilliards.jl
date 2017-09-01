"""
A Julia package for dynamical billiard systems in two dimensions.

The goals of the package is to provide a flexible, easy-to-use and intuitive framework for
fast implementation of billiard systems of arbitrary construction.
"""
module DynamicalBilliards

using StaticArrays
import Base.show

##########################################
# Core                                   #
##########################################
include("ParticlesObstacles.jl")
include("Propagation.jl")
include("StandardBilliards.jl")
include("RaySplitting.jl")
include("LyapunovSpectrum.jl")

##########################################
# Test Options                           #
##########################################
"""
    test_options(;print_info = false, long_tests = false)
Set if you want the long version of the tests and if you want information to be
printed during testing.
"""
function test_options(;print_info::Bool = false, long_tests::Bool = false)
  ENV["DYNAMICALBILLIARDS_PRINTTEST"] = print_info
  ENV["DYNAMICALBILLIARDS_LONGTEST"] = long_tests
end
##########################################
# Plotting Routines (loaded on demand)   #
##########################################
"""
    (DynamicalBilliards.) enableplotting()
Enable plotting for the package DynamicalBilliards.jl
"""
function enableplotting()
  dir = joinpath(dirname(@__FILE__), "plotting")
  for f in readdir(dir)
    include(joinpath(dir, f))
  end
end

end#module

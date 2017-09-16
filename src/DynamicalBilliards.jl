__precompile__()

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
include("RaySplitting.jl")
include("StandardBilliards.jl")
include("LyapunovSpectrum.jl")

##########################################
# Test Options                           #
##########################################
"""
    test_options(;print_info = false, long_tests = false)
Set if you want the long version of the tests and if you want information to be
printed during testing.
"""
function test_options(;print_info::Bool = true, long_tests::Bool = true)
  ENV["DYNAMICALBILLIARDS_PRINTTEST"] = print_info
  ENV["DYNAMICALBILLIARDS_LONGTEST"] = long_tests
end
##########################################
# Plotting Routines (loaded on demand)   #
##########################################
"""
    DynamicalBilliards.enableplotting()
Enable plotting for the package DynamicalBilliards.jl. Requires
`using PyPlot` to work properly.
"""
function enableplotting()
  dir = joinpath(dirname(@__FILE__), "plotting")
  for f in readdir(dir)
    include(joinpath(dir, f))
  end
end

#
# bt, ray = billiard_raysplitting_showcase(3, 1, 0.3, 0.2)
# p = randominside(bt)
# p.pos = SVector{2}(0.2115211414442486, 0.3892293033199159)
# p.vel = SVector{2}(0.907742279866133, 0.4195282509479374)
# evolve!(p, bt, 1000.0, ray)
# println("Done")

# enableplotting()
# plot_billiard(bt)
# animate_evolution(p, bt, 200, ray)

end#module

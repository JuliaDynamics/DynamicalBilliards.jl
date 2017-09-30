__precompile__()

"""
A Julia package for dynamical billiard systems in two dimensions.

The goals of the package is to provide a flexible, easy-to-use
and intuitive framework for
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

# dir = joinpath(dirname((dirname(@__FILE__))), "test")
# for f in readdir(dir)
#     f == "runtests.jl" && continue
#     include(joinpath(dir, f))
# end

# ω = big(0.2)
# (r, x, y) = big.([0.4, 1.0, 1.0])
# bt = billiard_sinai(r, x, y; setting="periodic")
# xmin, ymin, xmax, ymax = cellsize(bt)
# d = bt[5]
# c = d.c
# tt=1000.0
# invalid = 0
# minddist = min(x, y)
# p = randominside(ω, bt)

#tmin, i = next_collision(p, bt)

# lyapunov_spectrum(1)

end#module

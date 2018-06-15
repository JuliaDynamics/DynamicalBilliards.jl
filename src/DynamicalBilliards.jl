__precompile__()

"""
A Julia package for dynamical billiard systems in two dimensions.

The goals of the package is to provide a flexible, easy-to-use
and intuitive framework for
fast implementation of billiard systems of arbitrary construction.
"""
module DynamicalBilliards

using LinearAlgebra
using StaticArrays
import Base: show, eltype, getindex

const SV = SVector{2}
export SVector

cossin(a) = ((x, y) = sincos(a); (y, x))

##########################################
# Core                                   #
##########################################
include("billiards/particles.jl")
include("billiards/obstacles.jl")
include("billiards/billiardtable.jl")
include("billiards/standard_billiards.jl")

include("timeevolution/collisiontimes.jl")
include("timeevolution/propagation.jl")

include("highleveltimes.jl")
include("poincaresos.jl")
include("lyapunov_spectrum.jl")

include("mushroomtools.jl")
export MushroomTools

include("raysplitting.jl")

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

####################################################
# Plotting Routines (loaded when `Using PyPlot`)   #
####################################################
using Requires
@require PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee" begin
    dir = joinpath(@__DIR__, "plotting")
    for f in readdir(dir)
        include(joinpath(dir, f))
    end
end










end#module

__precompile__()

"""
A Julia package for dynamical billiard systems in two dimensions.

It provides a flexible, easy-to-use, and intuitive framework for
fast implementation of billiards of arbitrary construction.
"""
module DynamicalBilliards

using LinearAlgebra
using StaticArrays
using Elliptic
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

include("timeevolution/collisions.jl")
include("timeevolution/propagation.jl")
include("timeevolution/timeseries.jl")
include("timeevolution/highleveltimes.jl")

include("plotting_api.jl")

##########################################
# Advanced                               #
##########################################
include("boundary/boundarymap.jl")
include("boundary/phasespacetools.jl")

include("poincare.jl")
include("lyapunov_spectrum.jl")

include("mushroomtools.jl")
export MushroomTools

include("raysplitting.jl")
include("parallel.jl")
include("testing.jl")

###################
# Update messages #
###################
using Scratch
display_update = true
version_number = "4.1"
update_name = "update_v$(version_number)"

function __init__()
if display_update
    # Get scratch space for this package
    versions_dir = @get_scratch!("versions")
    if !isfile(joinpath(versions_dir, update_name))
        printstyled(
            stdout,
            """
            \nUpdate message: DynamicalBilliards v$(version_number)
            Plotting & animating is now inside DynamicalBilliards.jl again,
            utilizing Julia Package Extension systems. This means that Julia versions 1.9+
            are supported only.
            """;
            color = :light_magenta,
        )
        touch(joinpath(versions_dir, update_name))
    end
end
end


end#module

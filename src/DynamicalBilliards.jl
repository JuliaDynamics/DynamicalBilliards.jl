"""
A Julia package for dynamical billiard systems in two dimensions.

The goals of the package is to provide a flexible and intuitive framework for fast
implementation of billiard systems of arbitrary construction.
"""
module DynamicalBilliards

include("ParticlesObstacles.jl")
include("Propagation.jl")
include("StandardBilliards.jl")
include("PlotBilliards.jl")

export Particle, MagneticParticle, cyclotron, Disk, Circle, Antidot, FiniteWall, Obstacle,
       NullWall, PeriodicWall, resolvecollision!, randominside,
       billiard_rectangle, billiard_sinai, billiard_sinai_periodic,
       collisiontime, propagate!, evolve!, construct, plot_particle, plot_billiard,
       plot_evolution, plot_obstacle, plot_cyclotron, specular!, normalvec
end

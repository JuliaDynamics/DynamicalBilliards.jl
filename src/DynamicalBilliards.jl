module DynamicalBilliards

include("Collisions.jl")
include("Propagation.jl")
#include("StandardBilliards.jl")
include("PlotBilliards.jl")

export Particle, MagneticParticle, cyclotron, Disk, Circle, Antidot, FiniteWall, Obstacle,
       NullWall, PeriodicWall, normalvec, resolvecollision!, distance, randominside,
       billiard_rectangle, billiard_sinai, billiard_sinai_periodic,
       collisiontime, propagate!, evolve!, construct

end # module

export AbstractParticle, Particle, MagneticParticle,
cyclotron
####################################################
## Particles
####################################################
"""
    AbstractParticle
Particle supertype.
"""
abstract type AbstractParticle{T<:AbstractFloat} end
eltype(p::AbstractParticle{T}) where {T} = T


"""
```julia
Particle(ic::Vector{T}) #where ic = [x0, y0, φ0]
Particle(x0, y0, φ0)
```
Create a particle with initial conditions `x0, y0, φ0`. It propagates as
a straight line.

The field `current_cell` shows at which cell of a periodic billiard is the particle
currently located.
"""
mutable struct Particle{T<:AbstractFloat} <: AbstractParticle{T}
    pos::SVector{2,T}
    vel::SVector{2,T}
    current_cell::SVector{2,T}
    function Particle(
        pos::SVector{2,T}, vel::SVector{2,T}, cc::SVector{2,T}) where{T<:AbstractFloat}
        new{T}(pos, normalize(vel), cc)
    end
end

function Particle(ic::AbstractVector{S}) where {S<:Real}
    T = S<:Integer ? Float64 : S
    φ0 = ic[3]
    pos = SVector{2,T}(ic[1:2]); vel = SVector{2,T}(cos(φ0), sin(φ0))
    return Particle(pos, vel, SVector{2,T}(0,0))
end
Particle(x::Real, y::Real, φ::Real) = Particle(collect(promote(x,y,φ)))
Particle() = Particle(rand(), rand(), rand()*2π)
show(io::IO, p::Particle{T}) where {T} =
print(io, "Particle{$T}\n",
"position: $(p.pos+p.current_cell)\nvelocity: $(p.vel)")



"""
```julia
MagneticParticle(ic::AbstractVector{T}, ω::Real) #where ic = [x0, y0, φ0]
MagneticParticle(x0, y0, φ0, ω)
```
Create a *magnetic* particle with initial conditions `x0, y0, φ0` and angular
velocity `ω`. It propagates as a circle instead of a line.

The field `current_cell` shows at which cell of a periodic billiard is the particle
currently located.
"""
mutable struct MagneticParticle{T<:AbstractFloat} <: AbstractParticle{T}
    pos::SVector{2,T}
    vel::SVector{2,T}
    current_cell::SVector{2,T}
    omega::T
    function MagneticParticle(pos::SVector{2,T}, vel::SVector{2,T},
        current_cell::SVector{2,T}, ω::T) where {T<:AbstractFloat}
        if ω==0
            throw(ArgumentError("Angular velocity of magnetic particle cannot be 0."))
        end
        new{T}(pos, normalize(vel), current_cell, ω)
    end
end

function MagneticParticle(ic::AbstractVector{T}, ω::Real) where {T<:Real}
    φ0 = ic[3]
    S = T<:Integer ? Float64 : T
    pos = SVector{2,S}(ic[1:2]); vel = SVector{2,S}(cos(φ0), sin(φ0))
    return MagneticParticle(pos, vel, SVector{2,S}(0,0), convert(S,ω))
end
function MagneticParticle(x0::Real, y0::Real, φ0::Real, ω::Real)
    a = collect(promote(x0, y0, φ0, ω))
    MagneticParticle(a[1:3], a[4])
end
MagneticParticle() = MagneticParticle([rand(), rand(), rand()*2π], 1.0)

show(io::IO, p::MagneticParticle{T}) where {T} =
print(io, "MagneticParticle{$T}\n",
"position: $(p.pos+p.current_cell)\nvelocity: $(p.vel)\nang. velocity: $(p.omega)")


"""
    cyclotron(p::MagneticParticle, use_cell = false)
Return center and radius of circular motion performed by the particle based on
`p.pos` (or `p.pos + p.current_cell`) and `p.vel`.
"""
@inline function cyclotron(p::MagneticParticle{T}, use_cell = false) where {T}
    ω = p.omega; r = 1/ω
    pos = use_cell ? p.pos + p.current_cell : p.pos
    c::SV{T} = pos - r*SV{T}(p.vel[2], -p.vel[1])
    return c, abs(r)
end

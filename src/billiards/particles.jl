export AbstractParticle, Particle, MagneticParticle
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
Particle(pos::SVector, vel::SVector)
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

Base.copy(p::Particle) =
Particle(p.pos, p.vel, p.current_cell)

function Particle(ic::AbstractVector{S}) where {S<:Real}
    T = S<:Integer ? Float64 : S
    φ0 = ic[3]
    pos = SVector{2,T}(ic[1:2]); vel = SVector{2,T}(cossin(φ0)...)
    return Particle(pos, vel, SVector{2,T}(0,0))
end
Particle(x::Real, y::Real, φ::Real) = Particle(collect(promote(x,y,φ)))
Particle() = Particle(rand(), rand(), rand()*2π)
function Particle(pos::SV{T}, vel::SV{T}) where {T}
    S = T<:Integer ? Float64 : T
    return Particle(pos, vel, SVector{2,S}(0.0, 0.0))
end
show(io::IO, p::Particle{T}) where {T} =
print(io, "Particle{$T}\n",
"position: $(p.pos+p.current_cell)\nvelocity: $(p.vel)")



"""
```julia
MagneticParticle(ic::AbstractVector{T}, ω::Real) # where ic = [x0, y0, φ0]
MagneticParticle(x0, y0, φ0, ω)
MagneticParticle(pos::SVector, vel::SVector, ω)
MagneticParticle(p::AbstractParticle, ω)
```
Create a *magnetic* particle with initial conditions `x0, y0, φ0` and angular
velocity `ω`. It propagates as a circle instead of a line,
with radius `1/abs(ω)`.

The field `current_cell` shows at which cell of a periodic billiard is the particle
currently located.
"""
mutable struct MagneticParticle{T<:AbstractFloat} <: AbstractParticle{T}
    pos::SVector{2,T}
    vel::SVector{2,T}
    current_cell::SVector{2,T}
    omega::T
    r::T
    center::SVector{2, T}
    function MagneticParticle(pos::SVector{2,T}, vel::SVector{2,T},
        current_cell::SVector{2,T}, ω::T) where {T<:AbstractFloat}
        if ω==0
            throw(ArgumentError("Angular velocity of magnetic particle cannot be 0."))
        end
        r = 1/ω
        c::SV{T} = pos - r*SV{T}(vel[2], -vel[1])
        new{T}(pos, normalize(vel), current_cell, ω, abs(1/ω), c)
    end
end

function Base.getproperty(p::MagneticParticle, s::Symbol)
    if s == :ω
        return Base.getfield(p, :omega)
    else
        return Base.getfield(p, s)
    end
end

Base.copy(p::MagneticParticle{T}) where {T} =
MagneticParticle(p.pos, p.vel, p.current_cell, p.omega)

function MagneticParticle(ic::AbstractVector{T}, ω::Real) where {T<:Real}
    φ0 = ic[3]
    S = T<:Integer ? Float64 : T
    pos = SVector{2,S}(ic[1:2]); vel = SVector{2,S}(cossin(φ0)...)
    return MagneticParticle(pos, vel, SVector{2,S}(0,0), convert(S,ω))
end
function MagneticParticle(x0::Real, y0::Real, φ0::Real, ω::Real)
    a = collect(promote(x0, y0, φ0, ω))
    MagneticParticle(a[1:3], a[4])
end
MagneticParticle() = MagneticParticle([rand(), rand(), rand()*2π], 1.0)
function MagneticParticle(pos::SV{T}, vel::SV{T}, ω::T) where {T}
    S = T<:Integer ? Float64 : T
    return MagneticParticle(pos, vel, SVector{2,S}(0,0), convert(S,ω))
end

MagneticParticle(p::AbstractParticle, ω) = MagneticParticle(p.pos, p.vel, p.current_cell, ω)

show(io::IO, p::MagneticParticle{T}) where {T} =
print(io, "MagneticParticle{$T}\n",
"position: $(p.pos+p.current_cell)\nvelocity: $(p.vel)\nang. velocity: $(p.omega)")


"""
    find_cyclotron(p::MagneticParticle)
Return the center of cyclotron motion of the particle.
"""
@inline find_cyclotron(p::MagneticParticle{T}) where {T} =
    c::SV{T} = p.pos - (1/p.omega)*SV{T}(p.vel[2], -p.vel[1])

@inline cyclotron(p) = p.center, p.r
function cyclotron(p, use_cell)
    pos = use_cell ? p.pos + p.current_cell : p.pos
    return pos, p.r
end

"""
Say that Everything is a parametric type to allow for Float or BigFloat (but not Integers!).
So be careful when defining your own methods to make sure all functions are parametric!
"""
module DynamicalBilliards

include("ParticlesObstacles.jl")
include("Propagation.jl")
include("StandardBilliards.jl")
include("PlotBilliards.jl")

export Particle, MagneticParticle, cyclotron, Disk, Circle, Antidot, FiniteWall, Obstacle,
       NullWall, PeriodicWall, resolvecollision!, randominside,
       billiard_rectangle, billiard_sinai, billiard_sinai_periodic,
       collisiontime, propagate!, evolve!, construct

end

#= comments about constructors
this is not correct:
immutable FiniteWall{T<:AbstractFloat} <: Wall{T}
  sp::SVector{2,T}
  ep::SVector{2,T}
  normal::SVector{2,T}
  FiniteWall(sp::SVector{2,T}, ep::SVector{2,T}, normal::SVector{2,T}) = begin
    d = dot(normal, ep-sp)
    if d != zero(T)
      error("Normal vector is not actually normal to the wall")
    end
    new(sp, ep, normal)
  end
end
function FiniteWall{T<:AbstractFloat}(sp::Vector{T}, ep::Vector{T}, n::Vector{T})
  FiniteWall(SVector{2, T}(sp), SVector{2, T}(ep), SVector{2, T}(n))
end

I think I do not recognize how "inner constructors work",
because when I try to instatize it:

sp = [1.0,0.0]; ep = [1.0, 1.0]; n =[-1.0,0.0]
w = FiniteWall(sp, ep, n)

I get

MethodError: no method matching FiniteWall{T<:AbstractFloat}
(::StaticArrays.SVector{2,Float64}, ::StaticArrays.SVector{2,Float64},
::StaticArrays.SVector{2,Float64})
Closest candidates are:
  FiniteWall{T<:AbstractFloat}{T}(::Any) at sysimg.jl:53
 in FiniteWall{T<:AbstractFloat}(::Array{Float64,1}, ::Array{Float64,1},
 ::Array{Float64,1}) at .\console:15

 What I have to do is add a {T} at the outer constructor call:
 function FiniteWall{T<:AbstractFloat}(sp::Vector{T}, ep::Vector{T}, n::Vector{T})
  FiniteWall{T}(SVector{2, T}(sp), SVector{2, T}(ep), SVector{2, T}(n)) #note the added {T}
end

SINCE ALL FALLS BACK TO CALLING AN TYPE, IT MUST BE SPECIFIED AS TYPE{} SOMEWHERE!

Also. I should NOT add {T} to the Inner constructor. Why?
We break an inner constructor into two parts.

    The type FiniteWall{T}
    the function FiniteWall which is NOT Parametric type function!

When invoking an inner constructor, you are not yourself calling the function FiniteWall.
You are "calling" the type FiniteWall{T}, which then calls the function FiniteWall,
after consuming and subsituting in the T type parameter everywhere.

Notice that I had to add this {T} because I was calling a TYPE, not a FUNCTION!
If want to call a function, I MUST NOT put the {T}! E.g.:
function resolvecollision!{T<:AbstractFloat}(p::Particle{T}, w::FiniteWall{T})
  n = normalvec(w, p.pos)
  i = p.vel
  p.vel = 2*dot(n, i)*n - i
end
Here the normalvec is a parametric type function:
normalvec{T<:AbstractFloat}(wall::Wall{T}, pos) = wall.normal
However, I have to call it without the {T} !!!

In general it is easy to understand:
The Parametric Type syntax {T,N} can only be:
After a Type name, no matter where that is located (called or defined)
After a function name, ONLY at the definition point.
You cannot call a function f{T}(3.5) for example.
=#

#= using BigFloats:
just type: Particle{BigFloat}([1,2,3]) and everything will automatically
be done as a BigFloat.
Why?
You can think of the annotations ::T in the type definition as having
implicit convert(T, pos) attached (when you use the provided constructor). Example:


type Par{T<:AbstractFloat}
  pos::SVector{2,T}
  vel::SVector{2,T}
end

function Par{T<:AbstractFloat}(x::Vector{T},y::Vector{T})
  p = SVector{2,T}(x...)
  v = SVector{2,T}(y...)
  Par{T}(x,y)
end

Par{BigFloat}([3.5,2.5],[2.5,3.5])

Par{BigFloat}(BigFloat[3.500000000000000000000000000000000000 it works anyway
=#

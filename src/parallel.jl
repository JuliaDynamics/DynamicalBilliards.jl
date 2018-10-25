using Distributed
export parallelize

parallelize(f, bd, t, n::Int) = parallelize(f, bd, t, [randominside(bd) for i in 1:n])
parallelize(f, bd, t, n::Int, ω) = parallelize(f, bd, t, [randominside(bd, ω) for i in 1:n])



"""
    parallelize(f, bd::Billiard, t, particles; partype = :threads)
Parallelize function `f` across the available particles. The parallelization type can
be `:threads` or `:pmap`, which use threads or a worker pool initialized with `addprocs`
_before_ `using DynamicalBilliards`.

`particles` can be:
* A `Vector` of particles.
* An integer `n` optionally followed by an angular velocity `ω`.
  This uses [`randominside`](@ref).

The functions usable here are:
* [`meancollisiontime`](@ref)
* [`escapetime`](@ref)
* [`lyapunovspectrum`](@ref) (returns only the maximal exponents)
* [`boundarymap`](@ref) (returns vector of vectors of 2-vectors _and_ `arcintervals`)
"""
function parallelize(f, bd::Billiard, t, particles::Vector{<:AbstractParticle};
    partype = :threads)
    if partype == :threads
        return threads_pl(f, bd, t, particles)
    elseif partype == :pmap
        return pmap_pl(f, bd, t, particles)
    end
end

function threads_pl(f, bd, t, particles)
    ret = _retinit(f, particles)
    Threads.@threads for i in 1:length(particles)
        @inbounds ret[i] = _getval(f, particles[i], bd, t)
    end
    return ret
end

function pmap_pl(f, bd, t, particles)
    g(p) = _getval(f, p, bd, t)
    ret = pmap(g, particles)
    return ret
end

_retinit(f, p::Vector{<:AbstractParticle{T}}) where {T} = zeros(T, length(p))
_retinit(::typeof(boundarymap), p::Vector{<:AbstractParticle{T}}) where {T} =
    Vector{Vector{SV{T}}}(undef, length(p))

_getval(f, p, bd, t) = f(p, bd, t)
_getval(f::typeof(lyapunovspectrum), p, bd, t) = @inbounds f(p, bd, t)[1]
_getval(f::typeof(lyapunovspectrum!), p, bd, t) = @inbounds f(p, bd, t)[1]

# Methods for boundary map are trickier because of the weird call signature
# and return signature
function threads_pl(f::typeof(boundarymap), bd, t, particles)
    intervals = arcintervals(bd)
    ret = _retinit(f, particles)
    Threads.@threads for i in 1:length(particles)
        @inbounds ret[i] = f(particles[i], bd, t, intervals)[1]
    end
    return ret, intervals
end
function pmap_pl(f::typeof(boundarymap), bd, t, particles)
    intervals = arcintervals(bd)
    g(p) = f(p, bd, t, intervals)
    ret = pmap(g, particles)
    return ret, intervals
end

using Distributed
export parallelize

parallelize(f, bd, t, n::Int) = parallelize(f, bd, t, [randominside(bd) for i in 1:n])
parallelize(f, bd, t, n::Int, ω) = parallelize(f, bd, t, [randominside(bd, ω) for i in 1:n])

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
_getval(f, p, bd, t) = f(p, bd, t)
_getval(f::typeof(lyapunovspectrum), p, bd, t) = f(p, bd, t)[1]

export boundarymap, boundarymap_portion
export psos!, psos

#######################################################################################
## Boundary Map
#######################################################################################

#this function only exists because incidence_angle from raysplitting.jl only works
#if you pass the particle *before* collision, which I cannot do because of bounce!
function reflection_angle(p::AbstractParticle{T}, a::Obstacle{T})::T where {T}
    n = normalvec(a, p.pos)
    inverse_dot = clamp(dot(p.vel, n), -1.0, 1.0)
    φ = acos(inverse_dot)
    if cross2D(p.vel, n) < 0
        φ *= -1
    end
    return φ
end



"""
    arcintervals(bt::Billiard)
Generate an array of `SVector`s, with the `i`th `SVector` containing the arc
length intervals corresponding to the `i`th `Obstacle` in `bt`.

Used by [`boundarymap`](@ref) to compute arc lengths.
"""
function arcintervals(bt::Billiard{T, D}) where {T, D}
    intervals = Vector{SVector{2,T}}(undef, D)
    current = zero(T)
    for i ∈ 1:D
        l = totallength(bt[i]) + current
        intervals[i] = SVector{2,T}(current, l)
        current = l
    end
    return intervals
end



"""
    boundarymap(bt::Billiard, t, particles)
Compute the boundary map of the
billiard `bt` by evolving each particle for total amount `t` (either float for
time or integer for collision number). See below for the returned values.

`particles` can be:
* A single particle.
* A `Vector` of particles.
* An integer `n` optionally followed by an angular velocity `ω`.

In the last case `n` random particles are produced in the billiard table using
[`randominside`](@ref). If `ω` is also given, then the particles are magnetic.

The measurement direction of the arclengths of the individual obstacles
is dictated by `bt`, see [`Billiard`](@ref) for more.

Return
* the arclengths at the collisions `ξs`
* the sines of the incidence angles at the collisions `sφs`
* obstacle arclength `intervals`

If `particles` is not a single particle then both `ξs` and `sφs` are vectors
of `Vector`. The `i` inner vector corresponds to the results of the
`i` initial condition/particle.

*Notice* - this function only works for normal specular reflection. Random reflections
or ray-splitting will give unexpected results.
"""
function boundarymap(bt::Billiard{T}, t,
                     ps::Vector{<:AbstractParticle{T}}) where {T}

    params = Vector{T}[]
    sines = Vector{T}[]
    intervals = arcintervals(bt)
    for p ∈ ps
        pparams, psines = boundarymap(bt, t, p, intervals)
        push!(params, pparams)
        push!(sines, psines)
    end
    return params, sines, intervals
end

function boundarymap(bt::Billiard{T}, t, par::AbstractParticle{T},
                     intervals = arcintervals(bt)) where {T}

    p = deepcopy(par)
    pparams = T[]
    pangles = T[]
    count = zero(T)
    t_to_write = zero(T)

    while count < t
        i, tmin = bounce!(p,bt)
        t_to_write += tmin

        if typeof(bt[i]) <: PeriodicWall
            continue # do not write output if collision with with PeriodicWall
        else
            push!(pparams, arclength(p, bt[i]) + intervals[i][1])
            push!(pangles, reflection_angle(p, bt[i]))
            # set counter
            count += increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end
    end #time, or collision number, loop
    return pparams, sin.(pangles), intervals
end


boundarymap(bt::Billiard, t, n::Int) =
    boundarymap(bt, t, [randominside(bt) for i ∈ 1:n])

boundarymap(bt::Billiard, t, n::Int, ω::AbstractFloat) =
    boundarymap(bt, t, [randominside(bt, ω) for i ∈ 1:n])



function boundarymap_portion(bt::Billiard{T}, t, par::AbstractParticle{T}, δξ, δφ = δξ) where {T}
    p = deepcopy(par)

    count = zero(T)
    t_to_write = zero(T)

    d = Dict{SV{Int}, Int}()
    intervals = arcintervals(bt)
    while count < t
        i, tmin = bounce!(p,bt)
        t_to_write += tmin

        if typeof(bt[i]) <: PeriodicWall
            continue # do not write output if collision with with PeriodicWall
        else
            # get Birkhoff coordinates
            ξ = arclength(p, bt[i]) + intervals[i][1]
            sφ = sin(reflection_angle(p, bt[i]))

            # compute index & increment dictionary entry
            ind = SV{Int}(floor(Int, ξ/δξ), floor(Int, (sφ + 1)/δφ))

            d[ind] = get(d, ind, 0) + 1

            # set counter
            count += increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end
    end #time or collision number loop

    #calculate ratio of visited boxes
    total_boxes = ceil(Int, totallength(bt)/δξ) * ceil(Int, 2/δφ)
    ratio = length(keys(d))/total_boxes

    return ratio, d
end




######################################################################################
# Poincare sos
######################################################################################
propagate_posvel(pos, p::Particle{T}, t) where {T} =
    (SV{T}(pos[1] + p.vel[1]*t, pos[2] + p.vel[2]*t), p.vel)

function propagate_posvel(pos, p::MagneticParticle{T}, t) where {T}
    ω = p.omega
    φ0 = atan2(p.vel[2], p.vel[1])
    ppos = SV{T}(sin(ω*t + φ0)/ω - sin(φ0)/ω, -cos(ω*t + φ0)/ω + cos(φ0)/ω)
    psos_vel = SV{T}(cos(ω*t + φ0), sin(ω*t + φ0))
    return pos + ppos, psos_vel
end


"""
    psos(bt::Billiard, plane::InfiniteWall, t, particles)
Compute the Poincaré section of the `particles` with the given `plane`, by evolving
each one for time `t` (either integer or float) inside `bt`.

The `plane` can be an [`InfiniteWall`](@ref) of *any* orientation, however only
crossings of the `plane` such that `dot(velocity, normal) < 0` are allowed, with
`normal` the normal unit vector of the `plane`.

`particles` can be:
* A single particle.
* A `Vector` of particles.
* An integer `n` optionally followed by an angular velocity `ω`.

Return the positions `poss` and velocities `vels` at the instances of crossing
the `plane`. If given more than one particle, the result is a vector of vectors
of vectors.
"""
function psos(
    bt::Billiard{T}, plane::InfiniteWall{T}, t, par::AbstractParticle{T}) where {T}

    p = deepcopy(par)
    rpos = SV{T}[]
    rvel = SV{T}[]
    count = zero(t)

    # periodic = isperiodic(bt)
    # ismagnetic = typeof(p) <: MagneticParticle
    # ismagnetic && (absω = abs(p.omega))
    # t_to_write = zero(T)

    while count < t
        # compute collision times
        tmin::T, i::Int = next_collision(p, bt)

        # if periodic && typeof(bt[i]) <: PeriodicWall

        tplane = collisiontime(p, plane)

        # if tplane is smaller, I intersect the section
        if tplane ≥ 0 && tplane < tmin && tplane != Inf
            psos_pos, psos_vel = propagate_posvel(p.pos, p, tplane)
            if dot(psos_vel, plane.normal) < 0 # only write crossings with 1 direction
                push!(rpos, psos_pos); push!(rvel, psos_vel)
            end
        end

        tmin == Inf && break
        # Now "bounce" the particle normally:
        tmin = relocate!(p, bt[i], tmin)
        resolvecollision!(p, bt[i])

        count += increment_counter(t, tmin)
    end
    return rpos, rvel
end

function psos(bt::Billiard{T}, plane::InfiniteWall{T}, t,
    pars::Vector{<:AbstractParticle{T}}) where {T}

    retpos = Vector{SV{T}}[]
    retvel = Vector{SV{T}}[]

    for p ∈ pars
        rpos, rvel = psos(bt, plane, t, p)
        push!(retpos, rpos); push!(retvel, rvel)
    end
    return retpos, retvel
end

psos(bt, plane, t, n::Integer) = psos(bt, plane, t, [randominside(bt) for i=1:n])
psos(bt, plane, t, n::Integer, ω::Real) =
psos(bt, plane, t, [randominside(bt, ω) for i=1:n])

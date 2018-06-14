export boundarymap, boundarymap_portion, phasespace_portion
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
    arcintervals(bd::Billiard)
Generate an array of `SVector`s, with the `i`th `SVector` containing the arc
length intervals corresponding to the `i`th `Obstacle` in `bd`.

Used by [`boundarymap`](@ref) to compute arc lengths.
"""
function arcintervals(bd::Billiard{T, D}) where {T, D}
    intervals = Vector{SVector{2,T}}(undef, D)
    current = zero(T)
    for i ∈ 1:D
        l = totallength(bd[i]) + current
        intervals[i] = SVector{2,T}(current, l)
        current = l
    end
    return intervals
end



"""
    boundarymap(bd::Billiard, t, particles)
Compute the boundary map of the
billiard `bd` by evolving each particle for total amount `t` (either float for
time or integer for collision number). See below for the returned values.

`particles` can be:
* A single particle.
* A `Vector` of particles.
* An integer `n` optionally followed by an angular velocity `ω`.

In the last case `n` random particles are produced in the billiard table using
[`randominside`](@ref). If `ω` is also given, then the particles are magnetic.

The measurement direction of the arclengths of the individual obstacles
is dictated by `bd`, see [`Billiard`](@ref) for more.

Return
* the arclengths at the collisions `ξs`
* the sines of the incidence angles at the collisions `sφs`
* obstacle arclength `intervals`

If `particles` is not a single particle then both `ξs` and `sφs` are vectors
of `Vector`. The `i` inner vector corresponds to the results of the
`i` initial condition/particle.

The returned values of this function are can be used in conjuction with the
function `plot_boundarymap` (requires `using PyPlot`) to plot the boundary map
in an intuitive way.

*Notice* - this function only works for normal specular reflection. Random reflections
or ray-splitting will give unexpected results.
"""
function boundarymap(bd::Billiard{T}, t,
                     ps::Vector{<:AbstractParticle{T}}) where {T}

    params = Vector{T}[]
    sines = Vector{T}[]
    intervals = arcintervals(bd)
    for p ∈ ps
        pparams, psines = boundarymap(bd, t, p, intervals)
        push!(params, pparams)
        push!(sines, psines)
    end
    return params, sines, intervals
end

function boundarymap(bd::Billiard{T}, t, par::AbstractParticle{T},
                     intervals = arcintervals(bd)) where {T}

    p = deepcopy(par)
    pparams = T[]
    pangles = T[]
    count = zero(T)
    t_to_write = zero(T)

    while count < t
        i, tmin = bounce!(p,bd)
        t_to_write += tmin

        if typeof(bd[i]) <: PeriodicWall
            continue # do not write output if collision with with PeriodicWall
        else
            push!(pparams, arclength(p, bd[i]) + intervals[i][1])
            push!(pangles, reflection_angle(p, bd[i]))
            # set counter
            count += increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end
    end #time, or collision number, loop
    return pparams, sin.(pangles), intervals
end


boundarymap(bd::Billiard, t, n::Int) =
    boundarymap(bd, t, [randominside(bd) for i ∈ 1:n])

boundarymap(bd::Billiard, t, n::Int, ω::AbstractFloat) =
    boundarymap(bd, t, [randominside(bd, ω) for i ∈ 1:n])


"""
    boundarymap_portion(bd::Billiard, t, p::AbstractParticle, δξ, δφ = δξ)
Calculate the portion of the boundary map of the billiard `bd` covered by the
particle `p` when it is evolved for time `t` (float or integer).

The boundary map is partitioned into boxes of size `(δξ, δφ)` and as the particle
evolves visited boxes are counted. The returned ratio is this count divided
by the total boxes of size `(δξ, δφ)` needed to cover the boundary map.

**Important:** This portion **does not** equate the portion the particle's orbit covers
on the full, three dimensional phase-space. Use the function
[`phasespace_portion`](@ref) for that!
"""
function boundarymap_portion(bd::Billiard{T}, t, par::AbstractParticle{T}, δξ, δφ = δξ) where {T}
    p = deepcopy(par)

    count = zero(T)
    t_to_write = zero(T)

    d = Dict{SV{Int}, Int}()
    intervals = arcintervals(bd)
    while count < t
        i, tmin = bounce!(p,bd)
        t_to_write += tmin

        if typeof(bd[i]) <: PeriodicWall
            continue # do not write output if collision with with PeriodicWall
        else
            # get Birkhoff coordinates
            ξ = arclength(p, bd[i]) + intervals[i][1]
            sφ = sin(reflection_angle(p, bd[i]))

            # compute index & increment dictionary entry
            ind = SV{Int}(floor(Int, ξ/δξ), floor(Int, (sφ + 1)/δφ))

            d[ind] = get(d, ind, 0) + 1

            # set counter
            count += increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end
    end #time or collision number loop

    #calculate ratio of visited boxes
    total_boxes = ceil(Int, totallength(bd)/δξ) * ceil(Int, 2/δφ)
    ratio = length(keys(d))/total_boxes

    return ratio, d
end

#######################################################################################
## Phase space ratio
#######################################################################################
function real_coordinates(ξ, sφ, bd::Billiard{T}; return_obst::Bool = false) where T
    abs(sφ) > 1 && throw(DomainError())#"|sin φ| must not be larger than 1"))
    lower = zero(T)
    upper = lower
    for (i, obst) ∈ enumerate(bd)
        #println("testing $(obst.name)")
        upper = lower + totallength(obst)
        #println("\tbounds: $lower:$upper")
        if ξ <= upper
            ret = real_coordinates(ξ - lower, sφ, obst)
            return return_obst ? (ret..., i) : ret
        end
        lower = upper
        #println("\tNEXT!")
    end
    #println("was: $ξ\tmax: $upper")
    throw(DomainError())#"ξ is too large for this billiard!"))
end

"""
    phasespace_portion(bd::Billiard, t, p::AbstractParticle, δξ, δφ = δξ)
Calculate the portion of the phase space of the billiard `bd` covered by the
particle `p` when it is evolved for time `t` (float or integer).

This function extends [`boundarymap_portion`](@ref) using a novel approach. For
each visited box of the boundary map, [`bounce!`](@ref) attributes a third dimension
(the collision time, equal to collision distance) which expands the two dimensions
of the boundary map to the three dimensions of the phase space.

The true phase space portion is then the weighted portion of boxes visited by the
particle, divided by the total weighted sum of boxes. The weights of the boxes are
the collision times.
"""
function phasespace_portion(bd::Billiard{T}, t, par::AbstractParticle{T}, δξ, δφ = δξ) where {T}
    PT = typeof(par)

    r, dict = boundarymap_portion(bd, t, par, δξ, δφ)

    maxξ = ceil(Int, totallength(bd)/δξ)
    maxφ = ceil(Int, 2/δφ)

    dummy = PT <: MagneticParticle ? MagneticParticle(zeros(T, 3)..., p.omega) : Particle(zeros(T, 3)...)

    total = zero(T)
    visited = zero(T)

    for ξcell ∈ 1:maxξ-1, φcell ∈ 1:maxφ-1

        #get center position
        ξc = (ξcell - 0.5)*δξ
        φc = (φcell - 0.5)*δφ - 1

        #convert to real space
        pos, vel, i = real_coordinates(ξc, φc, bd, return_obst=true)

        #set dummy coordinates
        dummy.pos = pos
        dummy.vel = vel
        specular!(dummy, bd[i])
        #  1c) bounce! & look what happens
        τ, = next_collision(dummy,bd)
        # 2. add to total
        total += τ
        (haskey(dict, SV{Int64}(ξcell, φcell))) && (visited += τ)
    end


    return visited/total
end



######################################################################################
# Poincare sos
######################################################################################
propagate_posvel(pos, p::Particle{T}, t) where {T} =
    (SV{T}(pos[1] + p.vel[1]*t, pos[2] + p.vel[2]*t), p.vel)

function propagate_posvel(pos, p::MagneticParticle{T}, t) where {T}
    ω = p.omega
    φ0 = atan2(p.vel[2], p.vel[1])
    s0, c0 = sincos(φ0)
    sωφ0, cωφ0 = sincos(ω*t + φ0)
    ppos = SV{T}(sωφ0/ω - s0/ω, -cωφ0/ω + c0/ω)
    psos_vel = SV{T}(cωφ0, sωφ0)
    return pos + ppos, psos_vel
end


"""
    psos(bd::Billiard, plane::InfiniteWall, t, particles)
Compute the Poincaré section of the `particles` with the given `plane`, by evolving
each one for time `t` (either integer or float) inside `bd`.

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
    bd::Billiard{T}, plane::InfiniteWall{T}, t, par::AbstractParticle{T}) where {T}

    p = deepcopy(par)
    rpos = SV{T}[]
    rvel = SV{T}[]
    count = zero(t)

    # periodic = isperiodic(bd)
    # ismagnetic = typeof(p) <: MagneticParticle
    # ismagnetic && (absω = abs(p.omega))
    # t_to_write = zero(T)

    while count < t
        # compute collision times
        tmin::T, i::Int = next_collision(p, bd)

        # if periodic && typeof(bd[i]) <: PeriodicWall

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
        tmin = relocate!(p, bd[i], tmin)
        resolvecollision!(p, bd[i])

        count += increment_counter(t, tmin)
    end
    return rpos, rvel
end

function psos(bd::Billiard{T}, plane::InfiniteWall{T}, t,
    pars::Vector{<:AbstractParticle{T}}) where {T}

    retpos = Vector{SV{T}}[]
    retvel = Vector{SV{T}}[]

    for p ∈ pars
        rpos, rvel = psos(bd, plane, t, p)
        push!(retpos, rpos); push!(retvel, rvel)
    end
    return retpos, retvel
end

psos(bd, plane, t, n::Integer) = psos(bd, plane, t, [randominside(bd) for i=1:n])
psos(bd, plane, t, n::Integer, ω::Real) =
psos(bd, plane, t, [randominside(bd, ω) for i=1:n])

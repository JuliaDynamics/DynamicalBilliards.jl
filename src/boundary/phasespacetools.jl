################################################################################
##
##  This file contains boxcounting methods to analyse phase space volume
##  ratios in 2D billiard maps and 3D billiard flows.
##
################################################################################
export boundarymap_portion, phasespace_portion

## BOUNDARY MAP PORTION ########################################################

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
function boundarymap_portion(bd::Billiard{T}, t,
                             par::AbstractParticle{T}, δξ, δφ = δξ;
                             intervals = arcintervals(bd)) where {T}
    
    p = deepcopy(par)

    count = zero(T)
    t_to_write = zero(T)

    d = Dict{SV{Int}, Int}()
   
    while count < t
        i, tmin = bounce!(p,bd)
        t_to_write += tmin

        if typeof(bd[i]) <: PeriodicWall
            continue # do not write output if collision with with PeriodicWall
        else
            # get Birkhoff coordinates
            ξ = to_bcoords(p, bd[i]) + intervals[i][1]
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


## PHASE SPACE PORTION #########################################################

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
function phasespace_portion(bd::Billiard{T}, t,
                            par::AbstractParticle{T}, δξ, δφ = δξ) where {T}

    PT = typeof(par)

    r, dict = boundarymap_portion(bd, t, par, δξ, δφ,
                                  intervals = arcintervals(bd))
    
    ints =  arcintervals(bd)
    
    maxξ = ceil(Int, totallength(bd)/δξ)
    maxφ = ceil(Int, 2/δφ)

    dummy = PT <: MagneticParticle ? MagneticParticle(zeros(T, 3)..., p.omega) :
        Particle(zeros(T, 3)...)

    total = zero(T)
    visited = zero(T)

    for ξcell ∈ 1:maxξ-1, φcell ∈ 1:maxφ-1

        #get center position
        ξc = (ξcell - 0.5)*δξ
        φc = (φcell - 0.5)*δφ - 1

        #convert to real space
        pos, vel, i = from_bcoords(ξc, φc, bd, return_obstacle=true,
                                       intervals = ints)

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


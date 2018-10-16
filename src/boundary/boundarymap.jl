export to_bcoords, from_bcoords, arcintervals, totallength
export boundarymap
#######################################################################################
## Arclengths
#######################################################################################
"""
    totallength(o::Obstacle)
Return the total boundary length of `o`.
"""
@inline totallength(o::Wall) = norm(o.ep - o.sp)
@inline totallength(o::Semicircle) = π*o.r
@inline totallength(o::Circular) = 2π*o.r
@inline totallength(e::Ellipse) = 4e.arc

@inline totallength(bd::Billiard) = sum(totallength(x) for x in bd.obstacles)

"""
    arcintervals(bd::Billiard) -> s
Generate a vector `s`, with entries being the delimiters of the
arclengths of the obstacles of the billiard.
The arclength from `s[i]` to `s[i+1]` is the arclength spanned by
the `i`th obstacle.

`s` is used to transform an arc-coordinate `ξ` from local to global and
vice-versa. A local `ξ` becomes global by adding `s[i]` (where `i` is the
index of current obstacle). A global `ξ` becomes local by subtracting `s[i]`.

See also [`boundarymap`](@ref), [`to_bcoords`](@ref), [`from_bcoords`](@ref).
"""
function arcintervals(bd::Billiard{T, D}) where {T, D}
    intervals = SVector{D+1,T}(0, map(x->totallength(x), bd.obstacles)...)
    return cumsum(intervals)
end

################################################################################
## Coordinate systems
################################################################################
"""
    to_bcoords(pos, vel, o::Obstacle) -> ξ, sφ
Convert the real coordinates `pos, vel` to
boundary coordinates (also known as Birkhoff coordinates)
`ξ, sφ`, assuming that `pos` is on the obstacle.

`ξ` is the arc-coordinate, i.e. it parameterizes the arclength of the
obstacle. `sφ` is the
sine of the angle between the velocity vector and the vector normal
to the obstacle.

The arc-coordinate `ξ` is measured as:
* the distance from start point to end point in `Wall`s
* the arc length measured counterclockwise from the open face in `Semicircle`s
* the arc length measured counterclockwise from the rightmost point
  in `Circular`/`Ellipse`s

Notice that this function returns the *local* arclength. To get the global
arclength parameterizing an entire billiard, simply do
`ξ += arcintervals(bd)[i]` if the index of obstacle `o` is `i`.

See also [`from_bcoords`](@ref), which is the inverse
function.
"""
to_bcoords(p::AbstractParticle, o::Obstacle) = to_bcoords(p.pos, p.vel, o)

function to_bcoords(pos::SV, vel::SV, o::Obstacle)
    n = normalvec(o, pos)
    sinφ = cross2D(vel, n)
    ξ = _ξ(pos, o)
    return ξ, sinφ
end

# Internal function ξ is simply returning the boundary "arc-coordinate"
_ξ(pos::SV, o::Wall) = norm(pos - o.sp)

function _ξ(pos::SV{T}, o::Semicircle{T}) where {T<:AbstractFloat}
    #project pos on open face
    chrd = SV{T}(-o.facedir[2],o.facedir[1])
    d = (pos - o.c)/o.r
    x = dot(d, chrd)
    r =  acos(clamp(x, -1, 1))*o.r
    return r
end

function _ξ(pos::SV{T}, o::Circular{T}) where {T<:AbstractFloat}
    d = (pos - o.c)/o.r
    if d[2] > 0
        r = acos(clamp(d[1], -1, 1))
    else
        r = π + acos(-clamp(d[1], -1, 1))
    end
    return r*o.r
end

function _ξ(pos::SV{T}, o::Ellipse{T}) where {T<:AbstractFloat}
    θ = atan(pos[2] - o.c[2], pos[1] - o.c[1]) + π
    return ellipse_arclength(θ, o)
end


"""
    from_bcoords(ξ, sφ, o::Obstacle) -> pos, vel
Convert the boundary coordinates `ξ, sφ` on the obstacle to
real coordinates `pos, vel`.

Note that `vel` always points away from the obstacle.

This function is the inverse of [`to_bcoords`](@ref).
"""
function from_bcoords(ξ::T, sφ::T, o::Obstacle{T}) where {T}
    pos = real_pos(ξ, o)
    cφ = sqrt(1-sφ^2) # = cos(asin(sφ))
    n = normalvec(o, pos)
    vel = SV{T}(n[1]*cφ + n[2]*sφ, -n[1]*sφ + n[2]*cφ)

    return pos, vel
end

"""
    real_pos(ξ, o::Obstacle)
Converts the arclength coordinate `ξ` relative to the obstacle `o` into a real space
position vector.
"""
real_pos(ξ, o::Wall) = o.sp + ξ*normalize(o.ep - o.sp)

function real_pos(ξ, o::Semicircle{T}) where T
    sξ, cξ = sincos(ξ/o.r)
    chrd = SV{T}(-o.facedir[2], o.facedir[1])
    return o.c - o.r*(sξ*o.facedir - cξ*chrd)
end

real_pos(ξ, o::Circular{T}) where T = o.c .+ o.r * SV{T}(cossin(ξ/o.r))

function real_pos(ξ, o::Ellipse{T}) where T
    error("`from_bcoords` is not implemented for Ellipse. If you know of a way "*
    "to find a point (2-vector) on an ellipse given an arclength ξ, please let us know!")
    return nothing
end


"""
    from_bcoords(ξ, sφ, bd::Billiard, intervals = arcintervals(bd))
Same as above, but now `ξ` is considered to be the global arclength,
parameterizing the entire billiard, instead of a single obstacle.
"""
function from_bcoords(ξ, sφ, bd::Billiard{T}, intervals = arcintervals(bd)) where T

    for (i, obst) ∈ enumerate(bd)
        if ξ <= intervals[i+1]
            pos, vel  = from_bcoords(ξ - intervals[i], sφ, obst)
            return pos, vel, i
        end
    end
    throw(DomainError(ξ ,"ξ is too large for this billiard!"))
end



#######################################################################################
## Boundary Map
#######################################################################################
"""
    boundarymap(p, bd, t [,intervals]) → bmap, arclengths
Compute the boundary map of the particle `p` in the billiard `bd` by evolving
the particle for total amount `t` (either float for
time or integer for collision number).

Return a vector of 2-vectors `bmap` and also `arclengths(bd)`. The first entry of each
element of `bmap` is the arc-coordinate at collisions ``\\xi``, while the second  is the
sine of incidence angle ``\\sin(\\phi_n)``.

The measurement direction of the arclengths of the individual obstacles
is dictated by their order in `bd`. The sine of the angle is computed
*after* specular reflection has taken place.

The returned values of this function can be used in conjuction with the
function [`plot_boundarymap`](@ref) (requires `using PyPlot`) to plot the boundary map
in an intuitive way.

*Notice* - this function only works for normal specular reflection. Random reflections
or ray-splitting will give unexpected results.

See also [`to_bcoords`](@ref), [`boundarymap_portion`](@ref).
See [`parallelize`](@ref) for a parallelized version.
"""
function boundarymap(par::AbstractParticle{T}, bd::Billiard{T}, t,
                     intervals = arcintervals(bd)) where {T}

    p = copy(par)
    bmap = SV{T}[]
    if typeof(t) == Int
        sizehint!(bmap, t)
    end
    count = zero(t); t_to_write = zero(T)

    while count < t
        i, tmin = bounce!(p,bd)
        t_to_write += tmin

        if isperiodic(i, bd)
            continue # do not write output if collision with with PeriodicWall
        else
            ξ, sφ = to_bcoords(p, bd[i])
            @inbounds push!(bmap, SV{T}(ξ + intervals[i], sφ))
            # set counter
            count += increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end
    end #time, or collision number, loop
    return bmap, intervals
end

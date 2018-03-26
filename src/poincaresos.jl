using StaticArrays

export evolve_boundary!, evolve_boundary, poincaresection

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


#######################################################################################
## Helper functions to rearrange output of evolve_boundary
#######################################################################################


#WARNING: This code is ugly.
"""
    function shiftconstruct(bt::BilliardTable)
Use the `sortorder` field (see [`BilliardTable`](@ref)
and [`poincaresection`](@ref)) and [`totallength`](@ref) to
generate an array of `SVector`s, with the `i`th `SVector` containing the arc
length intervals corresponding to the `i`th `Obstacle` in `bt`.

Used by [`poincaresection`](@ref) to compute arc lengths.
"""
function shiftconstruct(bt::BilliardTable{T}) where {T}
    len = length(bt.sortorder)
    intervals = Array{SVector{2,T}}(len)
    signs = Array{Int}(len)
    current = zero(T)
    for i ∈ bt.sortorder
        absi = abs(i)
        l = totallength(bt[absi]) + current
        intervals[absi] = SVector{2,T}(current, l)
        signs[absi] = sign(i)
        current = l
    end
    return intervals, signs
end

#######################################################################################
## Main poincaresection function
#######################################################################################


#TODO:rewrite docstring
"""
```julia
poincaresection(bt::BilliardTable, t, ps::Vector{<:AbstractParticle})
poincaresection(bt::BilliardTable, t, n::Int [, ω])
```
Compute the poincare section (also called boundary map) of the
billiard table `bt` by evolving each particle for total amount `t` (either float for
time or integer for collision number).

If `n::Int` is given instead of `ps`,
then `n` random particles are produced in the billiard table. If `ω` is also
given, then the particles are magnetic.

The sorting of the arclengths of the individual obstacles is dictated by
the `sortorder` field of `bt`, see [`BilliardTable`](@ref) for more.

# HERE GOES PROPER DESCRIPTION OF THE RETURNED STUFF AND / OR WHAT THE METHOD DOES

Returns
* a vector of the arclengths at the collisions
* a vector of incidence angles at the collisions
* an array of intervals corresponding to the obstacle arc lengths ???
"""
function poincaresection(bt::BilliardTable{T}, t,
                         ps::Vector{<:AbstractParticle{T}}) where {T}

    params = T[]
    angles = T[]

    intervals, signs = shiftconstruct(bt)

    for p ∈ ps
        count = zero(T)
        t_to_write = zero(T)

        while count < t
            i, tmin = bounce!(p,bt)
            t_to_write += tmin

            if typeof(bt[i]) <: PeriodicWall
                continue # do not write output if collision with with PeriodicWall
            else
                if signs[i] > 0
                    push!(params, arclength(p, bt[i]) + intervals[i][1])
                else
                    push!(params, intervals[i][2] - arclength(p, bt[i]))
                end
                push!(angles, reflection_angle(p, bt[i]))
                # set counter
                count += increment_counter(t, t_to_write)
                t_to_write = zero(T)
            end
        end #time, or collision number, loop
    end

    return params, angles, intervals
end

poincaresection(bt::BilliardTable, t, n::Int) =
    poincaresection(bt, t, [randominside(bt) for i ∈ 1:n])

poincaresection(bt::BilliardTable, t, n::Int, ω::AbstractFloat) =
    poincaresection(bt, t, [randominside(bt, ω) for i ∈ 1:n])

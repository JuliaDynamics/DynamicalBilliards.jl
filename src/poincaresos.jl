export poincaresection

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

Used by [`poincaresection`](@ref) to compute arc lengths.
"""
function arcintervals(bt::Billiard{T, D}) where {T, D}
    intervals = Vector{SVector{2,T}}(D)
    current = zero(T)
    for i ∈ 1:D
        l = totallength(bt[i]) + current
        intervals[i] = SVector{2,T}(current, l)
        current = l
    end
    return intervals
end



"""
```julia
poincaresection(bt::Billiard, t, ps::Vector{<:AbstractParticle})
poincaresection(bt::Billiard, t, n::Int [, ω])
```
Compute the poincare section (also called boundary map) of the
billiard table `bt` by evolving each particle for total amount `t` (either float for
time or integer for collision number). See below for the returned values.

If `n::Int` is given instead of `ps`,
then `n` random particles are produced in the billiard table. If `ω` is also
given, then the particles are magnetic.

The sorting and measurement direction of the arclengths of the individual obstacles
is dictated by the `sortorder` field of `bt`, see [`Billiard`](@ref) for more.

## Returns
* the arclengths at the collisions `ξs`
* the incidence angles at the collisions `φs`
* obstacle arclength `intervals`

Both `ξs` and `φs` are vectors of `Vector`.
The `i` inner vectors correspond to the results of the `i` initial condition/particle.

The `intervals` is a vector of `SVector`. The `i` entry of `intervals` is the
arclength spanned by the `i` obstacle of the billiard table. The direction
of the measurement of the arclength is dictated by `bt.sortorder`, see
[`Billiard`](@ref) for more.
"""
function poincaresection(bt::Billiard{T}, t,
                         ps::Vector{<:AbstractParticle{T}}) where {T}

    params = Vector{T}[]
    angles = Vector{T}[]

    intervals = arcintervals(bt)

    for p ∈ ps
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
                if !bt.inverted[i]
                    push!(pparams, arclength(p, bt[i]) + intervals[i][1])
                else
                    push!(pparams, intervals[i][2] - arclength(p, bt[i]))
                end
                push!(pangles, reflection_angle(p, bt[i]))
                # set counter
                count += increment_counter(t, t_to_write)
                t_to_write = zero(T)
            end
        end #time, or collision number, loop

        push!(params, pparams)
        push!(angles, pangles)
    end

    return params, angles, intervals
end

poincaresection(bt::Billiard, t, n::Int) =
    poincaresection(bt, t, [randominside(bt) for i ∈ 1:n])

poincaresection(bt::Billiard, t, n::Int, ω::AbstractFloat) =
    poincaresection(bt, t, [randominside(bt, ω) for i ∈ 1:n])

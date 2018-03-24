using StaticArrays

export evolve_boundary!, evolve_boundary, psos

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
## Main evolve-like functions
#######################################################################################


"""
```julia
function evolve_boundary!(p::Particle, bt::BilliardTable, t)
```
"""
function evolve_boundary!(p::Particle{T},  bt::BilliardTable{T}, t) where {T<:AbstractFloat}
    if t <= 0
        throw(ArgumentError("`evolve!()` cannot evolve backwards in time."))
    end

    rt = T[]
    rparam = T[]
    rangle = T[]
    rindex = Int[]

    count = zero(T)
    t_to_write = zero(T)

    while count < t
        i, tmin = bounce!(p,bt)
        t_to_write += tmin

        if typeof(bt[i]) <: PeriodicWall
            continue # do not write output if collision with with PeriodicWall
        else
            push!(rparam, arclength(p, bt[i]))
            push!(rangle, reflection_angle(p, bt[i]))
            push!(rindex, i)
            push!(rt, t_to_write)
            # set counter
            count += increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end
    end #time, or collision number, loop

  return rt, rindex, rparam, rangle
end


#non-mutating version
function evolve_boundary(p::Particle{T}, bt::BilliardTable{T}, t) where {T<:AbstractFloat}
    p2 = deepcopy(p)
    return evolve_boundary!(p2, bt, t)
end

#######################################################################################
## Helper functions to rearrange output of evolve_boundary
#######################################################################################


#WARNING: This code is ugly.

function shiftconstruct(bt::BilliardTable{T}, sortorder::Array{Int}) where {T}
    intervals = Array{SVector{2,T}}(length(sortorder))
    signs = Array{Int}(length(sortorder))
    current = zero(T)
    for i ∈ sortorder
        absi = abs(i)
        l = totallength(bt[absi]) + current
        intervals[absi] = SVector{2,T}(current, l)
        signs[absi] = sign(i)
        current = l
    end
    return intervals, signs
end


"""
```julia
function psos(ps::Vector{Particle}, bt::BilliardTable, t,
                     sortorder::Vector{Int})
function psos(n::Int, bt::BilliardTable, t, sortorder::Vector{Int})
```
This function calls [`evolve_boundary'](@ref) and rearranges its output using the
    information provided in `sortorder`. `sortorder`should contain the indices of
    the `Obstacle`s in `bt` in the correct order, with the signs signifying the
    sign with which the individual `arclength`s are supposed to added up.

If `n` is given instead of ps, generates `n` random particles inside bt and
    evolves them.

Returns
* an array of arc length parameters
* an array of incidence angles
* an array of intervals corresponding to the obstacle arc lengths

"""
function psos(ps::Vector{Particle{T}}, bt::BilliardTable{T}, t,
                     sortorder::Vector{Int}) where {T}
    plotξs = T[]
    plotφs = T[]

    intervals, signs = shiftconstruct(bt, sortorder)
    for p ∈ ps
        ts, index, params, angles = evolve_boundary(p, bt, t)

        for (j,i) ∈ enumerate(index)
            if signs[i] > 0
                push!(plotξs, params[j] + intervals[i][1])
            else
                push!(plotξs, intervals[i][2] - params[j])
            end
            push!(plotφs,angles[j])
        end
    end

    return plotξs, plotφs, intervals
end

function psos(n::Int, bt::BilliardTable, t, sortorder::Vector{Int})
    ps = [randominside(bt) for i ∈ 1:n]
    return psos(ps, bt, t, sortorder)
end

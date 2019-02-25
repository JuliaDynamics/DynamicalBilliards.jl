export lyapunovspectrum!, lyapunovspectrum
export perturbationgrowth!, perturbationgrowth, perturbationevolution

const δqind = SV{Int}(1,2)
const δpind = SV{Int}(3,4)

@inline curvature(::Semicircle{T}) where {T} = -one(T)
@inline curvature(::Disk{T}) where {T} = one(T)

################################################################################
## SPECULAR (LINEAR)
################################################################################

#="""
    specular!(p::AbstractParticle, o::Obstacle, offset::MArray)
Perform specular reflection based on the normal vector of the Obstacle.
The function updates the position and velocity of the particle
together with the components of 4 offset vectors stored in the matrix
`offset` as columns.
"""=#
function specular!(p::Particle{T}, o::Circular{T},
                   offset::Vector{SVector{4, T}}) where {T<:AbstractFloat}
    n = normalvec(o, p.pos)
    ti = SV{T}(-p.vel[2],p.vel[1])

    cosa = -dot(n, p.vel)
    p.vel = p.vel + 2.0*cosa*n

    tf = SV{T}(-p.vel[2], p.vel[1])

    @inbounds for k in eachindex(offset)
        δqprev = offset[k][δqind]
        δpprev = offset[k][δpind]
        # Formulas from Dellago, Posch and Hoover, PRE 53, 2, 1996: 1485-1501 (eq. 27)
        # with norm(p) = 1
        δq  = δqprev - 2*dot(δqprev,n)*n
        δp  = δpprev - 2*dot(δpprev,n)*n - curvature(o)*2/o.r*dot(δqprev,ti)/cosa*tf
        ###
        offset[k] = vcat(δq, δp)
    end
end

function specular!(p::Particle{T}, o::Union{InfiniteWall{T},FiniteWall{T}},
                   offset::Vector{SVector{4, T}}) where {T<:AbstractFloat}

    n = normalvec(o, p.pos)
    specular!(p, o)

    for k in eachindex(offset)
        δqprev = offset[k][δqind]
        δpprev = offset[k][δpind]
        # Formulas from Dellago, Posch and Hoover, PRE 53, 2, 1996: 1485-1501 (eq. 20)
        δq  = δqprev - 2.0*dot(δqprev,n)*n
        δp  = δpprev - 2.0*dot(δpprev,n)*n
        ###
        offset[k] = vcat(δq, δp)
    end

end

################################################################################
## SPECULAR (MAGNETIC)
################################################################################

function specular!(p::MagneticParticle{T}, o::Circular{T},
                   offset::Vector{SVector{4, T}}) where {T<:AbstractFloat}
    n = normalvec(o, p.pos)

    cosa = -dot(n, p.vel)

    Rn = SV{T}(-n[2], n[1])
    ti = SV{T}(-p.vel[2],p.vel[1])

    # magterm = ω*δτ_c*(BR-RB)*p_i / δqprev_normal
    # i.e. everything that can be calculated without using the previous offset vector
    magterm = -2*p.omega*((dot(p.vel, Rn)/(-cosa))*n + Rn)

    #actual specular reflection should not occur before magterm is computed
    p.vel = p.vel + 2*cosa*n
    tf = SV{T}(-p.vel[2], p.vel[1])


    for k in eachindex(offset)
        δqprev = offset[k][δqind]
        δpprev = offset[k][δpind]

        δqprev_normal = dot(δqprev, n)

        # Formulas derived analogously to Dellago et al., PRE 53, 2, 1996: 1485-1501
        δq  = δqprev - 2*δqprev_normal*n
        δp  = δpprev - 2*dot(δpprev,n)*n - curvature(o)* 2/o.r * dot(δqprev,ti)/cosa*tf + magterm*δqprev_normal
        ###
        offset[k] = vcat(δq, δp)
    end

end

function specular!(p::MagneticParticle{T}, o::Union{InfiniteWall{T},FiniteWall{T}},
                   offset::Vector{SVector{4, T}}) where {T<:AbstractFloat}

    n = normalvec(o, p.pos)

    cosa = -dot(n, p.vel)
    Rn = SV{T}(-n[2], n[1])

    # magterm = ω*δτ_c*(BR-RB)*p_i / δqprev_normal
    # i.e. everything that can be calculated without using the previous offset vector
    magterm = -2.0*p.omega*((dot(p.vel, Rn)/(-cosa))*n + Rn)

    #actual specular reflection should not occur before magterm is calculated
    p.vel = p.vel + 2*cosa*n

    for k in eachindex(offset)
        δqprev = offset[k][δqind]
        δpprev = offset[k][δpind]

        δqprev_normal = dot(δqprev, n)

        # Formulas derived analogously to Dellago et al., PRE 53, 2, 1996: 1485-1501
        δq  = δqprev - 2.0*δqprev_normal*n
        δp  = δpprev - 2.0*dot(δpprev,n)*n + magterm*δqprev_normal
        ###
        offset[k] = vcat(δq, δp)
    end
end

################################################################################
## RESOLVECOLLISON
################################################################################

#="""
    resolvecollision!(p::AbstractParticle, o::Union{Disk, InfiniteWall}, offset::MArray)
Resolve the collision between particle `p` and obstacle `o` of type *Circular*,
updating the components of the offset vectors stored in the matrix `offset` as columns.
"""=#
resolvecollision!(p::Particle{T}, o::Obstacle{T},
                  offset::Vector{SVector{4, T}}) where {T} = specular!(p, o, offset)

resolvecollision!(p::Particle{T}, o::PeriodicWall{T},
                  offset::Vector{SVector{4, T}}) where {T} = resolvecollision!(p, o)


function resolvecollision!(p::MagneticParticle{T}, o::Obstacle{T},
    offset::Vector{SVector{4, T}}) where {T}

    specular!(p, o, offset)
    p.center = find_cyclotron(p)
    return
end

# this method only exists to avoid ambiguity for MagneticParticles colliding with
# periodic walls.
resolvecollision!(p::MagneticParticle{T}, o::PeriodicWall{T},
                  offset::Vector{SVector{4, T}}) where {T} = resolvecollision!(p, o)


################################################################################
## OFFSET PROPAGATION
################################################################################

"""
    propagate_offset!(offset::MArray{Tuple{4,4},T}, p::AbstractParticle)
Computes the linearized evolution of the offset vectors during propagation for a
time interval `t`
"""
function propagate_offset!(offset::Vector{SVector{4, T}}, t::T,     #linear case
                           p::Particle{T}) where T
    for k in eachindex(offset)
        δΓ = offset[k]
        offset[k] = SVector{4,T}(δΓ[1] + t*δΓ[3], δΓ[2] + t*δΓ[4], δΓ[3], δΓ[4])
    end
end

#magnetic
function propagate_offset!(offset::Vector{SVector{4, T}}, t::T,
                           p::MagneticParticle{T}) where T
    ω = p.omega
    sω, cω = sincos(ω*t)

    for k in eachindex(offset)
        δΓ = offset[k]
        offset[k] = SVector{4,T}(δΓ[1] + sω/ω*δΓ[3] + (cω - 1)/ω*δΓ[4],
                                 δΓ[2] + (1 - cω)/ω*δΓ[3] + sω/ω*δΓ[4],
                                 cω*δΓ[3] - sω*δΓ[4],
                                 sω*δΓ[3] + cω*δΓ[4])
    end
end

################################################################################
## PROPAGATION & RELOCATION
################################################################################

#="""
    propagate!(p::AbstractParticle{T}, newpos::SV{T}, t::T,
    offset::MArray{Tuple{4,4},T})
Propagate the particle `p` for given time `t`, changing appropriately the the
`p.pos` and `p.vel` fields together with the components of the offset vectors
stored in the `offset` matrix.
"""=#
function propagate!(p::AbstractParticle{T}, newpos::SV{T}, t::T,
                    offset::Vector{SVector{4, T}}) where {T<: AbstractFloat}

    propagate!(p, newpos, t)
    propagate_offset!(offset, t, p)
    return
end


#="""
    relocate(p::AbstractParticle, o::Obstacle, t, cp::SV{T}, offset::MArray)
Propagate the particle's position for time `t` (corrected) and update the components
of the `offset` matrix.
"""=#
function relocate!(p::AbstractParticle{T}, o::Obstacle{T}, tmin, cp::SV{T},
                   offset::Vector{SVector{4, T}}) where {T <: AbstractFloat}

    okay = relocate!(p, o, tmin, cp)
    propagate_offset!(offset, tmin, p)
    return okay
end


################################################################################
## HIGH-LEVEL FUNCTION
################################################################################
function lyapunovspectrum!(p::AbstractParticle{T}, bd::Billiard{T}, t;
                           warning::Bool = false) where {T<:AbstractFloat}

    if t ≤ 0.0
        throw(ArgumentError(
            "`lyapunovspectrum()` cannot evolve backwards in time."))
    end

    # intial offset vectors
    offset = [SVector{4, T}(1,0,0,0), SVector{4, T}(0,1,0,0),
              SVector{4, T}(0,0,1,0), SVector{4, T}(0,0,0,1)]
    λ = zeros(T, 4)

    timecount = zero(T)
    count = zero(t)

    # check for pinning before evolution
    if ispinned(p, bd)
        warning && @warn "Pinned particle!"
        return λ
    end

    ismagnetic = typeof(p) <: MagneticParticle

    while count < t
        # bouncing
        i::Int, tmin::T, cp::SV{T} = next_collision(p, bd)
        relocate!(p, bd[i], tmin, cp, offset)
        resolvecollision!(p, bd[i], offset)

        timecount += tmin
        count += increment_counter(t, tmin)

        # update cyclotron data
        ismagnetic && (p.center = find_cyclotron(p))

        # QR decomposition to get Lyapunov spectrum
        Q, R = qr(hcat(offset[1], offset[2], offset[3], offset[4]))
        offset[1], offset[2], offset[3], offset[4] = Q[:, 1], Q[:, 2], Q[:, 3], Q[:, 4]

        for i ∈ 1:4
            λ[i] += log(abs(R[i,i]))
        end

    end#time loop

    return λ./timecount
end

"""
    lyapunovspectrum([p::AbstractParticle,] bd::Billiard, t)
Returns the finite time lyapunov exponents (averaged over time `t`)
for a given particle in a billiard table using the method outlined in [1].
`t` can be either `Float` or `Int`, meaning total time or total amount
of collisions.

Returns zeros for pinned particles.

If a particle is not given, a random one is picked through [`randominside`](@ref).
See [`parallelize`](@ref) for a parallelized version.

[1] : Ch. Dellago *et al*, [Phys. Rev. E **53** (1996)](http://link.aps.org/doi/10.1103/PhysRevE.53.1485)
"""
lyapunovspectrum(p::AbstractParticle, args...) = lyapunovspectrum!(copy(p), args...)
lyapunovspectrum(bd::Billiard, args...) =
    lyapunovspectrum!(randominside(bd), bd, args...)


################################################################################
## Raw perturbation growth
################################################################################
function perturbationgrowth!(p::AbstractParticle{T}, bd::Billiard{T},
                             t) where {T<:AbstractFloat}

    offset = [SVector{4, T}(1,1,1,1)]

    count = zero(t)
    timecount = zero(T)

    Δ = Vector{SVector{4, T}}() # perturbation vectors
    tim = T[] # sample times
    obst = Int[] # obstacle indices

     # check for pinning before evolution
    if ispinned(p, bd)
        return tim, Δ, obst
    end
    ismagnetic = typeof(p) <: MagneticParticle

    while count < t
        off = offset[1]
        # bounce
        i::Int, tmin::T, cp::SV{T} = next_collision(p, bd)
        relocate!(p, bd[i], tmin, cp, offset)
        timecount += tmin

        # push perturbations before collision
        push!(tim, timecount)
        push!(obst, i)

        # push after propagation evolution
        push!(Δ, offset[1] ./ off)
        off = offset[1]

        resolvecollision!(p, bd[i], offset)
        # push perturbations after collision
        push!(tim, timecount)
        push!(obst, i)
        push!(Δ, offset[1] ./ off)

        ismagnetic && (p.center = find_cyclotron(p))

        # normalize perturbation vector
        offset[1] = offset[1] ./ norm(offset[1])

        count += increment_counter(t, tmin)

    end#time loop

    return tim, Δ, obst
end

"""
    perturbationgrowth([p,] bd, t) -> ts, Rs, is
Calculate the evolution of the perturbation vector `Δ` along the trajectory of `p` in
`bd` for total time `t`. `Δ` is initialised as `[1,1,1,1]`.

If a particle is not given, a random one is picked through [`randominside`](@ref).
Returns empty lists for pinned particles.

## Description
This function *safely* computes the time evolution of a perturbation vector using the
linearized dynamics of the system, as outlined by [1]. Because the dynamics are linear,
we can safely re-normalize the perturbation vector after every collision (otherwise the
perturbations grow to infinity).

Immediately before *and after* every collison, this function computes
* the current time.
* the element-wise ratio of Δ with its previous value
* the obstacle index of the current obstacle
and returns these in three vectors `ts, Rs, is`.

To obtain the *actual* evolution of the perturbation vector you can use
the function `perturbationevolution(Rs)` which simply does
```julia
Δ = Vector{SVector{4,Float64}}(undef, length(R))
Δ[1] = R[1]
for i in 2:length(R)
    Δ[i] = R[i] .* Δ[i-1]
end
```

[1] : Ch. Dellago *et al*, [Phys. Rev. E **53** (1996)](http://link.aps.org/doi/10.1103/PhysRevE.53.1485)
"""
perturbationgrowth(p::AbstractParticle, args...) = perturbationgrowth!(copy(p), args...)
perturbationgrowth(bd::Billiard, args...) =
    perturbationgrowth!(randominside(bd), bd, args...)

function perturbationevolution(R::Vector{SVector{4, T}}) where T
    Δ = Vector{SVector{4,T}}(undef, length(R))
    Δ[1] = R[1]
    @inbounds for i in 2:length(R)
        Δ[i] = R[i] .* Δ[i-1]
    end
    return Δ
end

@deprecate pertubationevolution perturbationevolution

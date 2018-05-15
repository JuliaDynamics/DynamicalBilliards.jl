export lyapunovspectrum!, lyapunovspectrum

const δqind = SV{Int}(1,2)
const δpind = SV{Int}(3,4)

#="""
    specular!(p::AbstractParticle, o::Obstacle, offset::MArray)
Perform specular reflection based on the normal vector of the Obstacle.
The function updates the position and velocity of the particle
together with the components of 4 offset vectors stored in the matrix
`offset` as columns.
"""=#
function specular!(p::AbstractParticle{T}, o::Circular{T},
                   offset::Vector{SVector{4, T}}) where {T<:AbstractFloat}
    n = normalvec(o, p.pos)
    ti = SV{T}(-p.vel[2],p.vel[1])

    cosa = dot(n, -p.vel)
    p.vel = p.vel + 2*cosa*n

    tf = SV{T}(-p.vel[2], p.vel[1])

    for k in 1:4
        δΓ = offset[k]
        # Formulas from Dellago, Posch and Hoover, PRE 53, 2, 1996: 1485-1501 (eq. 27)
        # with norm(p) = 1
        δq  = δΓ[δqind] - 2.*dot(δΓ[δqind],n)*n
        δp  = δΓ[δpind] - 2.*dot(δΓ[δpind],n)*n - curvature(o)*2./o.r*dot(δΓ[δqind],ti)/cosa*tf
        ###
        offset[k] = vcat(δq, δp)
    end
end

@inline curvature(::Semicircle) = -1
@inline curvature(::Disk) = +1


function specular!(p::AbstractParticle{T}, o::Union{InfiniteWall{T},FiniteWall{T}},
                   offset::Vector{SVector{4, T}}) where {T<:AbstractFloat}

    n = normalvec(o, p.pos)
    specular!(p, o)
    for k in 1:4
        δΓ = offset[k]
        # Formulas from Dellago, Posch and Hoover, PRE 53, 2, 1996: 1485-1501 (eq. 20)
        δq  = δΓ[δqind] -  2.*dot(δΓ[δqind],n)*n
        δp  = δΓ[δpind] - 2.*dot(δΓ[δpind],n)*n
        ###
        offset[k] = vcat(δq, δp)
    end
end

#="""
    resolvecollision!(p::AbstractParticle, o::Union{Disk, InfiniteWall}, offset::MArray)
Resolve the collision between particle `p` and obstacle `o` of type *Circular*,
updating the components of the offset vectors stored in the matrix `offset` as columns.
"""=#
function resolvecollision!(p::AbstractParticle{T}, o::Obstacle{T},
                offset::Vector{SVector{4, T}})::Void where {T}

    specular!(p, o, offset)
    return nothing
end

resolvecollision!(p::AbstractParticle{T}, o::PeriodicWall{T},
                  offset::Vector{SVector{4, T}}) where {T} = resolvecollision!(p, o)

"""
    propagate_offset!(offset::MArray{Tuple{4,4},T}, p::AbstractParticle)
Computes the linearized evolution of the offset vectors during propagation for a
time interval `t`
"""
#linear
function propagate_offset!(offset::Vector{SVector{4, T}}, t::T,
                           p::Particle{T}) where T
    for k in 1:4
        δΓ = offset[k]
        offset[k] = SVector{4,T}(δΓ[1] + t*δΓ[3], δΓ[2] + t*δΓ[4], δΓ[3], δΓ[4])
    end
end

#magnetic
#TODO: Use `sincos` for julia v0.7
function propagate_offset!(offset::Vector{SVector{4, T}}, t::T,
                           p::MagneticParticle{T}) where T
    ω = p.omega
    sω = sin(ω*t)
    cω = cos(ω*t)

    for k in 1:4
        δΓ = offset[k]
        offset[k] = SVector{4,T}(δΓ[1] + sω/ω*δΓ[3] + (cω - 1)/ω*δΓ[4],
                                 δΓ[2] + (1 - cω)/ω*δΓ[3] + sω/ω*δΓ[4],
                                 cω*δΓ[3] - sω*δΓ[4],
                                 sω*δΓ[3] + cω*δΓ[4])
    end
end

#="""
    propagate!(p::AbstractParticle{T}, t::T, offset::MArray{Tuple{4,4},T})
Propagate the particle `p` for given time `t`, changing appropriately the the
`p.pos` and `p.vel` fields together with the components of the offset vectors
stored in the `offset` matrix.
"""=#
function propagate!(p::AbstractParticle{T}, t::T,
                    offset::Vector{SVector{4, T}}) where {T<: AbstractFloat}

    propagate!(p, t)
    propagate_offset!(offset, t, p)
end


#="""
    relocate(p::AbstractParticle, o::Obstacle, t, offset::MArray) -> newt
Propagate the particle's position for time `t` (corrected) and update the components
of the `offset` matrix.
"""=#
function relocate!(p::AbstractParticle{T}, o::Obstacle{T}, tmin,
                   offset::Vector{SVector{4, T}}) where {T <: AbstractFloat}

    tmin = relocate!(p, o, tmin)
    propagate_offset!(offset, tmin, p)
    return tmin
end



"""
    lyapunovspectrum!(p::AbstractParticle, bt::Billiard, t)
Returns the finite time lyapunov exponents (averaged over time `t`)
for a given particle in a billiard table.
"""
function lyapunovspectrum!(p::AbstractParticle{T}, bt::Billiard{T}, t::T) where {T<:AbstractFloat}

    offset = [SVector{4, T}(1,0,0,0), SVector{4, T}(0,1,0,0),
              SVector{4, T}(0,0,1,0), SVector{4, T}(0,0,0,1)]

    ismagnetic = typeof(p) <: MagneticParticle
    if t <= 0.0
        error("`evolve!()` cannot evolve backwards in time.")
    end

    count = zero(T)
    λ = zeros(T, 4)

    while count < t
        tmin::T, i::Int = next_collision(p, bt)
        # set counter
        if count +  increment_counter(t, tmin) > t
            break
        end
        ###

        tmin = relocate!(p, bt[i], tmin, offset)
        resolvecollision!(p, bt[i], offset)
        ismagnetic && (p.center = find_cyclotron(p))
        count += increment_counter(t, tmin)

        Q, R = qr(hcat(offset[1], offset[2], offset[3], offset[4]))
        offset[1], offset[2], offset[3], offset[4] = Q[:, 1], Q[:, 2], Q[:, 3], Q[:, 4]
        for i ∈ 1:4
            λ[i] += log(abs(R[i,i]))
        end

    end#time loop

    return λ./count
end

"""
    lyapunovspectrum(p::AbstractParticle, bt::Billiard, t)
Non-mutating version of [`lyapunovspectrum!`](@ref)
"""
lyapunovspectrum(p::AbstractParticle, args...) = lyapunovspectrum!(deepcopy(p), args...)

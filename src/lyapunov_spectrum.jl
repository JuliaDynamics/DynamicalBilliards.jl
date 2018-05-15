export lyapunovspectrum!, lyapunovspectrum

const δqind = SV{Int}(1,2)
const δpind = SV{Int}(3,4)


##Auxiliar Functions ##
"""
```julia
gramschmidt(u::MArray{Tuple{4,4}, T})
```
Apply the Gram-Schmidt procedure to a 4x4 square matrix
"""
function gramschmidt(u::MArray{Tuple{4,4}, T}) where {T<:AbstractFloat}
    w = eye(MMatrix{4,4, T})
    w[:,1] = u[:,1];
    v1 = w[:,1]/norm(w[:,1])
    w[:,2] = u[:,2] - dot(u[:,2],v1)*v1;
    v2 = w[:,2]/norm(w[:,2]);
    w[:,3] = (u[:,3] - dot(u[:,3],v2)*v2 - dot(u[:,3],v1)*v1)
    v3 = w[:,3]/norm(w[:,3])
    w[:,4] = u[:,4] - dot(u[:,4],v3)*v3 - dot(u[:,4],v2)*v2 - dot(u[:,4],v1)*v1
    return w
end

#="""
    specular!(p::AbstractParticle, o::Obstacle, offset::MArray)
Perform specular reflection based on the normal vector of the Obstacle.
The function updates the position and velocity of the particle
together with the components of 4 offset vectors stored in the matrix
`offset` as columns.
"""=#
function specular!(p::AbstractParticle{T}, o::Circular{T},
                   offset::MArray{Tuple{4,4}, T}) where {T<:AbstractFloat}
    n = normalvec(o, p.pos)
    ti = [-p.vel[2],p.vel[1]]
    cosa = dot(n, -p.vel)
    p.vel = p.vel - 2*dot(n, p.vel)*n
    tf = [-p.vel[2], p.vel[1]]
    for k in 1:4
        x = offset[k]
        # Formulas from Dellago, Posch and Hoover, PRE 53, 2, 1996: 1485-1501 (eq. 27)
        # with norm(p) = 1
        a  = x[δpind] - 2.*dot(x[3:4],n)*n + curvature(o)*2./o.r*dot(x[1:2],ti)/cosa*tf
        b  = x[1:2] - 2.*dot(x[1:2],n)*n
        ###
        offset[k] = vcat(a, b)
    end
end

@inline curvature(::Semicircle) = +1
@inline curvature(::Disk) = -1


function specular!(p::AbstractParticle{T}, o::Union{InfiniteWall{T},FiniteWall{T}},
                   offset::MArray{Tuple{4,4}, T}) where {T<:AbstractFloat}

    n = normalvec(o, p.pos)
    specular!(p, o)
    for k in 1:4
        x = [offset[:,k]...]
        # Formulas from Dellago, Posch and Hoover, PRE 53, 2, 1996: 1485-1501 (eq. 20)
        x[1:2]  = x[1:2] -  2.*dot(x[1:2],n)*n
        x[3:4]  = x[3:4] - 2.*dot(x[3:4],n)*n
        ###
        offset[:,k] = x
    end
end

#="""
    resolvecollision!(p::AbstractParticle, o::Union{Disk, InfiniteWall}, offset::MArray)
Resolve the collision between particle `p` and obstacle `o` of type *Circular*,
updating the components of the offset vectors stored in the matrix `offset` as columns.
"""=#
function resolvecollision!(p::AbstractParticle{T},
                o::Union{Disk{T}, InfiniteWall{T}, FiniteWall{T}, Semicircle{T}},
                offset::MArray{Tuple{4,4}, T})::Void where {T<:AbstractFloat}

    specular!(p, o, offset)
    return nothing
end

resolvecollision!(p::AbstractParticle{T}, o::PeriodicWall{T},
                  offset::MArray{Tuple{4,4}, T}) where {T} = resolvecollision!(p, o)

"""
    propagate_offset!(offset::MArray{Tuple{4,4},T}, p::AbstractParticle)
Computes the linearized evolution of the offset vectors during propagation for a
time interval `t`
"""
#linear
function propagate_offset!(offset::MArray{Tuple{4,4},T}, t::T,
                           p::Particle{T}) where T
    for k in 1:4
        Γ = offset[k]
        offset[k] = SV{T}(Γ[1] + t*Γ[3], Γ[2] + t*Γ[4], Γ[3], Γ[4])
    end
end

#magnetic
#TODO: Use `sincos` for julia v0.7
function propagate_offset!(offset::MArray{Tuple{4,4},T}, t::T,
                           p::MagneticParticle{T}) where T
    ω = p.omega
    sω = sin(ω*t)
    cω = cos(ω*t)
    J = [1.0 0.0     sω/ω  (cω-1)/ω;
         0.0 1.0  (1-cω)/ω     sω/ω;
         0.0 0.0     cω         -sω;
         0.0 0.0     sω          cω]

         # DELETE J
         # INSTEAD FIND ANALYTICALLY WHAT DOES J*x does for a general x

    for k in 1:4
        Γ = offset[k]
        offset[k] = SV{T}(whatever you get)
    end
end

#="""
    propagate!(p::AbstractParticle{T}, t::T, offset::MArray{Tuple{4,4},T})
Propagate the particle `p` for given time `t`, changing appropriately the the
`p.pos` and `p.vel` fields together with the components of the offset vectors
stored in the `offset` matrix.
"""=#
function propagate!(p::AbstractParticle{T}, t::T,
                    offset::MArray{Tuple{4,4},T}) where {T<: AbstractFloat}

    propagate!(p, t)
    propagate_offset!(offset, t, p)
end


#="""
    relocate(p::AbstractParticle, o::Obstacle, t, offset::MArray) -> newt
Propagate the particle's position for time `t` (corrected) and update the components
of the `offset` matrix.
"""=#
function relocate!(p::AbstractParticle{T}, o::Obstacle{T}, tmin,
                   offset::MArray{Tuple{4,4}, T}) where {T <: AbstractFloat}

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

    offset = eye(MMatrix{4,4, T}) #The unit vectors in the 4 directions
    ismagnetic = typeof(p) <: MagneticParticle
    if t <= 0.0
        error("`evolve!()` cannot evolve backwards in time.")
    end

    count = zero(T)
    norms = ones(T, 1,4)#Where the norms of the offset vectors will be stored

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

        if typeof(bt[i]) <: Wall
            continue
        else
            offset = gramschmidt(offset)

            Q, R = qr(hcat(offset[1], v[2], v[3], v[4]))
            offset[1], v[2], v[3], v[4] = Q[:, 1], Q[:, 2], Q[:, 3], Q[:, 4]
            for i in 1:k
                λ[i] += log(abs(R[i,i]))
            end

            instantaneous_norms = [norm(offset[:,j]) for j in 1:4]
            norms = vcat(norms,instantaneous_norms')
            for j in 1:4
                offset[:,j] = offset[:,j]/norm(offset[:,j])
            end
        end
    end#time loop

    tmin = t - count
    propagate!(p, tmin, offset)
    instantaneous_norms = [norm(offset[:,j]) for j in 1:4]
    norms = vcat(norms,instantaneous_norms')
    for j in 1:4
        offset[:,j] = offset[:,j]/norm(offset[:,j])
    end

    exps = zeros(4)

    for k in 1:4
        exps[k] = sum(log.(norms[:,k]))/t
    end

    return exps
end

"""
    lyapunovspectrum(p::AbstractParticle, bt::Billiard, t)
Non-mutating version of [`lyapunovspectrum!`](@ref)
"""
lyapunovspectrum(p::AbstractParticle, args...) = lyapunovspectrum!(deepcopy(p), args...)

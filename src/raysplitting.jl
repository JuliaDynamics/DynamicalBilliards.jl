export isphysical, reset_billiard!
export RaySplitter, raysplit_indices
export law_of_refraction

#####################################################################################
# RaySplitter structures
#####################################################################################
"""
    RaySplitter(idxs, transmission, refraction [, newangular]; affect)
Return a `RaySplitter` instance, used to perform raysplitting.
`idxs` is a `Vector{Int}` with the indices of the obstacles
that this `RaySplitter` corresponds to.

`transmission`, `refraction` and `newangular` are **functions**. Let
`φ` be the angle of incidence and `ω` be the angular velocity
and `pflag` the propagation flag (before transmission).
The functions have the following signatures:

1. `transmission(φ, pflag, ω) -> T`, transmission probability.
2. `refraction(φ, pflag, ω) -> θ`, refraction angle.
   This angle is *relative* to the normal vector.
3. `newangular(ω, pflag) -> newω`, new angular velocity after transmission.

The above three functions use the **same convention**: the argument `pflag` is the
one the obstacle has **before transmission**. For example, if a particle is
outside an [`Antidot`](@ref) (with `pflag = true` here) and is transmitted inside
the `Antidot` (`pflag` becomes `false` here), then all three functions will be
given their second argument (the Boolean one) as `true`!

`affect` is a function, and denotes which obstacles
of the billiard are affected when transmission occurs at obstacle
`i` (for which obstacles should the field `pflag` be reversed).
Defaults to `idxs = (i) -> i`, i.e. only the colliding obstacle is affected.
If you want many obstacles to be affected you could write
`idxs = (i) -> SVector(2,3,5)`, etc.
Keep in mind that the only values of `i` that can be passed into this function
are the ones that are given in the argument `idxs`!
"""
struct RaySplitter{T, Φ, Ω, A}
    oidx::Vector{Int}
    transmission::T
    refraction::Φ
    newω::Ω
    affect::A
end

@inline default_affect(i) = i
@inline default_angular(ω, pflag) = ω

function RaySplitter(idxs, tr, ref, newangular = default_angular;
                                    affect = default_affect)
    for i ∈ idxs
        i ∈ affect(i) || throw(ArgumentError(
        "All indices that correspond to this RaySplitter must also be affected!"))
    end
    typeof(idxs) == Int && (idxs = [idxs])
    return RaySplitter(sort(idxs), tr, ref, newangular, affect)
end

#pretty print:
function Base.show(io::IO, ::MIME"text/plain", ray::RaySplitter)
    ps = 15
    angtext = ray.newω == default_angular ? "default" : "$(ray.newω)"
    afftext = ray.affect == default_affect ? "default" : "$(ray.affect)"
    print(io, "RaySplitter for indices $(ray.oidx)"*"\n",
    rpad(" transmission: ", ps)*"$(ray.transmission)\n",
    rpad(" refraction: ", ps)*"$(ray.refraction)\n",
    rpad(" new angular: ", ps)*angtext*"\n",
    rpad(" affect: ", ps)*afftext*"\n"
    )
end

function Base.show(io::IO, ray::RaySplitter)
    print(io, "RaySplitter for indices $(ray.oidx)")
end


"""
    raysplit_indices(bd::Billiard, raysplitters::Tuple)
Create a vector of integers. The `i`th entry tells you which entry of the
`raysplitters` tuple is associated with the `i`th obstacle of the billiard.

If the `i`th entry is `0`, this means that the obstacle does not do raysplitting.
"""
function raysplit_indices(bd::Billiard, raysplitters::Tuple)
    O = zeros(Int, length(bd.obstacles))
    for (k, rayspl) ∈ enumerate(raysplitters)
        O[rayspl.oidx] .= k
    end
    return O
end

allaffected(ray::RaySplitter) = union(ray.affect(i) for i in ray.oidx)

"""
    acceptable_raysplitter(raysplitters, bd::Billiard)
Return `true` if the given `raysplitters`
can be used in conjuction with given billiard `bd`.
"""
function acceptable_raysplitter(raysplitters::Tuple, bd::Billiard)

    for ray in raysplitters
        for i in (ray.oidx ∪ ray.affect(ray.oidx))
            if !supports_raysplitting(bd[i])
                error("Obstacle at index $i of given billiard table"*
                " does not have a field `pflag`"*
                " and therefore does not support ray-splitting."*
                " However a `RaySplitter` uses it!")
            end
        end
    end

    idx = raysplit_indices(bd, raysplitters)

    # Make sure that no indices are shared in the corresponding obstacles
    vecs = [ray.oidx for ray in raysplitters]
    if length(unique(vcat(vecs...))) != sum(length.(vecs))
        error("Different `RaySplitter`s cannot correspond to the same "*
        "obstacles!")
    end

    true
end

#####################################################################################
# Resolve collisions
#####################################################################################
timeprec_rayspl(::Particle{T}) where {T} = timeprec(T)
timeprec_rayspl(::MagneticParticle{T}) where {T} = timeprec_forward(T)

const CLAMPING_ANGLE = 0.1
angleclamp(φ::T) where {T} = clamp(φ, -π/2 + T(CLAMPING_ANGLE), π/2 - T(CLAMPING_ANGLE))

function incidence_angle(p::AbstractParticle{T}, a::Obstacle{T})::T where {T}
    # Raysplit Algorithm step 1: Determine incidence angle (0 < φ < π/4)
    n = normalvec(a, p.pos)
    inverse_dot = clamp(dot(p.vel, -n), -1.0, 1.0)
    φ = acos(inverse_dot)
    # Raysplit Algorithm step 2: get correct sign
    if cross2D(p.vel, n) < 0
        φ *= -1
    end
    return φ
end

function istransmitted(p::AbstractParticle{T}, a::Obstacle{T}, Tr::F) where {T, F}
    ω = typeof(p) <: MagneticParticle ? p.omega : zero(T)
    φ = incidence_angle(p, a)
    # Raysplit Algorithm step 3: check transmission probability
    trans = Tr(φ, a.pflag, ω) > rand()
end

# Raysplit Algorithm step 4: relocate the particle _inside_ the obstacle
# if ray-splitting happens (see `ineq` variable)
function relocate_rayspl!(p::AbstractParticle{T}, o::Obstacle{T},
    trans::Bool = false) where {T}

    ineq = 2trans - 1
    d = distance(p.pos, o)
    notokay = ineq*d > 0.0
    if notokay
        n = normalvec(o, p.pos)
        p.pos -= d * n
    end
    return notokay
end


function resolvecollision!(p::AbstractParticle{T}, bd::Billiard{T}, colidx::Int,
    trans::Bool, rayspl::RaySplitter) where {T<:AbstractFloat}

    a = bd[colidx]
    ismagnetic = typeof(p) <: MagneticParticle
    ω = ismagnetic ? p.omega : T(0)

    # Raysplit Algorithm step 5: recompute angle of incidence
    φ = incidence_angle(p, a)

    if trans #perform raysplitting
        # Raysplit Algorithm step 6: find transmission angle in relative angles
        theta = angleclamp(rayspl.refraction(φ, a.pflag, ω))

        # Raysplit Algorithm step 7: reverse the Obstacle propagation flag
        # for all obstacles dictated by the RaySplitter
        # (this also reverses `normalvec` !)
        for oi ∈ rayspl.affect(colidx)
            bd[oi].pflag = ! bd[oi].pflag
        end

        # Raysplit Algorithm step 8: find transmission angle in real-space angles
        n = normalvec(a, p.pos) #notice that this is reversed! It's the new normalvec!
        Θ = theta + atan(n[2], n[1])

        # Raysplit Algorithm step 9: Perform refraction
        p.vel = SVector{2,T}(cos(Θ), sin(Θ))

        # Raysplit Algorithm step 10: Set new angular velocity
        if ismagnetic
            ω = rayspl.newω(p.omega, !a.pflag)  # notice the exclamation mark
            p.omega = ω
            p.r = abs(1/ω)
        end
    else # No ray-splitting:
        #perform specular
        specular!(p, a)
        end
    return
end

# Ray-splitting version of bounce!
# raysplitters is a tuple of RaySplitter. raysidx is a vector that given the obstacle
# index it tells you which raysplitter to choose from the tuple OR to not
# do raysplitting at all (for returned index 0)
function bounce!(p::AbstractParticle{T}, bd::Billiard{T},
    raysidx::Vector{Int}, raysplitters::Tuple) where {T}

    i::Int, tmin::T, cp::SV{T} = next_collision(p, bd)

    if tmin == Inf
        return i, tmin, p.pos, p.vel

    elseif raysidx[i] != 0

        propagate!(p, cp, tmin)

        trans::Bool = istransmitted(p, bd[i], raysplitters[raysidx[i]].transmission)

        relocate_rayspl!(p, bd[i], trans)
        resolvecollision!(p, bd, i, trans, raysplitters[raysidx[i]])
    else
        relocate!(p, bd[i], tmin, cp)
        resolvecollision!(p, bd[i])
    end
    typeof(p) <: MagneticParticle && (p.center = find_cyclotron(p))
    return i, tmin, p.pos, p.vel
end

#####################################################################################
# Evolve raysplitting
#####################################################################################
@inline raysplit_indices(bd, ::Nothing) = nothing
@inline bounce!(p, bd, ::Nothing, ::Nothing) = bounce!(p, bd)
@inline acceptable_raysplitter(::Nothing, bd::Billiard) = true

evolve!(p, bd, t, ray::RaySplitter; warning = false) =
    evolve!(p, bd, t, (ray, ); warning = warning)
function evolve!(p::Matrix{T}, bd::Billiard{T}, t, raysplitters::Tuple;
    warning = false) where {T}

    if t <= 0
        throw(ArgumentError("`evolve!()` cannot evolve backwards in time."))
    end
    # Check if raysplitters are acceptable
    acceptable_raysplitter(raysplitters, bd)

    ismagnetic = typeof(p) <: MagneticParticle

    raysidx = raysplit_indices(bd, raysplitters)

    rt = T[]; push!(rt, 0)
    rpos = SVector{2,T}[]; push!(rpos, p.pos)
    rvel = SVector{2,T}[]; push!(rvel, p.vel)
    ismagnetic && (omegas = T[]; push!(omegas, p.omega))

    count = zero(t)
    t_to_write = zero(T)

    while count < t

        i, tmin, pos, vel = bounce!(p, bd, raysidx, raysplitters)
        t_to_write += tmin

        if ismagnetic && tmin == Inf
            warning && warn("Pinned particle in evolve! (Inf. col. t)")
            push!(rpos, pos); push!(rvel, vel)
            push!(rt, tmin); push!(omegas, p.omega)
            return (rt, rpos, rvel, omegas)
        end

        if isperiodic(bd) && i ∈ bd.peridx
            # Pinned particle:
            if ismagnetic && t_to_write ≥ 2π/absω
                warning && warn("Pinned particle in evolve! (completed circle)")
                push!(rpos, pos); push!(rvel, vel)
                push!(rt, tmin); push!(omegas, p.omega)
                return (rt, rpos, rvel, omegas)
            end
            #If not pinned, continue (do not write for PeriodicWall)
            continue
        else
            push!(rpos, p.pos + p.current_cell)
            push!(rvel, p.vel); push!(rt, t_to_write);
            ismagnetic && push!(omegas, p.omega)
            count += increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end

    end#time loop
    if ismagnetic
        return (rt, rpos, rvel, omegas)
    else
        return (rt, rpos, rvel)
    end
end


########################
# is physical, etc.
########################

function supports_raysplitting(obst::Obstacle)
  n = fieldnames(typeof(obst))
  in(:pflag, n)
end

"""
    reset_billiard!(bd)
Sets the `pflag` field of all ray-splitting obstacles of a billiard table
to `true`.
"""
function reset_billiard!(bd::Billiard)
    for obst in bd
        supports_raysplitting(obst) && (obst.pflag = true)
    end
end

"""
    isphysical(raysplitter(s))
Return `true` if the given `raysplitters` have physically
plausible properties.

Specifically, check if (φ is the incidence angle, θ the refraction angle):

* Critical angle means total reflection: If θ(φ) ≥ π/2 then Tr(φ) = 0
* Transmission probability is even function: Tr(φ) ≈ Tr(-φ) at ω = 0
* Refraction angle is odd function: θ(φ) ≈ -θ(-φ) at ω = 0
* Ray reversal is true: θ(θ(φ, pflag, ω), !pflag, ω) ≈ φ
* Magnetic conservation is true: (ω_new(ω_new(ω, pflag), !pflag) ≈ ω
"""
isphysical(ray::RaySplitter) = isphysical((ray,))
function isphysical(rays::Tuple)
  for (i, ray) in enumerate(rays)
    scatter = ray.refraction
    tr = ray.transmission
    om = ray.newω
    range = -1.5:0.01:1.5
    orange = -1.0:0.1:1.0
    for pflag in (true, false)
      for ω in orange
        for φ in range
          # Calculate refraction angle:
          θ = scatter(φ, pflag, ω)
          # Calculate transmission probability:
          T = tr(φ, pflag, ω)
          # Check critical angle:
          if θ >= π/2 && T > 0
            es = "Refraction angle ≥ π/2 and T > 0 !\n"
            es*= "For index = $i, tested with φ = $φ, pflag = $pflag, ω = $ω"
            println(es)
            return false
          end
            # Check symmetry:
            if ω==0
              if !isapprox(θ, -scatter(-φ, pflag, ω))
                es = "Scattering angle function is not odd!\n"
                es *="For index = $i, tested with φ = $φ, pflag = $pflag, ω = $ω"
                println(es)
                return false
              end
              if !isapprox(T, tr(-φ, pflag, ω))
                es = "Transmission probability function is not even!\n"
                es *="For index = $i, tested with φ = $φ, pflag = $pflag, ω = $ω"
                println(es)
                return false
              end
            end
            # Check ray-reversal:
            if !isapprox(scatter(θ, !pflag, ω), φ)
              es = "Ray-reversal does not hold!\n"
              es *="For index = $i, tested with φ = $φ, pflag = $pflag, ω = $ω"
              println(es)
              return false
            end
            if !isapprox(om(om(ω, pflag), !pflag), ω)
              es = "Magnetic reversal does not hold!\n"
              es *="For index = $i, tested with φ = $φ, pflag = $pflag, ω = $ω"
              println(es)
              return false
            end
        end#φ range
      end#ω range
    end#pflag range
  end#obstacle range
  return true
end


################################################################################

"""
    law_of_refraction(n1, n2 = 1.0) -> t, r

Create transmission and refraction functions `t, r` that follow Snell's
law, i.e. the transmission probability is set to 1.0 except
for the case of total internal reflection. 

`n1` is the index of refraction for the `pflag = false` side of an obstacle,
while `n2` is the index of refraction for `pflag = true`.
"""
function law_of_refraction(n1, n2 = 1.0)
    # n1 is "inside" the obstacle, n2 is "outside"

    if n1 < 1.0 || n2 < 1.0
        error("You have just given a physicist a headache.")
    end

    # Snell's law 
    refraction(φ, pflag, ω) =
        asin(clamp((pflag ? (n2/n1) : (n1/n2) )* sin(φ), -1.0, 1.0))
    
    # total internal reflection
    function transmission(φ, pflag, ω)
        nratio = pflag ? (n1/n2) : (n2/n1)
        
        if abs(nratio) < 1 &&  abs(φ) >= asin(nratio)
            return 0.0
        else
            return 1.0
        end
    end
    
    return transmission, refraction
end

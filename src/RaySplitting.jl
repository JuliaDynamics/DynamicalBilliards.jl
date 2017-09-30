export isphysical, acceptable_raysplitter, reset_billiard!

########################
# Resolve collisions
########################

±(x) = 2x - 1

function relocate_rayspl!(
    p::Particle{T}, o::Obstacle{T}, trans::Bool = false)::T where {T}

    ineq = ±(trans)
    newpos = p.pos; newt = zero(T)
    i = 1
    while ineq*distance(newpos, o) > 0
        newt += i*ineq*timeprec(T)
        newpos = propagate_pos(p.pos, p, newt)
        i *= 10
    end
    propagate!(p, newpos, newt)
    return newt
end

function relocate_rayspl!(
    p::MagneticParticle{T}, o::Obstacle{T}, trans::Bool = false)::T where {T}

    ineq = ±(trans)
    newpos = p.pos; newt = zero(T)
    i = 1
    while ineq*distance(newpos, o) > 0
        newt += ineq*timeprec_forward(T)
        newpos = propagate_pos(p.pos, p, newt)
        i *= 10
    end
    propagate!(p, newpos, newt)
    return newt
end

function incidence_angle(p::AbstractParticle{T}, a::Obstacle{T})::T where {T}
    # Raysplit Algorithm step 1: Determine incidence angle (0 < θ < π/4)
    n = normalvec(a, p.pos)
    inverse_dot = clamp(dot(p.vel, -n), -1.0, 1.0)
    φ = acos(inverse_dot)
    # Raysplit Algorithm step 2:
    if cross2D(p.vel, n) < 0
        φ *= -1
    end
    return φ
end

function istransmitted(p::Particle{T}, a::Obstacle{T}, Tr::Function) where {T}
    φ = incidence_angle(p, a)
    # Raysplit Algorithm step 3: check transmission probability
    trans = Tr(φ, a.pflag) > rand()
    return trans, φ
    # comment: more accurate would be to calculate incidence angle again
    # within relocate!(). But the angle changes so little that this must have
    # almost zero impact.
end
function istransmitted(p::MagneticParticle{T}, a::Obstacle{T}, Tr::Function) where {T}
    φ = incidence_angle(p, a)
    # Raysplit Algorithm step 3: check transmission probability
    trans = Tr(φ, a.pflag, p.omega) > rand()
    return trans, φ
    # comment: more accurate would be to calculate incidence angle again
    # within relocate!(). But the angle changes so little that this must have
    # almost zero impact.
end



# Resolve collision for ray-splitting with Normal
function resolvecollision!(p::AbstractParticle{T},
    a::Obstacle{T}, φ::T, trans::Bool,
    θ::Function, newω::Function = (ω, bool) -> ω) where {T<:AbstractFloat}

    ω = typeof(p) <: Particle ? zero(T) : p.omega

    if trans #perform raysplitting
        # Raysplit Algorithm step 4: find transmission angle in relative angles
        theta = θ(φ, a.pflag, ω)
        # Raysplit Algorithm step 5: reverse the Obstacle propagation flag
        a.pflag = !a.pflag
        # Raysplit Algorithm step 6: find transmission angle in real-space angles
        n = normalvec(a, p.pos) #notice that this is reversed! It's the new!
        Θ = theta + atan2(n[2], n[1])
        # Step 7+8: Relocation has been done already by the
        # call of `tmin = relocate!(p, bt[colobst_idx], tmin, trans)`
        # Sidenote: The above call must happen BEFORE flag reversal!

        # Raysplit Algorithm step 9: Perform refraction
        p.vel = SVector{2,T}(cos(Θ), sin(Θ))
        # Raysplit Algorithm step 10: Set new angular velocity
        if typeof(p) <: MagneticParticle
            p.omega = newω(p.omega, !a.pflag)  # notice the exclamation mark
        end
    else # No ray-splitting:
        #perform specular
        specular!(p, a)
        end
    return
end

###########
# Straight
###########

# evolve For Particle and Ray-Splitting:
function evolve!(p::Particle{T}, bt, t,
    ray::Dict) where {T}

    const debug = false

    if t <= 0
    throw(ArgumentError("`evolve!()` cannot evolve backwards in time."))
    end

    rt = T[]
    rpos = SVector{2,T}[]
    rvel = SVector{2,T}[]
    push!(rpos, p.pos)
    push!(rvel, p.vel)
    push!(rt, 0)

    count = zero(t)
    t_to_write = zero(T)

    while count < t
        # Declare these because `bt` is of un-stable type!
        tmin::T, i::Int = next_collision(p, bt)

        debug && println("Min. col. t with obst $(bt[i].name) = $tmin")

        if haskey(ray, i)
            propagate!(p, tmin)
            trans, φ = istransmitted(p, bt[i], ray[i][1])
            debug && println("Angle of incidence: $(φ), transmitted? $trans")
            if debug
                println("Currently, pflag is $(bt[colobst_idx].pflag)")
                println("After collision is resolved, it will be the opposite")
            end
            newt = relocate_rayspl!(p, bt[i], trans)
            resolvecollision!(p, bt[i], φ, trans, ray[i][2])
            t_to_write += tmin + newt
        else
            tmin = relocate!(p, bt[i], tmin)
            resolvecollision!(p, bt[i])
            t_to_write += tmin
        end

        debug && println()
        if typeof(bt[i]) <: PeriodicWall
            continue
        else
            push!(rpos, p.pos + p.current_cell)
            push!(rvel, p.vel)
            push!(rt, t_to_write)
            count += increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end

    end#time loop
    return (rt, rpos, rvel)
end

##########
# Magnetic
##########

# evolve For MagneticParticle and Ray-Splitting
function evolve!(p::MagneticParticle{T}, bt::Vector{<:Obstacle{T}},
t, ray::Dict; warning::Bool = false) where {T}

    if t <= 0
    error("`evolve!()` cannot evolve backwards in time.")
    end
    if isperiodic(bt) && T == BigFloat
        error("raysplitting + magnetic + BigFloat is not supported :(")
    end

    const debug = false

    omegas = T[]
    rt = T[]
    rpos = SVector{2,T}[]
    rvel = SVector{2,T}[]
    push!(rpos, p.pos)
    push!(rvel, p.vel)
    push!(rt, zero(T))
    push!(omegas, p.omega)

    count = zero(t)
    t_to_write = zero(T)

    while count < t
        # Declare these because `bt` is of un-stable type!
        tmin::T, i::Int = next_collision(p, bt)

        debug && println("Min. col. t with obst $(bt[i].name) = $tmin")

        if tmin == Inf
            warning && warn("Pinned particle in evolve! (Inf col t)")
            push!(rpos, rpos[end])
            push!(rvel, rvel[end])
            push!(rt, Inf)
            push!(omegas, p.omega)
            return (rt, rpos, rvel, omegas)
        end

        if haskey(ray, i)
            propagate!(p, tmin)
            trans, φ = istransmitted(p, bt[i], ray[i][1])
            debug && println("Angle of incidence: $(φ), transmitted? $trans")
            if debug
                println("Currently, pflag is $(bt[colobst_idx].pflag)")
                println("After collision is resolved, it will be the opposite")
            end
            newt = relocate_rayspl!(p, bt[i], trans)
            resolvecollision!(p, bt[i], φ, trans, ray[i][2], ray[i][3])
            t_to_write += tmin + newt
        else
            tmin = relocate!(p, bt[i], tmin)
            resolvecollision!(p, bt[i])
            t_to_write += tmin
        end


        # Write output only if the collision was not made with PeriodicWall
        if typeof(bt[i]) <: PeriodicWall
            # Pinned particle:
            if t_to_write >= 2π/abs(p.omega)
            warning && warn("Pinned particle in evolve! (completed circle)")
            push!(rpos, rpos[end])
            push!(rvel, rvel[end])
            push!(rt, Inf)
            push!(omegas, p.omega)
            return (rt, rpos, rvel, omegas)
            end
            #If not pinned, continue (do not write for PeriodicWall)
            continue
        else
            push!(rpos, p.pos + p.current_cell)
            push!(rvel, p.vel)
            push!(rt, t_to_write)
            push!(omegas, p.omega)
            count += increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end

    end#time loop
    return (rt, rpos, rvel, omegas)
end

function construct(t::Vector{T}, poss::Vector{SVector{2,T}},
vels::Vector{SVector{2,T}}, omegas::Vector{T}, dt=0.01) where T

    xt = [poss[1][1]]
    yt = [poss[1][2]]
    vxt= [vels[1][1]]
    vyt= [vels[1][2]]
    ts = [t[1]]
    ct = cumsum(t)

    for i in 2:length(t)
        ω = omegas[i-1]
        φ0 = atan2(vels[i-1][2], vels[i-1][1])
        x0 = poss[i-1][1]; y0 = poss[i-1][2]
        colt=t[i]

        t0 = ct[i-1]
        # Construct proper time-vector
        if colt >= dt
            timevec = collect(0:dt:colt)[2:end]
            timevec[end] == colt || push!(timevec, colt)
        else
            timevec = colt
        end

        for td in timevec
            push!(vxt, cos(ω*td + φ0))
            push!(vyt, sin(ω*td + φ0))
            push!(xt, sin(ω*td + φ0)/ω + x0 - sin(φ0)/ω)  #vy0 is sin(φ0)
            push!(yt, -cos(ω*td + φ0)/ω + y0 + cos(φ0)/ω) #vx0 is cos(φ0)
            push!(ts, t0 + td)
        end#collision time
    end#total time
    return xt, yt, vxt, vyt, ts
end

########################
# is physical, etc.
########################

function supports_raysplitting(obst::Obstacle)
  n = fieldnames(typeof(obst))
  in(:pflag, n)
end

"""
    reset_billiard!(bt)
Sets the `pflag` field of all ray-splitting obstacles of a billiard table
to `true`.
"""
function reset_billiard!(bt::Vector{<:Obstacle})
    for obst in bt
        supports_raysplitting(obst) && (obst.pflag = true)
    end
end

"""
    acceptable_raysplitter(raysplitter, bt)
Return `true` if the given ray-splitting dictionary `raysplitter`
can be used in conjuction with given billiard table `bt`.
"""
function acceptable_raysplitter(ray::Dict{Int, Any}, bt::Vector{Obstacle})
    for i in keys(ray)
        if !supports_raysplitting(bt[i])
            print("Obstacle at index $i of given billiard table")
            println("does not have a field `pflag`")
            println("and therefore does not support ray-splitting.")
            return false
        end
    end
    true
end

"""
    isphysical(raysplitter::Dict; only_mandatory = false)
Return `true` if the given ray-splitting dictionary has physically
plausible properties.

Specifically, check if (φ is the incidence angle, θ the refraction angle):

* Critical angle means total reflection: If θ(φ) ≥ π/2 then Tr(φ) = 0
* Transmission probability is even function: Tr(φ) ≈ Tr(-φ) at ω = 0
* Refraction angle is odd function: θ(φ) ≈ -θ(-φ) at ω = 0
* Ray reversal is true: θ(θ(φ, pflag, ω), !pflag, ω) ≈ φ
* Magnetic conservation is true: (ω_new(ω_new(ω, pflag), !pflag) ≈ ω

The first property is mandatory to hold for any setting and is always checked.
The rest are checked if `only_mandatory = false`.
"""
function isphysical(ray::Dict; only_mandatory = false)
  for i in keys(ray)
    scatter = ray[i][2]
    tr = ray[i][1]
    om = ray[i][3]
    range = -1.5:0.01:1.5
    orange = -1.0:0.1:1.0
    display_er = true
    for pflag in [true, false]
      for ω in orange
        for φ in range
          θ::Float64 = 0.0
          # Calculate refraction angle:
          try
            θ = scatter(φ, pflag, ω)
          catch er
            if display_er
              ws = "Got error message: $er\n"
              ws*= "while calculating the refraction angle with settings:\n"
              ws*= "index = $i, φ = $φ, pflag = $pflag, ω = $ω\n"
              ws*= "Similar warnings will be skipped as long as the Tr. prob. is 0."
              warn(ws)
            end
            display_er = false
            T = tr(φ, pflag, ω)
            if T!= 0
              println("Got error message: $er")
              println("while calculating the refraction angle with settings:")
              println("index = $i, φ = $φ, pflag = $pflag, ω = $ω")
              println("PROBLEM: Transmission prob. was not 0 for these settings!")
              return false
            else
              continue
            end
          end
          # Calculate transmission probability:
          T = tr(φ, pflag, ω)
          # Check critical angle:
          if θ >= π/2 && T > 0
            es = "Refraction angle >= π/2 and T > 0 !\n"
            es*= "For index = $i, tested with φ = $φ, pflag = $pflag, ω = $ω"
            println(es)
            return false
          end
          if !only_mandatory
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
          end
        end#φ range
      end#ω range
    end#pflag range
  end#obstacle range
  return true
end

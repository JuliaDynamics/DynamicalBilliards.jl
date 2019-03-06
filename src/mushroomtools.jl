"""
    MushroomTools
Module containing many functions helpful in simulating (perfect) mushroom
billiards, see [`billiard_mushroom`](@ref).
Contains stuff like initializing efficiently regular or chaotic particles
and functions that return the corresponding chaotic or regular phase-space volumes
or portions.
The functions [`V_3D_tot`](@ref) and [`V_3D_reg`](@results) use equations derived
in ref. [1].

Made by Lukas Hupe.

## References
[1] A. Barnett & T. Betcke, [Chaos **17**, 043125 (2007)](https://doi.org/10.1063/1.2816946).
"""
module MushroomTools

using DynamicalBilliards
using DynamicalBilliards: SV, cossin
using LinearAlgebra

#="""
    is_regular(p::Particle, o::Semicircle, w)
checks whether the particle `p` is on an integrable orbit in a mushroom billiard
with cap `o` and stem width `w`
"""=#
function is_regular(p::Particle{T}, l, w, r) where {T}
    δ = SV{T}(0.0, l) - p.pos
    vperp = [-p.vel[2],p.vel[1]]
    if norm(δ) > r
        #debug && println("not inside circle: pos = $(p.pos)")
        return false
    elseif p.pos[2] < l
        #debug && println("wrong half of circle")
        return false
    elseif abs(dot(vperp, δ)) < w/2
        #debug && println("not outside critical radius")
        return false
    else
        #debug && println("Perfectly happy with $(p.pos) and $(p.vel)")
        return true
    end
end

#Formulae from quantum mushroom paper
"""
    V_3D_tot(l,w,r)
Return the total phasespace volume (3D) of a [`billiard_mushroom`](@ref)
parameterized by `(l,w,r)`.
"""
V_3D_tot(l,w,r) = 2π*(l*w + (1/2)*π*r^2)
V_3D_reg(l,w,r) = 2π*r^2*(acos(w/(2r)) - w/(2r)*sqrt(1 - (w^2)/(4r^2)))
V_3D_cha(l,w,r) = V_3D_tot(l,w,r) - V_3D_reg(l,w,r)

"""
    g_r_3D(l, w, r)
Return the regular phasespace portion of the full (3D) phase-space of a
[`billiard_mushroom`](@ref) with stem length `l`, stem width `w` and cap radious `r`.
This result is known analytically, see [`MushroomTools`](@ref) for references.
"""
g_r_3D(l,w,r) = V_3D_reg(l,w,r)/V_3D_tot(l,w,r)

"""
    g_c_3D(l, w, r)
Return the chaotic phasespace portion of the full (3D) phase-space of a
[`billiard_mushroom`](@ref) with stem length `l`, stem width `w` and cap radious `r`.
This result is known analytically, see [`MushroomTools`](@ref) for references.
"""
g_c_3D(l,w,r) = 1 - g_r_3D(l,w,r)

"""
    V_2D_tot(l,w,r)
Return the total boundary map volume (2D) of a [`billiard_mushroom`](@ref)
parameterized by `(l,w,r)`.
"""
V_2D_tot(l, w, r) = 2(π*r + 2r + 2l)
V_2D_reg(l, w, r) = 2π*r*(1 - w/2r) + 2sqrt(4r*r - w*w) -2w*acos(w/2r)
V_2D_cha(l,w,r) = V_2D_tot(l,w,r) - V_2D_reg(l,w,r)

"""
    g_r_2D(l, w, r)
Return the regular phasespace portion of the boundary map (2D) of a
[`billiard_mushroom`](@ref) with stem length `l`, stem width `w` and cap radious `r`.
This result is known analytically, see [`MushroomTools`](@ref) for references.
"""
g_r_2D(l,w,r) = V_2D_reg(l,w,r)/V_2D_tot(l,w,r)

"""
    g_c_2D(l, w, r)
Return the chaotic phasespace portion of the boundary map (2D) of a
[`billiard_mushroom`](@ref) with stem length `l`, stem width `w` and cap radious `r`.
This result is known analytically, see [`MushroomTools`](@ref) for references.
"""
g_c_2D(l,w,r) = 1 - g_r_2D(l,w,r)

#applying Kac's Lemma to the stem
τ_r_2D_kac(l,w,r) = (V_2D_tot(l,w,r) - V_2D_reg(l,w,r)) / (2w + 4l)


#="""
    insidemushroom(pos::StaticVector, l, w, r)
Returns `true` if pos is within the mushroom parameterised by `l`,  `w` and `r`
"""=#
function insidemushroom(pos::SV{T}, l::T, w::T, r::T) where {T <: AbstractFloat}
    if pos[2] > 0
        if pos[2] <= l
            if abs(pos[1]) <= w/2
                return true
            end
        elseif pos[2] <= l + r
            if pos[1]^2 + (pos[2] - l)^2 <= r^2
                return true
            end
        end
    end
    return false
end

"""
    randin_mushroom(l, w, r [, ω])
Generate a random particle within the [`billiard_mushroom`](@ref)
parameterised by `l`, `w` and `r`.
If `ω` is given the particle is magnetic instead.

This function is much more efficient than [`randominside`](@ref).
"""
function randin_mushroom(l::T = 1.0, w::T = 0.2, r::T = 1.0) where {T <: AbstractFloat}
    pos, vel = _randin_mushroom(l, w, r)
    return Particle(pos, vel, zero(SV{T}))
end
function randin_mushroom(l::T, w::T, r::T, ω::T) where {T <: AbstractFloat}
    pos, vel = _randin_mushroom(l, w, r)
    return MagneticParticle(pos, vel, zero(SV{T}), T(ω))
end
function _randin_mushroom(l::T = 1.0, w::T = 0.2, r::T = 1.0) where {T <: AbstractFloat}
    x = (-r, r, 2r)
    y = (0, l + r, l + r)
    pos = SV{T}(rand(T)*x[3] + x[1], rand(T)*y[3] + y[1])
    while ! insidemushroom(pos, l, w, r)
        pos = SV{T}(rand(T)*x[3] + x[1], rand(T)*y[3] + y[1])
    end

    φ = T(2π * rand(T))
    vel = SV{T}(cossin(φ)...)

    return pos, vel
end

"""
    randomchaotic(l, w, r)
Generate a chaotic particle, i.e. not trapped in the cap.
"""
function randomchaotic(l, w, r)
    p = randin_mushroom(l, w, r)
    while is_regular(p, l, w, r)
        p = randin_mushroom(l, w, r)
    end
    return p
end

"""
    randomregular(l, w, r)
Generate a regular particle (i.e. trapped in the cap).
"""
function randomregular(l, w, r)
    p = randin_mushroom(l, w, r)
    while !is_regular(p, l, w, r)
        p = randin_mushroom(l, w, r)
    end
    return p
end

end # module

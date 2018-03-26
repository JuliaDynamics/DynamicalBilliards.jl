using PyPlot
export plot_psos

"""
    function plot_psos(ps::Vector{<:AbstractParticle}, bt::BilliardTable, t)

Plots the Poincaré surface of section (see [`psos`](@ref)) in boundary coordinates
by evolving the given particles `ps` for  `t` collisions if `t` is an integer or for
`t` units of time else.

    function plot_psos(n::Int, bt::BilliardTable, t)

Generates `n` random `Particles` inside `bt` and plots their PSOS

    function plot_psos(n::Int, ω::AbstractFloat, bt::BilliardTable, t)

Generates `n` random `MagneticParticles` with cyclotron frequency ω inside `bt` and
 plots their PSOS

"""
function plot_psos(ps::Vector{<:AbstractParticle{T}}, bt::BilliardTable{T}, t) where {T}
    ax = gca()

    ξs,φs,ints = psos(ps,bt,t)

    ax[:plot](ξs, φs, marker="o", ms = 1, color = "b", linestyle="None", mew=0.0,alpha=0.1)

    xmax = 0.0
    for intv ∈ ints
        for xval ∈ intv
            ax[:plot]([xval,xval], [-π/2, π/2], linewidth = 1.5, color = "C0")
            xmax = (xval > xmax)?xval:xmax
        end
    end
    ax[:set_xlim](0,xmax)
    ax[:set_ylim](-π/2,π/2)
    ax[:set_xlabel](L"arc length parameter $\xi$")
    ax[:set_ylabel](L"angle of incidence $\phi$")
end

plot_psos(n::Int, bt::BilliardTable, t) = plot_psos([randominside(bt) for i ∈ 1:n], bt, t)
plot_psos(n::Int, ω::AbstractFloat, bt::BilliardTable, t) = plot_psos([randominside(bt,ω) for i ∈ 1:n], bt, t)

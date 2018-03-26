using PyPlot
export plot_psos

"""
    plot_psos(ξs, φs, ints)

Plots the Poincaré surface of section (see [`psos`](@ref)) in boundary coordinates.

Do `ξs,φs,ints = poincaresection(...)` for the arguments.
"""
function plot_psos(ξs,φs,ints)
    ax = gca()

    ax[:plot](ξs, φs, marker="o", ms = 2.0, color = "b", linestyle="None",
    mew=0.0,alpha=0.1)

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

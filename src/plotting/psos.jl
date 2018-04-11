using PyPlot
export plot_psos

"""
    plot_psos(ξs, φs, intervals; kwargs...)

Plots the Poincaré surface of section in boundary coordinates.
The input arguments are the return values of `boundarymap`.

## Keyword Arguments
* `ax = PyPlot.gca()` : The axis to plot on.
* `color = "C0"` : The color to use for the plotted points. Can be either a
  color for `PyPlot.plot` or a vector of colors of length `length(ξs)`, in
  order to give each initial condition a different color.
* `ms = 1.0` : Marker size of the points.
* `bordercolor = "C3"` : The color of the vertical lines that denote the obstacle
  borders.
* Any other keyword argument is passed to `PyPlot.plot` which plots the points of
  the section.
"""
function plot_psos(ξs, φs, intervals; ax = PyPlot.gca(),
    color = "C0", bordercolor = "C3", ms = 1.0, kwargs...)

    # Plot PSOS
    for (i, (xis, phis)) in enumerate(zip(ξs, φs))
        c = typeof(color) <: AbstractVector ? color[i] : color

        ax[:plot](xis, phis; marker="o", color = c,
        linestyle="None", ms = ms, kwargs...)
    end

    # Plot obstacle limits
    xmax = 0.0
    for intv ∈ intervals
        for xval ∈ intv
            ax[:plot]([xval,xval], [-π/2, π/2], linewidth = 1.5, color = bordercolor,
            alpha = 0.5)
            xmax = (xval > xmax) ? xval : xmax
        end
    end
    ax[:set_xlim](0,xmax)
    ax[:set_ylim](-π/2,π/2)
    ax[:set_xlabel](L"arc length $\xi$")
    ax[:set_ylabel](L"angle of incidence $\phi$")
    return ax
end

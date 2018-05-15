using PyPlot
export plot_boundarymap

"""
    plot_boundarymap(ξs, sφs, intervals; kwargs...)

Plots the boundary map.
The input arguments are the return values of `boundarymap`.

## Keyword Arguments
* `ax = PyPlot.gca()` : The axis to plot on.
* `color = "C0"` : The color to use for the plotted points. Can be either a
  color for `PyPlot.plot` or a vector of colors of length `length(ξs)`, in
  order to give each initial condition a different color.
* `ms = 1.0` : Marker size of the points.
* `bordercolor = "C3"` : The color of the vertical lines that denote the obstacle
  borders.
* `obstacleindices = true`: show obstacle indices above plot
* Any other keyword argument is passed to `PyPlot.plot` which plots the points of
  the section.
"""
function plot_boundarymap(ξs, sφs, intervals; ax = PyPlot.gca(),
    color = "C0", bordercolor = "C3", ms = 1.0, obstacleindices = true, kwargs...)

    # Plot PSOS
    for (i, (xis, sphis)) in enumerate(zip(ξs, sφs))
        c = typeof(color) <: AbstractVector ? color[i] : color

        ax[:plot](xis, sphis; marker="o", color = c,
        linestyle="None", ms = ms, kwargs...)
    end

    # Plot obstacle limits
    xmax = 0.0
    for intv ∈ intervals
        for xval ∈ intv
            ax[:plot]([xval,xval], [-1, 1], linewidth = 1.5, color = bordercolor,
            alpha = 0.5)
            xmax = (xval > xmax) ? xval : xmax
        end
    end
    ax[:set_xlim](0,xmax)
    ax[:set_ylim](-1,1)
    ax[:set_xlabel](L"$\xi$")
    ax[:set_ylabel](L"$\sin(\phi)$")

    #number obstacles by index
    if obstacleindices
        #introduce twin axis
        ax2 = ax[:twiny]()

        #non-labelled major tics at every obstacle border
        ax2[:xaxis][:set_major_formatter](matplotlib[:ticker][:NullFormatter]())
        ax2[:set_xticks]([0,map(x->x[2],intervals)...])

        #zero-length minor tics in between borders, labelled with the obstacle index
        ax2[:xaxis][:set_minor_locator](matplotlib[:ticker][:FixedLocator](map(x->mean(x), intervals)))
        ax2[:xaxis][:set_minor_formatter](matplotlib[:ticker][:FixedFormatter](map(x->dec(x),1:length(intervals))))
        ax2[:tick_params](axis="x", which="minor", length=0)#

        ax2[:set_xlabel]("obstacle index")

        return (ax, ax2)
    end

    return ax
end

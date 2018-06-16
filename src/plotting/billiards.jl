using PyPlot, StaticArrays, InteractiveUtils
export plot_billiard

"""
    translation(obst::Obstacle, vector)
Create a copy of the given obstacle with its position
translated by `vector`.
"""
function translation end

for T in subtypes(Circular)
  @eval translation(d::$T, vec) = ($T)(d.c .+ vec, d.r)
end

for T in subtypes(Wall)
  @eval translation(w::$T, vec) = ($T)(w.sp + vec, w.ep + vec, w.normal)
end

periodicerror() = throw(ArgumentError(
"The billiard must be periodic, i.e. has at least two `PeriodicWall` obstacles"
))

function nonperiodic(bd::Billiard)
    toplot = Obstacle{eltype(bd)}[]
    for obst in bd
        typeof(obst) <: PeriodicWall && continue
        push!(toplot, obst)
    end
    return toplot
end


"""
    plot_billiard(bd::Billiard; ax = (figure(); gca()))
Plot all obstacles in `bd` using the default arguments, set
`xlim` and `ylim` to be 20% larger than `cellsize` and
set the axis aspect ratio to equal.

    plot_billiard(bd, xmin, ymin, xmax, ymax; hexagonal = false, ax = (figure(); gca()))
Plot the given **periodic** billiard `bd`, repeatedly
plotting from `(xmin, ymin)` to `(xmax, ymax)`.
Works for either rectangular periodic billiards, or hexagonal ones. Use
keyword `hexagonal` to denote which one you want.

    plot_billiard(bd, xt::Vector, yt::Vector; hexagonal = false, ax = (figure(); gca()))
Plot the given **periodic** billiard `bd` using the limits defined
by `xt` and `yt`.

Set the keyword argument `plot_orbit = false` to not
plot the orbit defined by `(xt, yt)` and only use the limits.
"""
function plot_billiard(bd::Billiard{T}; ax = (figure(); gca())) where {T}
    sca(ax)
    for obst in bd; plot_obstacle!(obst); end
    xmin, ymin, xmax, ymax = cellsize(bd)
    dx = xmax - xmin; dy = ymax - ymin
    xlim(xmin - 0.1dx, xmax + 0.1dx)
    ylim(ymin - 0.1dy, ymax + 0.1dy)
    ax[:set_aspect]("equal")
end

function plot_billiard(bd::Billiard, xmin, ymin, xmax, ymax;
    hexagonal = false, ax = (figure(); gca()))

    isperiodic(bd) || periodicerror()
    sca(ax)

    if hexagonal
        plot_periodic_hexagon(bd, xmin, ymin, xmax, ymax)
    else
        plot_periodic_rectangle(bd, xmin, ymin, xmax, ymax)
    end
end


function plot_periodic_rectangle(bd, xmin, ymin, xmax, ymax)

    # Cell limits:
    cellxmin, cellymin, cellxmax, cellymax = cellsize(bd)
    dcx = cellxmax - cellxmin
    dcy = cellymax - cellymin

    # Find displacement vectors
    dx = (floor((xmin - cellxmin)/dcx):1:ceil((xmax - cellxmax)/dcx)) * dcx
    dy = (floor((ymin - cellymin)/dcy):1:ceil((ymax - cellymax)/dcy)) * dcy

    # Plot displaced Obstacles
    toplot = nonperiodic(bt)
    for x in dx
        for y in dy
            disp = SVector(x,y)
            for obst in toplot
                plot_obstacle!(translation(obst, disp))
            end
        end
    end
    # Set limits etc.
    PyPlot.xlim(xmin, xmax)
    PyPlot.ylim(ymin, ymax)
    PyPlot.gca()[:set_aspect]("equal")
end

function plot_billiard(bd, xt::AbstractVector, yt::AbstractVector;
    hexagonal = false, ax = (figure(); gca()))

    xmin = floor(minimum(round.(xt,3))); xmax = ceil(maximum(round.(xt,3)))
    ymin = floor(minimum(round.(yt,3))); ymax = ceil(maximum(round.(yt,3)))

    plot_billiard(bd, xmin, ymin, xmax, ymax; hexagonal = hexagonal, ax = ax)

    if plot_orbit
        sca(ax)
        plot(xt, yt, color = "blue")
    end
end

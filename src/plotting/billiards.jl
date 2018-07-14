using InteractiveUtils
export plot_billiard

function nonperiodic(bd::Billiard)
    toplot = Obstacle{eltype(bd)}[]
    for obst in bd
        typeof(obst) <: PeriodicWall && continue
        push!(toplot, obst)
    end
    return toplot
end

periodicerror() = throw(ArgumentError(
"The billiard must be periodic, i.e. has at least two `PeriodicWall` obstacles"
))


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
function plot_billiard(bd::Billiard{T};
    ax = (PyPlot.figure(); PyPlot.gca())) where {T}
    PyPlot.sca(ax)
    for obst in bd; plot_obstacle(obst); end
    xmin, ymin, xmax, ymax = cellsize(bd)
    dx = xmax - xmin; dy = ymax - ymin
    ax[:set_aspect]("equal")
    if !isinf(xmin) && !isinf(xmax)
        PyPlot.xlim(xmin - 0.1dx, xmax + 0.1dx)
    end
    if !isinf(ymin) && !isinf(ymax)
        PyPlot.ylim(ymin - 0.1dy, ymax + 0.1dy)
    end
    return nothing
end

function plot_billiard(bd::Billiard, xmin, ymin, xmax, ymax;
    hexagonal = false, ax = (PyPlot.figure(); PyPlot.gca()))

    isperiodic(bd) || periodicerror()
    sca(ax)

    n = count(x -> typeof(x) <: PeriodicWall, bd)
    if hexagonal
        n != 6 && throw(ArgumentError("Hexagonally periodic billiards have "*
        "exactly 6 periodic walls arranged as a perfect hexagon. Use the "*
        "function `billiard_polygon(6, r, R, setting = \"periodic\")`."))
        plot_periodic_hexagon(bd, xmin, ymin, xmax, ymax)
    else
        n ∉ (2, 4) && throw(ArgumentError("Rectangular periodic billiards must have "*
        "exactly 2 or 4 periodic walls."))
        plot_periodic_rectangle(bd, xmin, ymin, xmax, ymax)
    end

    PyPlot.xlim(xmin, xmax)
    PyPlot.ylim(ymin, ymax)
    return nothing
end

function plot_billiard(bd, xt::AbstractVector, yt::AbstractVector;
    hexagonal = false, ax = (PyPlot.figure(); PyPlot.gca()), plot_orbit = true)

    xmin = minimum(xt); xmax = maximum(xt)
    ymin = minimum(yt); ymax = maximum(yt)

    plot_billiard(bd, xmin, ymin, xmax, ymax; hexagonal = hexagonal, ax = ax)

    if plot_orbit
        PyPlot.sca(ax)
        PyPlot.scatter(xt[1], yt[1], color = "black")
        PyPlot.plot(xt, yt, color = "C0")
    end

    cellxmin, cellymin, cellxmax, cellymax = cellsize(bd)
    dcx = cellxmax - cellxmin
    dcy = cellymax - cellymin
    PyPlot.xlim(xmin - dcx/2, xmax + dcx/2)
    PyPlot.ylim(ymin - dcy/2, ymax + dcy/2)
    PyPlot.gca()[:set_aspect]("equal")

    return nothing
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
    toplot = nonperiodic(bd)
    for x in dx
        for y in dy
            disp = SVector(x,y)
            for obst in toplot
                plot_obstacle(translate(obst, disp))
            end
        end
    end
end

function plot_periodic_hexagon(bd, xmin, ymin, xmax, ymax)

    v = 1
    while !(typeof(bd[v]) <: PeriodicWall); v += 1; end
    space = norm(bd[v].sp - bd[v].ep)*√3

    basis_a = space*SVector(0.0, 1.0)
    basis_b = space*SVector(√3/2, 1/2)
    basis_c = space*SVector(√3, 0.0)

    PyPlot.xlim(xmin - space/2, xmax + space/2)
    PyPlot.ylim(ymin, ymax)
    PyPlot.gca()[:set_aspect]("equal")

    # Cell limits:
    cellxmin, cellymin, cellxmax, cellymax = cellsize(bd)
    dcx = cellxmax - cellxmin
    dcy = cellymax - cellymin
    jmin = Int((ymin - cellymin - dcy/2)÷space) - 1
    jmax = Int((ymax - cellymax + dcy/2)÷space) + 1
    imin = Int((xmin - cellxmin - dcx/2)÷(√3*space)) - 1
    imax = Int((xmax - cellxmax + dcx/2)÷(√3*space)) + 1

    obstacles = nonperiodic(bd)
    for d in obstacles
        for j ∈ jmin:jmax
            for i ∈ imin:imax
                plot_obstacle(translate(d, j*basis_a + i*basis_c))
                plot_obstacle(translate(d, j*basis_a + i*basis_c + basis_b))
            end
        end
    end
end

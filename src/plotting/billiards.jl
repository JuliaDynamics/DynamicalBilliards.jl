using PyPlot, StaticArrays
export plot_billiard, billiard_julia


"""
```julia
plot_billiard(bt::Billiard)
```
Plot all obstacles in `bt` using the default arguments, set
`xlim` and `ylim` to be 10% larger than `cellsize` and
set the axis aspect ratio to equal.

```julia
plot_billiard(bt, xmin, ymin, xmax, ymax)
```
Plot the given (periodic) billiard `bt` on the current PyPlot figure, repeatedly
plotting from `(xmin, ymin)` to `(xmax, ymax)`. Only works for rectangular billiards.

```julia
plot_billiard(bt, xt::Vector, yt::Vector; plot_orbit = true)
```
Plot the given (periodic) billiard `bt` along with a particle trajectory defined
by `xt` and `yt`, on the current PyPlot figure. Only works for rectangular billiards.

Sets limits automatically. Set the keyword argument `plot_orbit = false` to not
plot the orbit defined by `(xt, yt)`.
"""
function plot_billiard(bt::Billiard{T}) where {T}
  for obst in bt
    plot_obstacle(obst)
  end
  xmin, ymin, xmax, ymax = cellsize(bt)
  dx = xmax - xmin; dy = ymax - ymin
  PyPlot.xlim(xmin - 0.1dx, xmax + 0.1dx)
  PyPlot.ylim(ymin - 0.1dy, ymax + 0.1dy)
  PyPlot.gca()[:set_aspect]("equal")
end


function plot_billiard(bt, xmin, ymin, xmax, ymax)
  # Cell limits:
  cellxmin, cellymin, cellxmax, cellymax = cellsize(bt)
  dcx = cellxmax - cellxmin
  dcy = cellymax - cellymin
  # Obstacles to plot:
  toplot = Obstacle{eltype(bt)}[]
  for obst in bt
    typeof(obst) <: PeriodicWall && continue
    push!(toplot, obst)
  end
  # Find displacement vectors (they will multiply dcx, dcy)
  dx = (floor((xmin - cellxmin)/dcx):1:ceil((xmax - cellxmax)/dcx))*dcx
  dy = (floor((ymin - cellymin)/dcy):1:ceil((ymax - cellymax)/dcy))*dcy
  # Plot displaced Obstacles
  for x in dx
    for y in dy
      disp = SVector(x,y)
      for obst in toplot
        plot_obstacle(translation(obst, disp))
      end
    end
  end
  # Set limits etc.
  PyPlot.xlim(xmin, xmax)
  PyPlot.ylim(ymin, ymax)
  PyPlot.gca()[:set_aspect]("equal")
end

function plot_billiard(bt, xt::AbstractVector, yt::AbstractVector; plot_orbit = true)
  xmin = floor(minimum(round.(xt,3))); xmax = ceil(maximum(round.(xt,3)))
  ymin = floor(minimum(round.(yt,3))); ymax = ceil(maximum(round.(yt,3)))
  if plot_orbit
    plot(xt, yt, color = "blue")
  end

  plot_billiard(bt, xmin, ymin, xmax, ymax)
end


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


"""
```julia
billiard_julia(; plotit = true)
```
Return the awesome "Julia-logo" billiard shown in the introduction
of DynamicalBilliards.jl.

By default it also plots the billiard in a new `PyPlot.figure()` using the correct colors.
"""
function billiard_julia(; plotit = true)

  bt = billiard_rectangle()

  r = 0.165
  ewidth = 6.0
  redcent = [0.28, 0.32]
  red = Disk(redcent, r, "Red dot")
  purple = Disk([1 - redcent[1], redcent[2]], r, "Purple dot")
  green = Disk([0.5, 1 - redcent[2]], r, "Green dot")
  push!(bt, red, purple, green)

  if plotit == true
    PyPlot.figure()
    for w in bt
      plot_obstacle(w; color = (0,0,0, 1), linewidth = 3.0)
    end
    plot_obstacle(red; edgecolor = (203/255, 60/255, 51/255),
    facecolor = (213/255, 99/255, 92/255), linewidth = ewidth)
    plot_obstacle(purple; edgecolor = (149/255, 88/255, 178/255),
    facecolor = (170/255, 121/255, 193/255), linewidth = ewidth)
    plot_obstacle(green, edgecolor = (56/255, 152/255, 38/255),
    facecolor = (96/255, 173/255, 81/255), linewidth = ewidth)

    # particle edge color
    # darkblue = (64/255, 99/255, 216/255)
    # lightblue = (102/255, 130/255, 223/255)

    PyPlot.axis("off")
    PyPlot.tight_layout()
    PyPlot.gca()[:set_aspect]("equal")
    PyPlot.xlim(-0.1,1.1)
    PyPlot.ylim(-0.1,1.1)
  end

  return bt
end

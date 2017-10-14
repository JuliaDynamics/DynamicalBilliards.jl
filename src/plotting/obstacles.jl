using PyPlot, StaticArrays
export plot_obstacle

"""
```julia
plot_obstacle(obst::Obstacle; kwargs...)
```
Plot given obstacle on the current `PyPlot` figure.

The default arguments for each type of obstacle have been chosen for maximum clarity and
consistency.

The `kwargs...` given by the user are keywords passed directly into PyPlot's
constructors. For `Wall` obstacles, kwargs are passed into `PyPlot.plot()`. For
`Disk` obstacles, kwargs are passed into `PyPlot.plt[:Circle]()`.
"""
function plot_obstacle(d::Disk; kwargs...)
  circle1 = PyPlot.plt[:Circle](d.c, d.r;
  edgecolor = (0,0.6,0), facecolor = (0, 0.6,0, 0.5), kwargs...)
  PyPlot.gca()[:add_artist](circle1)
  PyPlot.show()
end

function plot_obstacle(d::RandomDisk; kwargs...)
  circle1 = PyPlot.plt[:Circle](d.c, d.r;
  edgecolor = (0.8,0.8,0), facecolor = (0.8, 0.8, 0, 0.5), linewidth = 2.0, kwargs...)
  PyPlot.gca()[:add_artist](circle1)
  PyPlot.show()
end

function plot_obstacle(d::Antidot; kwargs...)
  circle1 = PyPlot.plt[:Circle](d.c, d.r;
  edgecolor = (0.8,0.0,0), linewidth = 2.0, facecolor = (0.6, 0.0, 0, 0.1),
  linestyle="dashed", kwargs...)
  PyPlot.gca()[:add_artist](circle1)
  PyPlot.show()
end

function plot_obstacle(w::Wall; kwargs...)
  PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
  color=(0,0.6,0), linewidth = 2.0, ms=0, kwargs...)
  PyPlot.show()
end

function plot_obstacle(w::FiniteWall; kwargs...)
  if w.isdoor
    PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
    color="black", linestyle = "-", linewidth = 2.0, ms=0, kwargs...)
    PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
    color=(0, 0.9, 0.9), linestyle = "--", linewidth = 2.0, ms=0, kwargs...)
  else
    PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
    color=(0,0.6,0), linewidth = 2.0, ms=0, kwargs...)
  end
  PyPlot.show()
end

function plot_obstacle(w::RandomWall; kwargs...)
  PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
  color=(0.8,0.8,0), linewidth = 2.0, ms=0, kwargs...)
  PyPlot.show()
end

function plot_obstacle(w::SplitterWall; kwargs...)
  PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
  color=(0.8,0.0,0), linewidth = 3.0, linestyle="dashed", kwargs...)
  PyPlot.show()
end

function plot_obstacle(w::PeriodicWall; kwargs...)
  PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
  color="purple", linewidth = 1.0, linestyle="dotted", alpha = 0.5, kwargs...)
  PyPlot.show()
end

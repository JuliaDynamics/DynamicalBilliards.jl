export plot_cyclotron

"""
    plot_cyclotron(p::MagneticParticle; use_cell=true, kwargs...)
Plot the circle traced by the free particle motion. Optionally use `p.current_cell` for
the particle's position. The user provided `kwargs...`
are passed onto `PyPlot.plt[:Circle]`.
"""
function plot_cyclotron(p::MagneticParticle; use_cell=true, kwargs...)
  ω = p.omega
  pos = use_cell ? p.pos + p.current_cell : p.pos
  c = pos - (1/ω)*[p.vel[2], -p.vel[1]]
  r = abs(1/ω)
  circle1 = PyPlot.plt[:Circle](c, r,
  edgecolor = (0.0,0.0,0.8, 0.5), linewidth = 2.0, facecolor = (0., 0.0, 0.8, 0.05),
  kwargs...)
  PyPlot.gca()[:add_artist](circle1)
  PyPlot.show()
end

"""
    plot(p::AbstractParticle; use_cell=true, kwargs...)
Plot given particle on the current `PyPlot` axes. Optionally use `p.current_cell` for
the particle's position. Given `kwargs...` are passed onto `PyPlot.scatter`.

The particle is represented as a small ball (`PyPlot.scatter`) and a small arrow
(`PyPlot.quiver`).
All `kwargs...` are given to `scatter` but if a keyword argument `color` is given,
it is also passed to `quiver`.
"""
function plot(p::AbstractParticle) end

function plot(p::AbstractParticle{T}; use_cell=true, kwargs...) where {T}
  pos = use_cell ? p.pos + p.current_cell : p.pos
  kwargs = Dict(kwargs)
  # Set same color for arrow and point:
  if haskey(kwargs, :color)
    c = kwargs[:color]
  else
    c = (0,0,0)
  end
  # Plot position:
  s1 = PyPlot.scatter(pos...; color=c, s= 30.0, kwargs...)
  # Plot velocity:
  q1 = PyPlot.quiver(pos..., 0.08p.vel...; angles = "xy", scale = 1, width = 0.005, color=c)
  return s1, q1
end

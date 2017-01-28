using PyPlot

export plot_obstacle, plot_particle, plot_billiard, plot_cyclotron,
animate_evolution
####################################################
## Plot Billiards
####################################################
"""
```julia
plot_obstacle(obst::Obstacle; kwargs...)
```
Plot given obstacle on the current `PyPlot` figure.

The default arguments for each type of obstacle have been chosen for maximum clarity and
visual beauty.

The `kwargs...` given by the user are keywords passed directly into PyPlot's
constructors (like e.g. linewidth = 2.0). Be careful to only pass keywords
that are accepted by `PyPlot`.
"""
function plot_obstacle(d::Disk; kwargs...)
  circle1 = PyPlot.plt[:Circle](d.c, d.r;
  edgecolor = (0,0.6,0), facecolor = (0, 0.6,0, 0.5), kwargs...)
  PyPlot.gca()[:add_artist](circle1)
  PyPlot.show()
end


function plot_obstacle(d::Antidot; kwargs...)
  circle1 = PyPlot.plt[:Circle](d.c, d.r;
  edgecolor = (0.8,0.0,0), linewidth = 2.0, facecolor = (0.6, 0.0, 0, 0.1),
  kwargs...)
  PyPlot.gca()[:add_artist](circle1)
  PyPlot.show()
end

function plot_obstacle(w::FiniteWall; kwargs...)
  PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
  color=(0,0.6,0), linewidth = 2.0, ms=0, kwargs...)
  PyPlot.show()
end

function plot_obstacle(w::SplitterWall; kwargs...)
  PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
  color=(0.8,0.0,0), linewidth = 3.0, kwargs...)
  PyPlot.show()
end

function plot_obstacle(w::PeriodicWall; kwargs...)
  PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]];
  color="purple", linewidth = 1.0, alpha = 0.5, linestyle="dashed", kwargs...)
  PyPlot.show()
end

"""
```julia
plot_billiard(bt::Vector{Obstacle})
```

Plot all obstacles in `bt` using the default arguments, set
`xlim` and `ylim` to be 10% larger than `cellsize` and
set the axis aspect ration to equal.
"""
function plot_billiard(bt::Vector{Obstacle})
  for obst in bt
    plot_obstacle(obst)
  end
  xmin, ymin, xmax, ymax = cellsize(bt)
  dx = xmax - xmin; dy = ymax - ymin
  PyPlot.xlim(xmin - 0.1dx, xmax + 0.1dx)
  PyPlot.ylim(ymin - 0.1dy, ymax + 0.1dy)
  PyPlot.gca()[:set_aspect]("equal")
end


####################################################
## Plot Particle
####################################################
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
```julia
plot_particle(p::AbstractParticle; use_cell=true, kwargs...)
```

Plot given particle on the current `PyPlot` figure. Optionally use `p.current_cell` for
the particle's position.

The particle is represented as a small ball (`scatter()`) and a small arrow (`quiver()`).
The user provided `kwargs...` are passed only onto the `scatter()` call. However,
if a keyword argument `color` is given, it is also passed to `quiver()`.
"""
function plot_particle(p::AbstractParticle; use_cell=true, kwargs...)
  pos = use_cell ? p.pos + p.current_cell : p.pos
  kwargs = Dict(kwargs)
  # Set same color for arrow and point:
  if haskey(kwargs, :color)
    c = kwargs[:color]
  else
    c = (0,0,0)
  end
  # Plot position:
  s1 = scatter(pos...; color=c, s= 20.0, kwargs...)
  # Plot velocity:
  q1 = quiver(pos..., p.vel...,
  units="dots", angles = "xy", scale = 0.02, width = 2.0, color=c)
  return s1, q1
end

####################################################
## Animate Particle
####################################################

"""
```julia
animate_evolution(p, bt, colnumber[, ray-splitter];
sleeptime = 0.1, col_to_plot = 5, orbit_color = (0,0,1), savefigs = false, savename = "")
```

Animate the evolution of the particle, plotting the orbit from collision to collision.

Notice the difference with `evolve!()`: No time is given here; instead a number of
collisions is passed.

### Arguments
* `p::AbstractParticle` : Either standard or magnetic.
* `bt::Vector{Obstacle}` : The billiard table.
* `colnumber::Int` : Number of collisions to evolve the particle for.
* `ray-splitter::Dict{Int, Vector{Function}}` : (Optional) Ray-splitting dictionary
  that enables ray-splitting processes during evolution.
### Keyword Arguments
* `sleeptime` : Time passed to `sleep()` between each collision.
* `col_to_plot` : How many previous collisions are shown during the animation.
* `savefigs` : Save .png figures to enable the creation of animation afterwards.
  **WARNING:** currently the .gif production has to be made by the user!
* `savename` : Name (*including path!*) of the figures to be produced. The ending
  "_ i.png" will be attached to all.
"""
function animate_evolution(p::MagneticParticle, bt, colnumber;
  sleeptime = 0.1, col_to_plot = 5, orbit_color = (0,0,1), savefigs = false, savename = "")

  sleeptime == 0 && (sleeptime = 1e-6)
  ω = p.omega
  ε = eps()
  i=0
  xdata = Vector{Float64}[]
  ydata = Vector{Float64}[]

  while i < colnumber

    xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, ε)...)

    if i < col_to_plot
      push!(xdata, xt)
      push!(ydata, yt)
    else
      shift!(xdata); shift!(ydata)
      push!(xdata, xt); push!(ydata, yt)
    end

    xpd = Float64[]
    for el in xdata; append!(xpd, el); end

    ypd = Float64[]
    for el in ydata; append!(ypd, el); end

    if i == 0
      line, = plot(xpd, ypd, color = orbit_color)
    end
    line[:set_xdata](xpd)
    line[:set_ydata](ypd)
    point, quiv = plot_particle(p)
    if savefigs
      s = savename*"_$(i+1).png"
      savefig(s, dpi = 60, bbox_inches="tight")
    end


    sleep(sleeptime)
    if i < colnumber - 1
      point[:remove]()
      quiv[:remove]()
    end
    i+=1
  end
end

# Magnetic + Ray-splitting
function animate_evolution(p::MagneticParticle, bt,
  colnumber, rayspl::Dict{Int, Vector{Function}};
  sleeptime = 0.1, col_to_plot = 5, color = (0,0,1), savefigs = false, savename = "")

  sleeptime == 0 && (sleeptime = 1e-6)
  ω = p.omega
  ε = eps()
  i=0
  xdata = Vector{Float64}[]
  ydata = Vector{Float64}[]

  while i < colnumber

    xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, ε, rayspl)...)

    if i < col_to_plot
      push!(xdata, xt)
      push!(ydata, yt)
    else
      shift!(xdata); shift!(ydata)
      push!(xdata, xt); push!(ydata, yt)
    end

    xpd = Float64[]
    for el in xdata; append!(xpd, el); end

    ypd = Float64[]
    for el in ydata; append!(ypd, el); end

    if i == 0
      line, = plot(xpd, ypd, color = orbit_color)
    end
    line[:set_xdata](xpd)
    line[:set_ydata](ypd)
    point, quiv = plot_particle(p)
    if savefigs
      s = savename*"_$(i+1).png"
      savefig(s, dpi = 60, bbox_inches="tight")
    end


    sleep(sleeptime)

    if i < colnumber - 1
      point[:remove]()
      quiv[:remove]()
    end
    i+=1
  end
end

# Straight
function animate_evolution(p::Particle, bt, colnumber;
  sleeptime = 0.1, col_to_plot = 5, orbit_color = (0,0,1), savefigs = false, savename = "")

  sleeptime == 0 && (sleeptime = 1e-6)
  ε = eps()
  i=0
  xdata = Vector{Float64}[]
  ydata = Vector{Float64}[]

  while i < colnumber

    xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, ε)...)

    if i < col_to_plot
      push!(xdata, xt)
      push!(ydata, yt)
    else
      shift!(xdata); shift!(ydata)
      push!(xdata, xt); push!(ydata, yt)
    end

    xpd = Float64[]
    for el in xdata; append!(xpd, el); end

    ypd = Float64[]
    for el in ydata; append!(ypd, el); end

    if i == 0
      line, = plot(xpd, ypd, color = orbit_color)
    end
    line[:set_xdata](xpd)
    line[:set_ydata](ypd)
    point, quiv = plot_particle(p)
    if savefigs
      s = savename*"_$(i+1).png"
      savefig(s, dpi = 60, bbox_inches="tight")
    end


    sleep(sleeptime)
    if i < colnumber - 1
      point[:remove]()
      quiv[:remove]()
    end
    i+=1
  end
end


# Straight + Ray-splitting
function animate_evolution(p::Particle, bt, colnumber, rayspl::Dict{Int, Vector{Function}};
  sleeptime = 0.1, col_to_plot = 5, color = (0,0,1), savefigs = false, savename = "")

  sleeptime == 0 && (sleeptime = 1e-6)
  ε = eps()
  i=0
  xdata = Vector{Float64}[]
  ydata = Vector{Float64}[]

  while i < colnumber

    xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, ε, rayspl)...)

    if i < col_to_plot
      push!(xdata, xt)
      push!(ydata, yt)
    else
      shift!(xdata); shift!(ydata)
      push!(xdata, xt); push!(ydata, yt)
    end

    xpd = Float64[]
    for el in xdata; append!(xpd, el); end

    ypd = Float64[]
    for el in ydata; append!(ypd, el); end

    if i == 0
      line, = plot(xpd, ypd, color = color)
    end
    line[:set_xdata](xpd)
    line[:set_ydata](ypd)
    point, quiv = plot_particle(p)
    if savefigs
      s = savename*"_$(i+1).png"
      savefig(s, dpi = 60, bbox_inches="tight")
    end


    sleep(sleeptime)
    if i < colnumber - 1
      point[:remove]()
      quiv[:remove]()
    end
    i+=1
  end
end

using PyPlot

export plot_obstacle, plot_particle, plot_billiard, plot_cyclotron,
plot_evolution
####################################################
## Plot Billiards
####################################################
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
  PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]],
  color=(0,0,0), linewidth = 2.0, ms=0)
  PyPlot.show()
end

function plot_obstacle(w::SplitterWall; kwargs...)
  PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]],
  color=(0.8,0.0,0), linewidth = 3.0)
  PyPlot.show()
end

function plot_obstacle(w::PeriodicWall; kwargs...)
  PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]],
  color="purple", linewidth = 1.0, alpha = 0.5, linestyle="dashed", kwargs...)
  PyPlot.show()
end

function plot_billiard(bt::Vector{Obstacle})
  for obst in bt
    plot_obstacle(obst)
  end
  PyPlot.gca()[:set_aspect]("equal")
end


####################################################
## Plot Particle Evolution
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



function plot_evolution(p::MagneticParticle, bt, colnumber = 50;
  sleeptime = 0.1, col_to_plot = 5, color = (0,0,1), savefigs = false, savename = "")

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

    point[:remove]()
    quiv[:remove]()
    i+=1
  end
end

function plot_evolution(p::Particle, bt, colnumber = 50;
  sleeptime = 0.1, col_to_plot = 5, color = (0,0,1), savefigs = false, savename = "")

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

    point[:remove]()
    quiv[:remove]()
    i+=1
  end
end

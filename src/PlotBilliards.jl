#using PyPlot

export plot_obstacle, plot_particle, plot_billiard, plot_cyclotron,
plot_evolution
####################################################
## Plot Billiards
####################################################
function plot_obstacle(d::Disk; color = "green")
  circle1 = PyPlot.plt[:Circle](d.c, d.r, alpha=0.3, color=color, lw=0.0)
  PyPlot.gca()[:add_artist](circle1)
  PyPlot.show()
end

function plot_obstacle(d::Antidot; color = "red")
  circle1 = PyPlot.plt[:Circle](d.c, d.r, color=color, fill=false, lw=1.0)
  PyPlot.gca()[:add_artist](circle1)
  PyPlot.show()
end

function plot_obstacle(w::FiniteWall; color = "black")
  PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]], color=color, linewidth = 3.0)
  PyPlot.show()
end

function plot_obstacle(w::SplitterWall; color = "red")
  PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]], color=color, linewidth = 2.0)
  PyPlot.show()
end

function plot_obstacle(w::PeriodicWall; color = "purple")
  PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]], color=color,
  linewidth = 2.0, alpha = 0.5, linestyle="dashed")
  PyPlot.show()
end

function plot_billiard(bt::Vector{Obstacle})
  for obst in bt
    plot_obstacle(obst)
  end
  gca()[:set_aspect]("equal")
end


####################################################
## Plot Particle Evolution
####################################################
function plot_cyclotron(p::MagneticParticle; use_cell=true)
  ω = p.omega
  pos = use_cell ? p.pos + p.current_cell : p.pos
  c = pos - (1/ω)*[p.vel[2], -p.vel[1]]
  r = abs(1/ω)
  circle1 = PyPlot.plt[:Circle](c, r, alpha=0.1, color = "red")
  PyPlot.gca()[:add_artist](circle1)
  PyPlot.show()
end

function plot_particle(p::AbstractParticle; color = (0.0, 0.0, 0.0), use_cell=true)
  pos = use_cell ? p.pos + p.current_cell : p.pos
  s1 = scatter(pos..., color=color)
  q1 = quiver(pos..., p.vel...,
  units="dots", angles = "xy", scale = 0.02, width = 2.0, color=color)
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

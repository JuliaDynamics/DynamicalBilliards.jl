using PyPlot

####################################################
## Plot Particle Evolution
####################################################
function plot_obstacle(d::Disk; color = "green")
  circle1 = PyPlot.plt[:Circle](d.c, d.r, alpha=0.3, color=color, lw=0.0)
  PyPlot.gca()[:add_artist](circle1)
  PyPlot.show()
end

function plot_obstacle(d::Circle; color = "black")
  circle1 = PyPlot.plt[:Circle](d.c, d.r, color=color, fill=false, lw=3.0)
  PyPlot.gca()[:add_artist](circle1)
  PyPlot.show()
end

function plot_obstacle(w::FiniteWall; color = "black")
  PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]], color=color, linewidth = 3.0)
  PyPlot.show()
end

function plot_obstacle(w::PeriodicWall; color = "purple")
  PyPlot.plot([w.sp[1],w.ep[1]],[w.sp[2],w.ep[2]], color=color,
  linewidth = 2.0, alpha = 0.5, linestyle="dashed")
  PyPlot.show()
end

function plot_billiard{T<:AbstractFloat}(bt::Vector{Obstacle{T}})
  for obst in bt
    plot_obstacle(obst)
  end
  gca()[:set_aspect]("equal")
end



####################################################
## Plot Particle Evolution
####################################################
function plot_cyclotron(ω::AbstractFloat, p::AbstractParticle; use_cell=true)
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


#REDO
# function plot_evolve!(p::Particle, bt::Vector{Obstacle}, ttotal::AbstractFloat,
#   dt::AbstractFloat = 0.1; which_figure = PyPlot.figure(figsize=(10,10)))
#
#   fig = PyPlot.figure(which_figure[:number])
#   PyPlot.cla()
#   index_t0 = length(p.t) #plot from index_t0 to end
#   evolve!(p, bt, ttotal, dt)
#
#   #first plot billiard
#   for obst in bt
#     plot_obstacle(obst)
#   end
#   #plot particle:
#   plot(p.xoft[index_t0:end], p.yoft[index_t0:end], color="blue")
#   scatter(p.pos..., color="blue")
#   quiver(p.pos..., p.vel..., color="blue")
# end

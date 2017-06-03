using PyPlot
export animate_evolution

"""
```julia
animate_evolution(p, bt, colnumber[, ray-splitter]; kwargs)
```
Animate the evolution of the particle, plotting the orbit from collision to collision.

Notice the difference with `evolve!()`: No time is given here; instead a number of
collisions is passed.

## Arguments
* `p::AbstractParticle` : Either standard or magnetic.
* `bt::Vector{Obstacle}` : The billiard table.
* `colnumber::Int` : Number of collisions to evolve the particle for.
* `ray-splitter::Dict{Int, Any}` : (Optional) Ray-splitting dictionary
  that enables ray-splitting processes during evolution.
## Keyword Arguments
* `sleeptime` : Time passed to `sleep()` between each collision.
* `col_to_plot` : How many previous collisions are shown during the animation.
* `savefigs::Bool` : If `true` save .png figures to enable the creation of animation afterwards.
  (currently the .gif production has to be made by the user!)
* `savename` : Name (**including path!**) of the figures to be produced. The ending
  "\_i.png" will be attached to all figures.
* `particle_kwargs` : Either a Dict{Symbol, Any} or a vector of Tuple{Symbol, Any}.
  Keywords passed into `plot_particle()`.
* `orbit_kwargs` : Either a Dict{Symbol, Any} or a Vector of Tuple{Symbol, Any}.
  Keywords passed into `PyPlot.plot()` which plots the orbit of the particle (`line` object).
"""
function animate_evolution(p::AbstractParticle, bt, colnumber;
  sleeptime = 0.1, col_to_plot = 5, savefigs = false, savename = "",
  particle_kwargs = nothing, orbit_kwargs = nothing)

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
      if orbit_kwargs != nothing
        line, = plot(xpd, ypd; orbit_kwargs...)
      else
        line, = plot(xpd, ypd; color = "blue")
      end
    end
    line[:set_xdata](xpd)
    line[:set_ydata](ypd)

    if particle_kwargs != nothing
      point, quiv = plot_particle(p; particle_kwargs...)
    else
      point, quiv = plot_particle(p)
    end

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
function animate_evolution(p::AbstractParticle, bt,
  colnumber, rayspl::Dict;
  sleeptime = 0.1, col_to_plot = 5, orbit_color = (0,0,1),
  savefigs = false, savename = "", particle_color = (0,0,0),
  particle_kwargs = nothing, orbit_kwargs = nothing)

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
      if orbit_kwargs != nothing
        line, = plot(xpd, ypd; orbit_kwargs...)
      else
        line, = plot(xpd, ypd; color = "blue")
      end
    end
    line[:set_xdata](xpd)
    line[:set_ydata](ypd)

    if particle_kwargs != nothing
      point, quiv = plot_particle(p; particle_kwargs...)
    else
      point, quiv = plot_particle(p)
    end

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

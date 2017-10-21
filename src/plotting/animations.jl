using PyPlot
export animate_evolution!

"""
```julia
animate_evolution!(p, bt, colnumber[, ray-splitter]; kwargs...)
```
Animate the evolution of the particle, plotting the orbit from collision to collision.

### Arguments
  * `p::AbstractParticle` : The particle to be evolved (gets mutated!).
  * `bt::Vector{Obstacle}` : The billiard table.
  * `colnumber::Int` : Number of collisions to evolve the particle for.
  * `ray-splitter::Dict{Int, Any}` : (Optional) Ray-splitting dictionary
      that enables ray-splitting processes during evolution.
### Keyword Arguments
  * `newfig = true` : Creates a new figure at the function call, and plots
    the billiard table in that figure.
  * `sleeptime` : Time passed to `sleep()` between each collision.
  * `col_to_plot` : How many previous collisions are shown during the animation.
  * `particle_kwargs` : Either a Dict{Symbol, Any} or a vector of Tuple{Symbol, Any}.
    Keywords passed into `plot_particle()`.
  * `orbit_kwargs` : Either a Dict{Symbol, Any} or a Vector of Tuple{Symbol, Any}.
    Keywords passed into `PyPlot.plot()` which plots the orbit of the particle
    (`line` object).
  * `savefigs::Bool` : If `true` save .png figures of each frame of the animation
    A direct movie (like creating a .mp4) of the animation cannot be made automatically,
    since the animation process mutates the particle.
  * `savename` : Name (*including path*) of the figures to be produced. The ending
    "\_i.png" will be attached to all figures.
"""
function animate_evolution!(p::AbstractParticle, bt, colnumber, raysplit = Dict();
    sleeptime = 0.1, col_to_plot = 5, savefigs = false, savename = "",
    particle_kwargs = nothing, orbit_kwargs = nothing, newfig = true)

    if newfig == true
        fig = figure()
        plot_billiard(bt)
    end

    sleeptime == 0 && (sleeptime = 1e-3)
    i=0
    xdata = Vector{Float64}[]
    ydata = Vector{Float64}[]

    while i < colnumber

    if raysplit == Dict()
        xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 1)...)
    else
        xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 1, rayspl)...)
    end

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
        savefig(s, bbox_inches="tight")
    end


    sleep(sleeptime)
    if i < colnumber - 1
        point[:remove]()
        quiv[:remove]()
    end
    i+=1
    end
end

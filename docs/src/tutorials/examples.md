This section has some examples of usage of `DynamicalBilliards.jl`, with some brief
comments.

## Julia-logo Billiard
The "Julia-logo" billiard, accessed by `billiard_julia()` simply wraps this code:
```julia
using DynamicalBilliards, DynamicalBilliardsPlotting, PyPlot

bt = Obstacle{Float64}[]
bt = billiard_rectangle()
for w in bt
  plot_obstacle(w; color = (0,0,0, 1), linewidth = 3.0)
end

r = 0.165
ewidth = 6.0
redcent = [0.28, 0.32]
red = Disk(redcent, r, "red")
plot_obstacle(red; edgecolor = (203/255, 60/255, 51/255),
facecolor = (213/255, 99/255, 92/255), linewidth = ewidth)

purple = Disk([1 - redcent[1], redcent[2]], r, "purple")
plot_obstacle(purple; edgecolor = (149/255, 88/255, 178/255),
facecolor = (170/255, 121/255, 193/255), linewidth = ewidth)

green = Disk([0.5, 1 - redcent[2]], r, "green")
plot_obstacle(green, edgecolor = (56/255, 152/255, 38/255),
facecolor = (96/255, 173/255, 81/255), linewidth = ewidth)

push!(bt, red, purple, green)
# particle colors
darkblue = (64/255, 99/255, 216/255)
lightblue = (102/255, 130/255, 223/255)

p = randominside(bt, 2.0)

PyPlot.axis("off")
PyPlot.tight_layout()
PyPlot.gca()[:set_aspect]("equal")
PyPlot.xlim(-0.1,1.1)
PyPlot.ylim(-0.1,1.1)

okwargs = Dict(:linewidth => 2.0, :color => lightblue)
pkwargs = Dict(:color => darkblue, :s => 150.0)

# create and save the animation:
sname = "C:\\***\\anim"
animate_evolution(p, bt, 200; col_to_plot = 7,
particle_kwargs = pkwargs, orbit_kwargs = okwargs,
savefigs = true, savename = sname)
```
and produces:

![Julia-logo billiard animation](http://i.imgur.com/EtKof48.gif)


## Mean Free Path of the Lorentz Gas
```julia
using DynamicalBilliards
bt = billiard_lorentz(0.2) #alias for billiard_sinai(setting = "periodic")
mfp = 0.0
for i in 1:1000
  p = randominside(bt)
  ct, poss, vels = evolve!(p, bt, 10000.0)
  #skip first two entries because they are not "full" collisions:
  mfp += mean(ct[3:end])
end
mfp /= 1000
```
gives the value of `2.1899...` which is very close to the analytic result:

``\text{m.f.p.} =  \frac{1-\pi r^2 }{2r} \stackrel{r=0.2}{=} 2.18584 ``

which you can find for example [here](http://www.cmls.polytechnique.fr/perso/golse/Surveys/FGIcmp03.pdf).

## Semi-Periodic Billiard
`DynamicalBilliards.jl` allows for your system to be periodic in only some specific
directions. For example, the following code produces a billiard that is periodic
in only the x-direction:

```julia
using DynamicalBilliards
o = 0.0; x = 2.0; y=1.0
bt = Obstacle{Float64}[]

sp = [o,o]; ep = [o, y]; n = [x,o]
leftw = PeriodicWall(sp, ep, n, "Left periodic boundary")
sp = [x,o]; ep = [x, y]; n = [-x,o]
rightw = PeriodicWall(sp, ep, n, "Right periodic boundary")

sp = [o,y]; ep = [x, y]; n = [o,-y]
topw2 = FiniteWall(sp, ep, n, "Top wall")
sp = [o,o]; ep = [x, o]; n = [o,y]
botw2 = FiniteWall(sp, ep, n, "Bottom wall")
push!(bt, leftw, rightw, topw2, botw2)

r = 0.25
d = Disk([0.5, 0.5], r)
d2 = Disk([1.5, 0.5], r/2)
push!(bt, d, d2)

p = randominside(bt, 0.5)
xt, yt, vxt, vyt, t = construct(evolve!(p, bt, 25)...)
plot_billiard(bt, xt, yt)
plot_particle(p)
```
Result:

![Semi-Periodic Billiard](http://i.imgur.com/Dbxmq8y.png)

## Ray-Splitting
The following code produces an animation of a Ray-Splitting billiard:
```julia
using DynamicalBilliards, PyPlot

bt = billiard_rectangle(2, 1)
sw = SplitterWall([1.0, 0.0], [1,1], [-1,0], true)
push!(bt, sw)
a1 = Antidot([0.5, 0.5], 0.3)
push!(bt, a1)
a2 = Antidot([1.5, 0.5], 0.2)
push!(bt, a2)

sa = (θ, where, ω) -> where ? 1.25*θ : 0.8*θ
Tp = (p) -> (θ, where, ω) -> begin
  if where
    abs(θ) < π/2/1.25 ? p*exp(-(θ)^2/2(π/8)^2) : 0.0
  else
    (p)*exp(-(θ)^2/2(π/4)^2)
  end
end
newo = ((x, bool) -> bool ? -2.0x : -0.5x)
newo2 = ((x, bool) -> bool ? -x : -x)
rayspl = Dict(
5 => (Tp(0.9), sa, newo2),
6 => (Tp(0.7), sa, newo),
7 => (Tp(0.65), sa, newo))

p = randominside(bt, 1.0)
plot_billiard(bt)

savedir = "C:\\***\\anim"
animate_evolution(p, bt, 200, rayspl, savefigs = true, savename = savedir)
```
result:

![Ray-splitter animation](http://i.imgur.com/89s0fon.gif)

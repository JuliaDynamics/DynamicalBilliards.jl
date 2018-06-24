This section has some examples of usage of `DynamicalBilliards`, with some brief
comments.

## Julia-logo Billiard
The "Julia-logo" billiard animation was made with:
```julia
using DynamicalBilliards, PyPlot
figure()
bt = billiard_rectangle()
# Plot walls:
for w in bt
  plot_obstacle(w; color = (0,0,0, 1), linewidth = 3.0)
end

# Create and plot the 3 disks:
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

# Create billiard
bt = Billiard(bt.obstacles..., red, purple, green)

# Set axis
axis("off")
tight_layout()
gca()[:set_aspect]("equal")
xlim(-0.1,1.1)
ylim(-0.1,1.1)

# Create a particle
p = randominside(bt, 2.0)
# particle colors
darkblue = (64/255, 99/255, 216/255)
lightblue = (102/255, 130/255, 223/255)

okwargs = Dict(:linewidth => 2.0, :color => lightblue)
pkwargs = Dict(:color => darkblue, :s => 150.0)

# create the animation:
animate_evolution(p, bt, 200; col_to_plot = 7,
particle_kwargs = pkwargs, orbit_kwargs = okwargs, newfig = false)
```
and produces:

![Julia-logo billiard animation](http://i.imgur.com/EtKof48.gif)


## Mean Free Path of the Lorentz Gas
```@example tut3
using DynamicalBilliards
bt = billiard_lorentz(0.2) #alias for billiard_sinai(setting = "periodic")
mfp = meancollisiontime(randominside(bt), bt, 1000000.0)
```
The result is very close to the analytic result:

``\text{m.f.p.} =  \frac{1-\pi r^2 }{2r} \stackrel{r=0.2}{=} 2.18584 ``

which you can find for example [here](http://www.cmls.polytechnique.fr/perso/golse/Surveys/FGIcmp03.pdf).

## Semi-Periodic Billiard
`DynamicalBilliards.jl` allows for your system to be periodic in only some specific
directions. For example, the following code produces a billiard that is periodic
in only the x-direction:

```@example tut3
using DynamicalBilliards, PyPlot
o = 0.0; x = 2.0; y=1.0
bt = Obstacle{Float64}[]

sp = [o,o]; ep = [o, y]; n = [x,o]
leftw = PeriodicWall(sp, ep, n, "Left periodic boundary")
sp = [x,o]; ep = [x, y]; n = [-x,o]
rightw = PeriodicWall(sp, ep, n, "Right periodic boundary")

sp = [o,y]; ep = [x, y]; n = [o,-y]
topw2 = InfiniteWall(sp, ep, n, "Top wall")
sp = [o,o]; ep = [x, o]; n = [o,y]
botw2 = InfiniteWall(sp, ep, n, "Bottom wall")

r = 0.25
d = Disk([0.5, 0.5], r)
d2 = Disk([1.5, 0.5], r/2)
push!(bt, d, d2)

bd = Billiard(leftw, rightw, topw2, botw2, d, d2)

p = randominside(bd)
p.pos = [0.311901, 0.740439]
p.vel = [0.548772, 0.835972]
xt, yt, vxt, vyt, t = construct(evolve(p, bd, 25)...)
plot_billiard(bd, xt, yt)
scatter(xt[1], yt[1])
scatter(xt[end], yt[end], color = "black")
ylim(0,y)
savefig("xperiodic.svg"); nothing # hide
```
![](xperiodic.svg)

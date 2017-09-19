![DynamicalBilliards Logo: The Julia billiard](http://i.imgur.com/NKgzYrt.gif)

`DynamicalBilliards.jl` is an easy-to-use, modular and extendable Julia package for
Dynamical Billiard systems in two dimensions.

!!! tip "Logo billiard"
    Checkout the [Julia billiard animation](#julia-billiard-animation)
    code to see how to create the animated billiard of our logo.

A dynamical billiard is a system where a particle is propagating as a straight line from obstacle to obstacle, performing specular reflection at the boundary of the obstacles. Billiard systems have been used extensively in mathematics, nonlinear dynamics and chaos and played an important role in the development of nonlinear science.

The above will be referred to as the "standard" billiard. The [wikipedia page](https://en.wikipedia.org/wiki/Dynamical_billiards) has many examples of different types of billiards. Also, the [scholarpedia](http://www.scholarpedia.org/article/Dynamical_billiards) entry is a good read on the subject. In general, the standard billiard can be extended in many ways.

The types that are currently offered by this package, besides the standard one, are magnetic and ray-splitting billiards. In a magnetic billiard the particle's orbit is a circle (like electrons in a perpendicular magnetic field). In ray-splitting (a.k.a. semiclassical) billiards the particle may propagate *through* an obstacle, given some arbitrary transmission and refraction law.

This package does not support finite-sized particles and, as a result, there is
also no support for collision between particles.

---

## Installation

This package is registered, simply use `Pkg.add("DynamicalBilliards")` to install it.
The [stable documentation](https://juliadynamics.github.io/DynamicalBilliards.jl/stable/) accompanies the version installed with `Pkg.add()`.

!!! note "Plotting"
    Plotting in `DynamicalBilliards` is done through `PyPlot` and it is available on-demand only. Simply use the function `DynamicalBilliards.enableplotting()` and it will define and bring into scope all the relevant names. Notice that you must be able to `using PyPlot` for plotting to work. If you are not sure about how to install PyPlot,
    simply run the commands:
    `ENV["PYTHON"]=""; Pkg.add("PyCall"); Pkg.add("PyPlot"); using PyPlot;`

After the first installation, it is advised to run the tests to be sure that
everything works as expected.
```julia
using DynamicalBilliards
DynamicalBilliards.test_options(print_info = true, long_tests = true)
Pkg.test("DynamicalBilliards")
```
The long tests take on average 4 minutes to complete, while the shorts need some
seconds.

---

## Usage
It is highly suggested to first visit the [Basic Usage](/basic/basic_usage) section,
for a general overview of how to use `DynamicalBilliards.jl`. The very important discussion about [numerical precision](/physics/#numerical-precision) is done in the [Physics](/physics) page.

The following tutorials offer detailed descriptions for various aspects of `DynamicalBilliards.jl`:

- [How to define your custom Billiard Table](/tutorials/billiard_table)
- [Using Ray-Splitting billiards](/tutorials/ray-splitting)
- [Visualizing the billiard table and animating the particle's evolution](/tutorials/visualizing)
- [Creating your own Obstacle Type](/tutorials/own_obstacle)
- [Examples page](/tutorials/examples)

The [Library](/basic/library) section has the docstrings of all exported names in convenient groups and the [Physics](/physics) briefly discusses physical aspects of billiard systems as well as inner workings of the package.

---
## Contributing
Everyone is welcomed to contribute to `DynamicalBilliards.jl`! If you have some new
algorithm, types of Obstacles or anything new to add, do not hesitate!

If you would like to help but do not have anything new to contribute, please go ahead
and take a look at the [GitHub issues page](https://github.com/JuliaDynamics/DynamicalBilliards.jl/issues) of the package.
Most of the issues are easy to solve and are there specifically for people that would
like to contribute.

---

## Julia Billiard Animation
The animation of a particle inside a "Julia" billiard was generated with the code:
```julia
using DynamicalBilliards
DynamicalBilliards.enableplotting()

bt = billiard_julia(plotit = true)
p = randominside(bt)

darkblue = (64/255, 99/255, 216/255)
lightblue = (102/255, 130/255, 223/255)
okwargs = Dict(:linewidth => 2.0, :color => lightblue)
pkwargs = Dict(:color => darkblue, :s => 150.0)

sname = "C:\\****\\anim"

animate_evolution(p, bt, 200; col_to_plot = 7,
particle_kwargs = pkwargs, orbit_kwargs = okwargs,
savefigs = true, savename = sname)

# use gifmaker.me to merge all figures into one .gif
# in a future update, automatic support will be added!
```

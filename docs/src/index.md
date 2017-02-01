![DynamicalBilliards Logo: The Julia billiard](http://i.imgur.com/NKgzYrt.gif)

`DynamicalBilliards.jl` is an easy-to-use, modular and extendable Julia package for
Dynamical Billiard systems in two dimensions.

A dynamical billiard is a system where a particle is propagating as a straight line from obstacle-to-obstacle, performing specular reflection at the boundary of the obstacles. Billiard systems have been used extensively in chaos and nonlinear dynamics and played an important role in the development of nonlinear science.

The "standard" billiard described above can be extended in many ways. The [wikipedia page](https://en.wikipedia.org/wiki/Dynamical_billiards) has many examples of different types of billiards. The types that are currently offered by this package, besides the standard one, are magnetic and ray-splitting billiards. In a magnetic billiard the particle's orbit is a circle (like electrons in a perpendicular magnetic field). In ray-splitting (aka semiclassical) billiards the particle may propagate *through* an obstacle, given some arbitrary transmission and refraction law.

This package does not support finite-sized particles and, as a result, there is
also no support for collision between particles.

## Installation

---

This package is registered, simply use `Pkg.add("DynamicalBilliards")` to install it.
The [stable documentation](https://datseris.github.io/DynamicalBilliards.jl/stable/) accompanies the version installed with `Pkg.add()`.

!!! warning "PyPlot dependency"
    This package has a dependency on `PyPlot` for its plotting features, because of its maturity, detailed documentation
    and vast library of features. If you are not sure about how to install PyPlot,
    simply run the commands:
    `ENV["PYTHON"]=""; Pkg.add("PyCall"); Pkg.add("PyPlot"); using PyPlot;`

If you want to use the
latest features, compatible with the [latest documentation](https://datseris.github.io/DynamicalBilliards.jl/latest/), use `Pkg.checkout("DynamicalBilliards")`.

After the first installation, it is advised to run the short tests to be sure that
everything works as expected. This will only take about 2 minutes:
```julia
using DynamicalBilliards
test_options(print_info = true)
Pkg.test("DynamicalBilliards")
```
If you do not want to see what tests are done, do not use any keywords.
If you use this package for scientific research, you should run the long tests at least once.
To do this, pass the keyword argument `long_tests = true` to the `test_options` function.
These tests take on average 20 minutes to complete.

## Usage

---

For a crash course on how to use `DynamicalBilliards.jl`, you should visit the [Basic Usage](/basic/basic_usage) section.

If however, you want to make the most out if, the following tutorials offer detailed descriptions:
- [How to define your custom Billiard Table](/tutorials/billiard_table)
- [Using Ray-Splitting billiards](/tutorials/ray-splitting)
- [Visualizing the billiard table and animating the particle's evolution](/tutorials/visualizing)
- [Creating your own Obstacle Type](/tutorials/own_obstacle)
- [Examples page](/tutorials/examples)

The [Library](/basic/library) section has the docstrings of all exported names in convenient groups.

## Example Animation
The animation of a particle inside a "Julia" billiard was generated with the code:
```
using DynamicalBilliards

jbt = billiard_julia()
p = randominside(bt)

darkblue = (64/255, 99/255, 216/255)
lightblue = (102/255, 130/255, 223/255)
okwargs = Dict(:linewidth => 2.0, :color => lightblue)
pkwargs = Dict(:color => darkblue, :s => 150.0)

sname = "C:\\****\\anim"

animate_evolution(p, bt, 200; col_to_plot = 7,
particle_kwargs = pkwargs, orbit_kwargs = okwargs,
savefigs = true, savename = sname)
```

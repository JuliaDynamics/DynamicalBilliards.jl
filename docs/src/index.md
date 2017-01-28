# Dynamical Billiards

`DynamicalBilliards.jl` is an easy-to-use, modular and extendable Julia package for Dynamical Billiard systems in two dimensions.


![Example of a dynamical billiard with magnetic field](http://i.imgur.com/OasQRyQ.gif)

A dynamical billiard is a system where a particle is propagating as a straight line from obstacle-to-obstacle, performing specular reflection at the boundary of the obstacles. Billiard systems have been used extensively in chaos and nonlinear dynamics and played an important role in the development of nonlinear science. 

The "standard" billiard described above can be extended in many ways. The [wikipedia page](https://en.wikipedia.org/wiki/Dynamical_billiards) has many examples of different types of billiards. The types that are currently offered by this package, besides the standard one, are magnetic and ray-splitting billiards. In a magnetic billiard the particle's orbit is a circle (like electrons in a perpendicular magnetic field). In ray-splitting (aka semiclassical) billiards the particle may propagate *through* an obstacle, given some arbitrary transmission and refraction law.

## Installation

---

This package is currently under the registration process. When this process is over, you can install the package using `Pkg.add("DynamicalBilliards")`.

After first installation, it is advised to run the short tests to be sure that everything works as expected. This will only take 2 minutes:
```julia
using DynamicalBilliards
test_options(print_info = false)
Pkg.test("DynamicalBilliards")
```
If you want to see what tests are done, use `print_info = true` (`false` is actually the default value). If you use this package for research purposes, you should run the long tests at least one. To do this, pass the keyword argument `long_tests = true` to the `test_options` function. These tests take on average 20 minutes to complete.

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

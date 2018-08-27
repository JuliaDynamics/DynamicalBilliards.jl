using DynamicalBilliards

using Documenter, PyPlot

# cd(Pkg.dir("DynamicalBilliards")*"/docs")
# mkpath("src/anim")
#
# download("https://github.com/JuliaDynamics/Tutorials-and-Resources/blob/master/billiard_animations/penta.mp4?raw=true",
#         "src/anim/penta.mp4")
# download("https://github.com/JuliaDynamics/Tutorials-and-Resources/blob/master/billiard_animations/inverse.mp4?raw=true",
#         "src/anim/inverse.mp4")
# download("https://github.com/JuliaDynamics/Tutorials-and-Resources/blob/master/billiard_animations/ray.mp4?raw=true",
#         "src/anim/ray.mp4")

using Pkg
pkg"add StaticArrays"

makedocs(modules=[DynamicalBilliards], doctest=false, root = @__DIR__)
close("all")

deploydocs(
    deps   = Deps.pip("Tornado>=4.0.0,<5.0.0", "mkdocs",
    "mkdocs-material" ,"python-markdown-math", "pygments", "pymdown-extensions"),
    repo   = "github.com/JuliaDynamics/DynamicalBilliards.jl.git",
    julia  = "nightly",
    osname = "linux"
)

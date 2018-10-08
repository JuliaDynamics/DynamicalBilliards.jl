using DynamicalBilliards

using Documenter, PyPlot

using Pkg
pkg"add StaticArrays"

ioff()
makedocs(modules=[DynamicalBilliards], doctest=false, root = @__DIR__)
close("all")

deploydocs(
    deps   = Deps.pip("mkdocs==0.17.5", "mkdocs-material==2.9.4",
    "python-markdown-math", "pygments", "pymdown-extensions"),
    repo   = "github.com/JuliaDynamics/DynamicalBilliards.jl.git",
    julia  = "1.0",
    osname = "linux"
)

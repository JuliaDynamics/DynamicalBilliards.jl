using DynamicalBilliards

using Documenter, PyPlot

makedocs(modules=[DynamicalBilliards], doctest=false, root = @__DIR__)
close("all")

deploydocs(
    deps   = Deps.pip("Tornado>=4.0.0,<5.0.0", "mkdocs",
    "mkdocs-material" ,"python-markdown-math", "pygments", "pymdown-extensions"),
    repo   = "github.com/JuliaDynamics/DynamicalBilliards.jl.git",
    julia  = "0.6",
    osname = "linux"
)

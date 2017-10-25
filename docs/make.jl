using DynamicalBilliards

using Documenter

makedocs(modules=[DynamicalBilliards], doctest=false)

deploydocs(
    deps   = Deps.pip("mkdocs==0.16.3",
    "mkdocs-material" ,"python-markdown-math", "pygments", "pymdown-extensions"),
    repo   = "github.com/JuliaDynamics/DynamicalBilliards.jl.git",
    julia  = "0.6",
    osname = "linux"
)

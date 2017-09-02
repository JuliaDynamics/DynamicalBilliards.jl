using Documenter, DynamicalBilliards

makedocs(modules=[DynamicalBilliards], doctest=false)

deploydocs(
    deps   = Deps.pip("mkdocs",
    "mkdocs-material" ,"python-markdown-math", "pygments", "pymdown-extensions"),
    repo   = "github.com/Datseris/DynamicalBilliards.jl.git",
    julia  = "0.6",
    osname = "linux"
)

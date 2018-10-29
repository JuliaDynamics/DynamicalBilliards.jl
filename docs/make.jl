using DynamicalBilliards

using Documenter, PyPlot, DocumenterMarkdown

ioff()
makedocs(modules=[DynamicalBilliards], doctest=false, root = @__DIR__, format = :markdown)
close("all")

if !Sys.iswindows()
    deploydocs(
        # deps   = Deps.pip("mkdocs==0.17.5", "mkdocs-material==2.9.4",
        deps   = Deps.pip("mkdocs", "mkdocs-material",
        "python-markdown-math", "pygments", "pymdown-extensions"),
        repo   = "github.com/JuliaDynamics/DynamicalBilliards.jl.git",
        target = "site",
        make = () -> run(`mkdocs build`)
    )
end

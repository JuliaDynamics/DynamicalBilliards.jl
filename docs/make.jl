using Documenter, DynamicalBilliards
 
makedocs(modules=[DynamicalBilliards], doctest=true)
 
deploydocs(deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/Datseris/DynamicalBilliards.git",
    julia  = "0.5.0",
    osname = "windows")

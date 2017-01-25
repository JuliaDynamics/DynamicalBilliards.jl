using Documenter, DynamicalBilliards
 
makedocs(modules=[DynamicalBilliards], doctest=false)

deploydocs(
    deps   = Deps.pip("mkdocs", "mkdocs-material" ,"python-markdown-math"),
    repo   = "github.com/USER_NAME/PACKAGE_NAME.jl.git",
    julia  = "0.5",
    osname = "windows"
)

cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
CI && Pkg.instantiate()

using DynamicalBilliards
using Documenter, PyPlot, DocumenterTools

# %%
using DocumenterTools: Themes
# download the themes
for file in ("juliadynamics-lightdefs.scss", "juliadynamics-darkdefs.scss", "juliadynamics-style.scss")
    download("https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/$file", joinpath(@__DIR__, file))
end
# create the themes
for w in ("light", "dark")
    header = read(joinpath(@__DIR__, "juliadynamics-style.scss"), String)
    theme = read(joinpath(@__DIR__, "juliadynamics-$(w)defs.scss"), String)
    write(joinpath(@__DIR__, "juliadynamics-$(w).scss"), header*"\n"*theme)
end
# compile the themes
Themes.compile(joinpath(@__DIR__, "juliadynamics-light.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-light.css"))
Themes.compile(joinpath(@__DIR__, "juliadynamics-dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))


ioff()


makedocs(
modules=[DynamicalBilliards],
doctest=false,
depth = 1,
sitename= "DynamicalBilliards.jl",
root = @__DIR__,
format = Documenter.HTML(
    prettyurls = CI,
    assets = [
        "assets/logo.ico",
        asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css),
        ],
    ),
pages = [
    "Introduction" => "index.md",
    "High Level API" => "basic/high_level.md",
    "Visualizing & Animating" => "visualizing.md",
    "Phase Spaces" => "basic/phasespaces.md",
    "Ray-Splitting" => "ray-splitting.md",
    "Lyapunov Exponents" => "lyapunovs.md",
    "MushroomTools" => "mushroomtools.md",
    "Physics" => "physics.md",
    "Internals" => "basic/low_level.md",
    "Tutorials" => [
        "Defining a Billiard" => "tutorials/billiard_table.md",
        "Defining your own Obstacles" => "tutorials/own_obstacle.md",
        "Examples" => "tutorials/examples.md",
    ]
],
)

close("all")
ion()

if CI
    deploydocs(
        repo = "github.com/JuliaDynamics/DynamicalBilliards.jl.git",
        target = "build",
        push_preview = true
    )
end

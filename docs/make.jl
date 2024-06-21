cd(@__DIR__)

using DynamicalBilliards
using CairoMakie

timeseries! = DynamicalBilliards.timeseries!

import Downloads
Downloads.download(
    "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/build_docs_with_style.jl",
    joinpath(@__DIR__, "build_docs_with_style.jl")
)
include("build_docs_with_style.jl")

using Literate
infile = joinpath(@__DIR__, "src", "billiards_visualizations.jl")
outdir = joinpath(@__DIR__, "src")
Literate.markdown(infile, outdir; credit = false)

pages = [
    "Introduction" => "index.md",
    "High Level API" => "basic/high_level.md",
    "Visualizing & Animating" => "billiards_visualizations.md",
    "Phase Spaces" => "basic/phasespaces.md",
    "Ray-Splitting" => "ray-splitting.md",
    "Lyapunov Exponents" => "lyapunovs.md",
    "MushroomTools" => "mushroomtools.md",
    "Physics" => "physics.md",
    "Internals" => "basic/low_level.md",
    "Tutorials" => [
        "Defining a Billiard" => "tutorials/billiard_table.md",
        "Defining your own Obstacles" => "tutorials/own_obstacle.md",
    ]
]

build_docs_with_style(pages, DynamicalBilliards;
    authors = "George Datseris <datseris.george@gmail.com>",
    expandfirst = ["index.md"],
    # warnonly = [:doctest, :missing_docs, :cross_references],
)

![DynamicalBilliards v3.0 Logo: The Julia billiard](https://github.com/JuliaDynamics/JuliaDynamics/blob/master/videos/billiards/DynamicalBilliards_logo_animated.gif?raw=true)

A Julia package for dynamical billiard systems in two dimensions.
The goals of the package is to provide a flexible and intuitive framework for fast implementation of billiard systems of arbitrary construction.

[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaDynamics.github.io/DynamicalBilliards.jl/dev)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaDynamics.github.io/DynamicalBilliards.jl/stable)
[![CI](https://github.com/JuliaDynamics/DynamicalBilliards.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/DynamicalBilliards.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/JuliaDynamics/DynamicalBilliards.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaDynamics/DynamicalBilliards.jl)
[![citation](http://joss.theoj.org/papers/753469f6b18c9c38127a7727d13c87cd/status.svg)](http://joss.theoj.org/papers/753469f6b18c9c38127a7727d13c87cd)

If you have used this package for research that resulted in a publication, please be kind enough to cite the papers listed in the [CITATION.bib](CITATION.bib) file.

## Features

Please see the [documentation](https://JuliaDynamics.github.io/DynamicalBilliards.jl/dev) for list of features, tutorials and installation instructions.

## Acknowledgements

This package is mainly developed by George Datseris. However, this development would not have been possible without significant help from other people:

1. [Lukas Hupe](https://github.com/lhupe)(@lhupe) Contributed the lyapunov spectrum calculation for magnetic propagation, implemented the boundary map function and did other contributions in bringing this package to version 2.0 (see [here](https://github.com/JuliaDynamics/DynamicalBilliards.jl/projects/1)).
1. [Diego Tapias](https://github.com/dapias) (@dapias) Contributed the lyapunov spectrum calculation method for straight propagation.
1. [David. P. Sanders](https://github.com/dpsanders) (@dpsanders) and [Ragnar Fleischmann](https://www.ds.mpg.de/person/20199/118124) contributed in fruitful discussions about the programming and physics of Billiard systems all-around.
2. [Christopher Rackauckas](https://github.com/ChrisRackauckas) (@ChrisRackauckas) helped set-up the continuous integration, testing, documentation publishing and all around package development-related concepts.
3. [Tony Kelman](https://github.com/tkelman) (@tkelman) helped significantly in the package publication process, especially in making it work correctly without destroying METADATA.jl.

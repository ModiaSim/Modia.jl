# TinyModia
 
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://modiasim.github.io/TinyModia.jl/stable/index.html)
[![The MIT License](https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square)](https://github.com/ModiaSim/TinyModia.jl/blob/master/LICENSE.md)

TinyModia is part of [ModiaSim](https://modiasim.github.io/docs/). 

TinyModia is a minimalistic environment in form of a Julia package to model and simulate physical systems (electrical, mechanical, thermo-dynamical, etc.) described by differential and algebraic equations. A user defines a model on a high level with model components (like a mechanical body, an electrical resistance, or a pipe) that are physically connected together. A model component is constructed by "expression = expression" equations. The defined model is symbolically processed (for example, equations might be analytically differentiated) with algorithms from package [ModiaBase.jl](https://github.com/ModiaSim/ModiaBase.jl). From the transformed model a Julia function is generated that is used to simulate the model with integrators from [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). 
The basic type of the floating point variables in the generated function is usually `Float64`, 
but can be set to any type `T<:AbstractFloat`, for example 
`Float32, DoubleFloat, Measurement{Float64}, StaticParticles{Float64,100}`.

## Installation
 
The package is registered and is installed with (Julia >= 1.5 is required):

```julia
julia> ]add TinyModia
```

It is recommended to also add the following packages, in order that all tests and examples can be executed in your standard environment:

```julia
julia> ]add ModiaPlot, Unitful, DifferentialEquations, Measurements, MonteCarloMeasurements, Distributions
```

## Main Developers

- [Hilding Elmqvist](mailto:Hilding.Elmqvist@Mogram.net), [Mogram](http://www.mogram.net/).

- [Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/),
  [DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en).

License: MIT (expat)

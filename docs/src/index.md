# TinyModia Documentation

[TinyModia](https://github.com/ModiaSim/TinyModia.jl) is a minimalistic environment in form of a Julia package to model and simulate physical systems (electrical, mechanical, thermo-dynamical, etc.) described by differential and algebraic equations. A user defines a model on a high level with model components (like a mechanical body, an electrical resistance, or a pipe) that are physically connected together. A model component is constructed by "expression = expression" equations. The defined model is symbolically processed (for example, equations might be analytically differentiated) with algorithms from package [ModiaBase.jl](https://github.com/ModiaSim/ModiaBase.jl). From the transformed model a Julia function is generated that is used to simulate the model with integrators from [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). 

The basic type of the floating point variables is usually `Float64`, but can be set to any
type `T<:AbstractFloat`, for example `Float32, DoubleFloats.DoubleFloat, Measurements.Measurement{Float64}, MonteCarloMeasurements.StaticParticles{Float64,100}`
via `@instantiateModel(..., FloatType = xxx)`.

  
## Installation

The package is currently being registered. During this phase, it is installed as
(Julia >= 1.5 is required):

```julia
julia> ]add https://github.com/ModiaSim/TinyModia.jl#main
```

It is recommended to also add the following packages, in order that all tests and examples can be executed:

```julia
julia> ]add Unitful, Measurements, MonteCarloMeasurements, Distributions
```

## Release Notes

### Version 0.7.0

- Initial version, based on code developed for Modia 0.6 and ModiaMath 0.6.


## Main developers

- Hilding Elmqvist (Mogram AB)
- Martin Otter ([DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)


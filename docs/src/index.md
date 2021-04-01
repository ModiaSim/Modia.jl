# TinyModia Documentation

[TinyModia](https://github.com/ModiaSim/TinyModia.jl) is a minimalistic environment in form of a Julia package to model and simulate physical systems (electrical, mechanical, thermo-dynamical, etc.) described by differential and algebraic equations. A user defines a model on a high level with model components (like a mechanical body, an electrical resistance, or a pipe) that are physically connected together. A model component is constructed by "expression = expression" equations. The defined model is symbolically processed (for example, equations might be analytically differentiated) with algorithms from package [ModiaBase.jl](https://github.com/ModiaSim/ModiaBase.jl). From the transformed model a Julia function is generated that is used to simulate the model with integrators from [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).
The basic type of the floating point variables is usually `Float64`, but can be set to any
type `FloatType<:AbstractFloat` via `@instantiateModel(..., FloatType = xxx)`, for example
it can be set to `Float32, DoubleFloat, Measurement{Float64}, StaticParticles{Float64,100}`.

## Installation

The package is registered and is installed with (Julia >= 1.5 is required):

```julia
julia> ]add TinyModia
```

It is recommended to also add the following packages, in order that all tests and examples can be executed in your standard environment:

```julia
julia> ]add ModiaPlot, Unitful, DifferentialEquations, Measurements, MonteCarloMeasurements, Distributions
```

## Release Notes

### Version 0.7.2

- Missing dependency of Test package added.

### Version 0.7.1

- Variable constructor `Var(...)` introduced. For example:
  `v = input | Var(init = 1.2u"m")`. 
  For details see section [A.1 Var constructor](@ref).

- Functions are called in the scope where macro [`@instantiateModel`](@ref) is called.

- New arguments of function [`simulate!`](@ref):
  - Parameter and init/start values can be changed with argument `merge`.
  - A simulation can be checked with argument `requiredFinalStates`.
  - Argument `logParameters` lists the parameter and init/start values used for the simulation.
  - Argument `logStates` lists the states, init, and nominal values used for the simulation.

- `end` in array ranges is supported, for example `v[2:end]`.

- New (small) model library `TinyModia/models/HeatTransfer.jl`.

- [TinyModia Tutorial](@ref) improved.

- [Functions](@ref) docu improved.

### Version 0.7.0

- Initial version, based on code developed for Modia 0.6 and ModiaMath 0.6.


## Main developers

- [Hilding Elmqvist](mailto:Hilding.Elmqvist@Mogram.net), [Mogram](http://www.mogram.net/).

- [Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/),
  [DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en).
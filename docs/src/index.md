# Modia Documentation

[Modia](https://github.com/ModiaSim/Modia.jl) is a minimalistic environment in form of a Julia package to model and simulate physical systems (electrical, mechanical, thermo-dynamical, etc.) described by differential and algebraic equations. A user defines a model on a high level with model components (like a mechanical body, an electrical resistance, or a pipe) that are physically connected together. A model component is constructed by "expression = expression" equations. The defined model is symbolically processed (for example, equations might be analytically differentiated) with algorithms from package [ModiaBase.jl](https://github.com/ModiaSim/ModiaBase.jl). From the transformed model a Julia function is generated that is used to simulate the model with integrators from [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).
The basic type of the floating point variables is usually `Float64`, but can be set to any
type `FloatType<:AbstractFloat` via `@instantiateModel(..., FloatType = xxx)`, for example
it can be set to `Float32, DoubleFloat, Measurement{Float64}, StaticParticles{Float64,100}`.

## Installation

The package is registered and is installed with (Julia >= 1.5 is required):

```julia
julia> ]add Modia
```

Furthermore, one or more of the following packages should be installed in order 
to be able to generate plots:

```julia
julia> ]add ModiaPlot_PyPlot        # if plotting with PyPlot desired
        add ModiaPlot_GLMakie       # if plotting with GLMakie desired
        add ModiaPlot_WGLMakie      # if plotting with WGLMakie desired
        add ModiaPlot_CairoMakie    # if plotting with CairoMakie desired
```

It is recommended to also add the following packages, in order that all tests and examples can be executed in your standard environment:

```julia
julia> ]add Unitful, DifferentialEquations, Measurements
        add MonteCarloMeasurements, Distributions
```

## Release Notes


### Version 0.8.0

- Improved scalability by using OrderedDicts instead of named tuples for models, variables and parameter modifications.
- Speed improvements for structural and symbolic algorithms.
- Added support for state events, time events and synchronous operators.
- Added support for mixed linear equation systems having Real and Boolean unknowns.
- Added support for user-defined components defined by structs and functions
  (multibody modeling with Modia3D is based on this feature).
  This makes it possible to utilize algorithms specialized for a component.
- Added support for numerical and analytic linearization.
- Added support for propagation of parameters (e.g. deep in a model, the value of a parameter can be defined as a function of some top level parameter and this parameter is changed before simulation starts).
- New small model libraries Translational.jl and PathPlanning.jl added.
- Result storage changed: `sol = simulate!(...)` calls internally `sol = solve(..)` from   DifferentialEquations.jl. `sol` contains time and the states at the communication time grid and
  at events. This is now kept in simulate(..), so the return value of simulate!(..) can be exactly used as if `solve(..)` would have been used directly.
- The plot(..) command now supports the following underlying plot packages: 
  [PyPlot](https://github.com/JuliaPy/PyPlot.jl),
  [GLMakie](https://github.com/JuliaPlots/GLMakie.jl),
  [WGLMakie](https://github.com/JuliaPlots/WGLMakie.jl), and
  [CairoMakie](https://github.com/JuliaPlots/CairoMakie.jl).
  It is also possible to select `NoPlot`, to ignore `plot(..)` calls 
  or `SilenNoPlot` to ignore `plot(..)` calls silently. The latter is useful for `runtests.jl`.
  Note, often [PyPlot](https://github.com/JuliaPy/PyPlot.jl) is the best choice.

Changes that are **not backwards compatible** to version 0.7.x:

- Models are OrderedDicts and no longer NamedTuples.

- simulate!(..): 
  - If FloatType=Float64 and no algorithm is defined, then Sundials.CVODE\_BDF() is used
    instead of the default algorithm of DifferentialEquations as in 0.7. The reason is that Modia models
    are usually large and expensive to evaluate and have often stiff parts, so that multi-step
    methods are often by far the best choice. CVODE_BDF() seems to be a good choice in many applications
    (another algorithm should be used, if there are many events, highly oscillatory vibrations, or if all states are non-stiff). 
  - The default value of `stopTime` is equal to `startTime` (which has a default value of 0.0 s), and is no longer 1.0 s.

- Plotting is defined slightly differently (`@useModiaPlot`, instead of `using ModiaPlot`).


### Version 0.7.3

- Evaluation and propagation of parameter expressions (also in simulate!(..., merge=Map(...))).
  Propagation of start/init values of states is not yet supported.

- State events supported.


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

- New (small) model library `Modia/models/HeatTransfer.jl`.

- Modia Tutorial improved.

- [Functions](@ref) docu improved.

### Version 0.7.0

- Initial version, based on code developed for Modia 0.6 and ModiaMath 0.6.


## Main developers

- [Hilding Elmqvist](mailto:Hilding.Elmqvist@Mogram.net), [Mogram](http://www.mogram.net/).

- [Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/),
  [DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en).
# Modia Documentation

[Modia](https://github.com/ModiaSim/Modia.jl) is an environment in form of a Julia package to model and simulate physical systems (electrical, mechanical, thermo-dynamical, etc.) described by differential and algebraic equations. A user defines a model on a high level with model components (like a mechanical body, an electrical resistance, or a pipe) that are physically connected together. A model component is constructed by **`expression = expression` equations** or by Julia structs/functions, such as the pre-defined Modia 3D-mechanical components. The defined model is symbolically processed (for example, equations might be analytically differentiated) with algorithms from package [ModiaBase.jl](https://github.com/ModiaSim/ModiaBase.jl). From the transformed model a Julia function is generated that is used to simulate the model with integrators from [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).
The basic type of the floating point variables is usually `Float64`, but can be set to any
type `FloatType<:AbstractFloat` via `@instantiateModel(..., FloatType = xxx)`, for example
it can be set to `Float32, DoubleFloat, Measurement{Float64}, StaticParticles{Float64,100}`.

Modia includes a multibody program and 3D shapes for visualization and collision handling. It is then, for example, possible to model the 3D mechanical part of a robot with Modia multibody components and the electrical motors and gearboxes that are driving the joints with equation-based Modia components. Collision handling with elastic response calculation is performed for shapes that are defined with a contact material and have a convex geometry or are approximated by the convex hull of a concave geometry.
The multibody program supports currently tree-structured multibody systems, but does not (yet) support kinematic loops.

Example videos:

- [YouBot robots with gripping](https://modiasim.github.io/Modia3D.jl/resources/videos/YouBotsGripping.mp4)
- [Billiard table with 16 balls](https://modiasim.github.io/Modia3D.jl/resources/videos/Billard16Balls.mp4)


## Installation

The package is registered and is installed with (Julia >= 1.5 is required):

```julia
julia> ]add ModiaBase, ModiaLang, Modia
```

Modia exports all exported symbols of [ModiaLang](https://github.com/ModiaSim/ModiaLang.jl), [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl) and of [Unitful](https://github.com/PainterQubits/Unitful.jl)
and the essential exported symbols fo [Modia3D](https://github.com/ModiaSim/Modia3D.jl).

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
julia> ]add Measurements, MonteCarloMeasurements, Distributions
```

## Release Notes


### Version 0.5.0

- Using and reexporting ModiaLang 0.8.1 (see [release notes](https://modiasim.github.io/ModiaLang.jl/stable/)).
- Using and partially reexporting Modia3D 0.5.0 (see [release notes](file:///D:/otter/home/.julia/dev/Modia3D/docs/build/index.html#Release-Notes)).
- New plot package interface via [ModiaResult](https://github.com/ModiaSim/ModiaResult.jl). Additional support for PyPlot, WGLMakie, CairoMakie (besides GLMakie).


### Version 0.4.0

- Initial version of new Modia design.


## Main developers

### ModiaLang, ModiaBase

- [Hilding Elmqvist](mailto:Hilding.Elmqvist@Mogram.net), [Mogram](http://www.mogram.net/).
- [Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/),
  [DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en).
  
### Modia3D

[Andrea Neumayr](mailto:andrea.neumayr@dlr.de),
[Martin Otter](https://rmc.dlr.de/sr/de/staff/martin.otter/) and
[Gerhard Hippmann](mailto:gerhard.hippmann@dlr.de),\
[DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en)
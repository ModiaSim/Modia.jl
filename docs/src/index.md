# Modia Documentation

[Modia](https://github.com/ModiaSim/Modia.jl) is a Julia package for modeling and simulation of multidomain engineering systems described by differential equations, algebraic equations, and (space-discretized) partial differential equations. It shares many powerful features of the [Modelica language](https://www.modelica.org/modelicalanguage), together with features not available in Modelica. 

A user defines a model on a high level with model components (like an electrical resistance, a rotational inertia, a rod with heat transfer, a PI controller etc.) that are physically connected together. A model component is constructed by "expression = expression" equations. The defined model is symbolically transformed to ODE form dx/dt = f(x,t). For example, equations might be analytically differentiated, ODE states selected, linear equation systems numerically solved when evaluating the transformed model. From the transformed model a Julia function is generated that is used to simulate the model with integrators from [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl). Simulation results can be plotted with [ModiaPlot](https://github.com/ModiaSim/ModiaPlot.jl), that provides a convenient interface to [GLMakie](https://github.com/JuliaPlots/GLMakie.jl) line plots.

Other packages from the Julia eco-systems that are specially supported:

- [Unitful](https://github.com/PainterQubits/Unitful.jl) to define and process physical units.
- [Measurements](https://github.com/JuliaPhysics/Measurements.jl) to perform simulations with uncertain parameters and initial values with linear propagation theory.
- [MonteCarloMeasurements](https://github.com/baggepinnen/MonteCarloMeasurements.jl) to perform simulations with uncertain parameters and initial values with particle theory.


## Installation

Modia and ModiaPlot are registered and are installed with (Julia >= 1.5 is required):

```julia
julia> ]add Modia, ModiaPlot
```

It is recommended to also add the following packages, in order that all tests and examples can be executed in your standard environment:

```julia
julia> ]add Measurements 
        add MonteCarloMeasurements
        add Distributions
        add Interpolations
```

## Release Notes

### Version 0.7.2

- Initial version, based on code developed for Modia 0.6, ModiaMath 0.6,
  and TinyModia 0.7.2.


## Main developers

- [Hilding Elmqvist](mailto:Hilding.Elmqvist@Mogram.net), [Mogram](http://www.mogram.net/).

- [Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/),
  [DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en).
# Modia.jl

[![Travis](https://travis-ci.org/ModiaSim/Modia.jl.svg?branch=master)](https://travis-ci.org/ModiaSim/Modia.jl)
[![AppVoyer](https://ci.appveyor.com/api/projects/status/github/ModiaSim/Modia.jl?svg=true)](https://ci.appveyor.com/project/MartinOtter/modia-jl)
[![codecov.io](http://codecov.io/github/ModiaSim/Modia.jl/coverage.svg?branch=master)](http://codecov.io/github/ModiaSim/Modia.jl?branch=master)
[![The MIT License](https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square)](https://github.com/ModiaSim/Modia.jl/blob/master/LICENSE)

## Introduction

[Modia](https://modiasim.github.io/) is a domain specific extension of [Julia](http://julialang.org/ "Julia") for **modeling and simulation of physical systems**.

Papers and presentations about Modia:

- [Overview about Modia](http://link.springer.com/chapter/10.1007%2F978-3-319-47169-3_15) ([ISoLA conference Oct. 2016](http://www.isola-conference.org/isola2016/); slides in [pdf](https://modiasim.github.io/Modia.jl/slides/Systems-Modeling-and-Programming-Slides.pdf) format).
- [Overview and new features in Modia](https://www.modelica.org/events/modelica2017/proceedings/html/submissions/ecp17132693_ElmqvistHenningssonOtter.pdf) ([12th International Modelica Conference, May 2017](https://www.modelica.org/events/modelica2017/proceedings/html/index.html); slides in [pdf](https://modiasim.github.io/Modia.jl/slides/Innovations-for-Future-Modelica.pdf) format).
- [New algorithms in Modia](https://www.modelica.org/events/modelica2017/proceedings/html/submissions/ecp17132565_OtterElmqvist.pdf) ([12th International Modelica Conference, May 2017](https://www.modelica.org/events/modelica2017/proceedings/html/index.html); slides in [pptx](https://modiasim.github.io/Modia.jl/slides/Modelica2017-DAAE-Transformation.pptx) and [pdf](https://modiasim.github.io/Modia.jl/slides/Modelica2017-DAAE-Transformation.pdf) format).
- [Modia: A Domain Specific Extension of Julia for Modeling and Simulation](http://juliacon.org/2017/talks#talk-21) ([juliaCon 2017, June 2017](http://juliacon.org/2017/); recording at [YouTube](https://youtu.be/hVg1eL1Qkws) and slides in [pdf](https://modiasim.github.io/Modia.jl/slides/Modia-JuliaCon-2017.pdf) format).

Modia is designed to model and simulate physical systems (electrical, mechanical, thermo-dynamical, etc.) described by differential and algebraic equations. A user defines a model on a high level with model components (like a mechanical body, an electrical resistance, or a pipe) that are physically connected together. A model component is constructed by "expression = expression" equations. The defined model is symbolically processed (for example, equations might be analytically differentiated), JIT compiled and simulated with [Sundials IDA solver](http://computation.llnl.gov/projects/sundials/ida) with the KLU sparse matrix package. By this approach it's possible and convenient to build models with hundred thousands of equations describing the dynamics of a car, an airplane, a power plant, etc. and simulate them. The authors used previous experience from the design of the modeling language [Modelica](https://www.modelica.org/). Modia will also be used to design and evaluate features for future Modelica versions.

Component models are defined by `@model` macros. Such models contain definition of variables with various attributes such as start values, min, max, SI unit, etc. An `@equations` macro is used to define the equations of such a component. Coupling between components is expressed using a connect statement involving groups of variables. The semantics is either to constrain connected variables to be equal or to constrain a set of variables to sum to zero, for example to model Kirchhoff's current law.

## Installation

**Modia** is registered in METADATA.jl and can be installed with Pkg.add.
The latest released version (0.2.3) is the last one with support for Julia >= 0.6. Trunk and later versions support Julia >=1.0.

```julia
# Julia >= 0.6:
julia> Pkg.add("Modia")
julia> Pkg.add("ModiaMath")

# Julia >= 0.7:
julia> ]add Modia ModiaMath
```

Modia uses [PyPlot](https://github.com/JuliaPy/PyPlot.jl) for plotting.
If `PyPlot` is not available in your current Julia environment
an information message is printed and all `plot(..)` calls are ignored.

In order that plot windows are displayed, you need to add `PyPlot` to your current environment
via `]add PyPlot`. Often this automatic installation fails and it is recommended to follow
the instructions
[Installing PyPlot in a robust way](https://github.com/ModiaSim/ModiaMath.jl/wiki/Installing-PyPlot-in-a-robust-way).


## Use

### To define a model:

```julia
  using Modia
  @model FirstOrder begin
     x = Variable(start=1)   # start means x(0)
     T = Parameter(0.5)      # Time constant
     u = 2.0                 # Same as Parameter(2.0)
  @equations begin
     T*der(x) + x = u        # der() means time derivative
     end
  end;
```

### To simulate a model:
```julia
  result = simulate(FirstOrder, 2);
  @show result["x"][end];
  ModiaMath.plot(result, "x")
```

## Examples

The schematics below are screenshots of [Modelica models](https://www.modelica.org/). These models have been converted to Modia and the examples below execute these models. Note, in Modia there is not (yet) a graphical definition of models.

### Current Controller
![Current Controller](https://github.com/ModiaSim/Modia.jl/blob/master/docs/CurrentController.png "Multi-domain model: Current Controller")

[`examples/CurrentController.jl`](examples/CurrentController.jl)

### Cauer Low Pass Filter
![Cauer Low Pass Filter](https://github.com/ModiaSim/Modia.jl/blob/master/docs/CauerLowPassFilter.png "Electrical model: Cauer Low Pass Filter")

[`examples/CauerLowPassFilter.jl`](examples/CauerLowPassFilter.jl)


### To run examples:
```julia
  using Modia
  include("$(Modia.ModiaDir)/examples/runexamples.jl")
```

### To run tests:
```julia
  using Modia
  include("$(Modia.ModiaDir)/test/runtests.jl")
```

---
### Examples (Jupyter Notebooks)
- [Electrical low pass filter](https://github.com/ModiaSim/Modia.jl/blob/master/docs/LPfilter.ipynb).
- [Rectifier - State Events](https://github.com/ModiaSim/Modia.jl/blob/master/docs/Rectifier.ipynb).

## Status

The version released now is partial since certain prototype functionalities needs to be generalized and refactored. See below. Such prototype features are enabled by flags in the simulate command. See examples.

## Available functionalities

- Model instantiation (including handling of modifiers and extends)
- Model flattening
- Redeclarations
- Handling of T, size, start, state and flow attributes
- Size deduction of variables and equations
- Array equations
- Index reduction
- BLT
- Symbolic transformations of equations (solving, differentiating)
- Algebraic loop handling
- Logging of transformation steps
- Simulation with ModiaMath

## To Do
### Enhance and refactor prototype codes for:

- Alias handling
- Handle overdetermined equations
- Introduction of partial and block attribute to models
- Automatic state selection
- Arrays of components
- Complex data type
- Event handling
- Synchronous equations
- Sparse Jacobian handling
- Impulse handling
- API to dynamically change model topology
- More models converted from Modelica Standard Library

### To implement

- Tearing of systems of equations
- Improved code generation for really large models
- Allowing change of parameters without recompilation
- Taking min and max attributes into account
- Handling of rotation matrices involved in algebraic loops
- Linearization
- ...

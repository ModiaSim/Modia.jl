# ModiaLang.jl
 
[![The MIT License](https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square)](https://github.com/ModiaSim/ModiaLang.jl/blob/master/LICENSE.md)

ModiaLang is part of [ModiaSim](https://modiasim.github.io/docs/).

ModiaLang is usually used via [Modia](https://github.com/ModiaSim/Modia.jl).

ModiaLang provides the core equation-based language of Modia, transformation of a model to ODE form dx/dt = f(x,t) and a thin interface to [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl).

## Installation
 
Typically, a user installs [Modia](https://github.com/ModiaSim/Modia.jl) and does not need
to install ModiaLang separately. 

If needed, ModiaLang is installed with (Julia >= 1.5 is required):

```julia
julia> ]add ModiaLang
```

## Main Developers

- [Hilding Elmqvist](mailto:Hilding.Elmqvist@Mogram.net), [Mogram](http://www.mogram.net/).

- [Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/),
  [DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en).

License: MIT (expat)

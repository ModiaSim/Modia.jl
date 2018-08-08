#=
Modia main module.

* Developers: Hilding Elmqvist, Mogram AB, Toivo Henningsson, Lund and Martin Otter, DLR 
* Copyright (c) 2016-2018: Hilding Elmqvist, Toivo Henningsson, Martin Otter
* License: MIT (expat)
=#


import ModiaMath

doc"""
Modia - Dynamic Modeling and Simulation in Julia

To define a model:

```julia
  @model FirstOrder begin
     x = Variable(start=1)   # start means x(0)
     T = Parameter(0.5)      # Time constant
     u = 2.0                 # Same as Parameter(2.0)
  @equations begin
     T*der(x) + x = u        # der() means time derivative
     end
  end;
```

To simulate a model:
```julia
  result = simulate(FirstOrder, 2);
  @show result["x"][end];
  ModiaMath.plot(result, "x")
```

To run examples:
```julia
  include("$(Modia.ModiaDir)/examples/runexamples.jl")
```

To run tests:
```julia
  include("$(Modia.ModiaDir)/test/runtests.jl")
```
For more information, see (https://github.com/ModiaSim/Modia.jl/blob/master/README.md)
"""
module Modia

const Version = "0.2.0-beta.2"
const Date = "2018-08-08"

#println(" \n\nWelcome to Modia - Dynamic MODeling and Simulation in julIA")
print(" \n\nWelcome to ")
print("Mod")
print_with_color(:red, "ia", bold=true)
print(" - ")
print_with_color(:light_black, "Dynamic ")
print("Mod")
print_with_color(:light_black, "eling and Simulation with Jul")
print_with_color(:red, "ia", bold=true)
println()

println("Version $Version ($Date)")
println("Type \"?Modia\" for help.\n\n")

const ModiaDir = Pkg.dir("Modia")
export ModiaDir

export @model, simulateModel, simulate, checkSimulation, simulateMultiModeModel
export Variable, Float, Boolean, Integ, Str, Parameter, Var, Par, undefined
export Variability, constant, parameter, discrete, continuous
export Property, general, symmetric, orthogonal, rotationGroup3D
export skew, skewCoords
export allInstances
export @component, addComponent!

include("language/ModiaLogging.jl")
include("language/Instantiation.jl")
include("language/Execution.jl")

include("symbolic/DAEquations/Synchronous.jl") # Before models/Electric, etc

include("models/models.jl")  # Before symbolic because MultiBody is used in BasicStructuralTransform

include("symbolic/symbolic.jl")

using .Instantiation

include("models/ModiaBase.jl")

end

#=
Modia main module.

* Developers: Hilding Elmqvist, Mogram AB, Toivo Henningsson, Lund and Martin Otter, DLR
* Copyright (c) 2016-2019: Hilding Elmqvist, Toivo Henningsson, Martin Otter
* License: MIT (expat)
=#


#import ModiaMath

"""
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

const Version = "0.3.0"
const Date = "2019-04-07"


#println(" \n\nWelcome to Modia - Dynamic MODeling and Simulation in julIA")
print(" \n\nWelcome to ")
print("Mod")
@static if VERSION < v"0.7.0-DEV.2005"
    print_with_color(:red, "ia", bold=true)
    print(" - ")
    print_with_color(:light_black, "Dynamic ")
    print("Mod")
    print_with_color(:light_black, "eling and Simulation with Jul")
    print_with_color(:red, "ia", bold=true)
else
    printstyled("ia", bold=true, color=:red)
    print(" - ")
    printstyled("Dynamic ", color=:light_black)
    print("Mod")
    printstyled("eling and Simulation with Jul", color=:light_black)
    printstyled("ia", bold=true, color=:red)
end

println()
println("Version $Version ($Date)")
println("Type \"?Modia\" for help.\n\n")

const ModiaDir = dirname(Base.@__DIR__)

export ModiaDir

export @model, simulateModel, simulate, checkSimulation, simulateMultiModeModel
export Variable, Float, Boolean, Integ, Str, Parameter, Var, Par, undefined
export Variability, constant, parameter, discrete, continuous
export Property, general, symmetric, orthogonal, rotationGroup3D
export skew, skewCoords
export allInstances
export @component, addComponent!

# Import packages that are used in examples and tests
# (in order that there are no requirements on the environment
#  in which the examples and tests are executed).
import DataStructures
import ModiaMath
import StaticArrays
import Unitful

@static if VERSION >= v"0.7.0-DEV.2005"
    import LinearAlgebra
    import Markdown
    import SparseArrays
    import Test
end


# Include all sub-modules of Modia
include("language/ModiaLogging.jl")
include("language/Instantiation.jl")
include("language/Execution.jl")

include("symbolic/DAEquations/Synchronous.jl") # Before models/Electric, etc
include("symbolic/symbolic.jl")

using .Instantiation

include("models/ModiaBase.jl")
include("models/models.jl")  # Before symbolic because MultiBody is used in BasicStructuralTransform



end

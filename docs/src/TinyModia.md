# Introduction to TinyModia

In this chapter an equation and object based language called **TinyModia** is defined that is used to evaluate and demonstrate the basic algorithms.

## 2.1 Equation oriented TinyModia models

To define models, a constructor `Model` taking a comma separated list of name/value pairs is used:

A simple differential equation

```math
T \cdot \frac{dx}{dt} + x = 2
```

can be defined as:

```julia
    using TinyModia

    SimpleModel = Model(
        T = 0.2,
        equation = :(T * der(x) + x = 2),
    )
```
The model consist of a definition of a parameter `T` and one equation. An equation can have a Julia expression on both sides of the equal sign. This model will be symbolically solved for the derivative `der(x)` before simulation, so the following equation will be used for the integration:

```math
\frac{dx}{dt} = (2 - x) / T
```

A low pass filter block with input `u` and output `y`

```math
\begin{aligned}
T \cdot \frac{dx}{dt} + x &= u\\
     y &= x \\
x(t_0) &= 0
\end{aligned}
```
can be defined as:

```julia
    using ModiaBase

    LowPassFilter = Model(
        T = 0.2,
        inputs = :[u],
        outputs = :[y],
        init = Map(x=0),
        equation = :[T * der(x) + x = u],
        y = :x,
    )
```
The `init` key defines the initial condition of the state `x` to 0. A constructor `Map` is used. If an equation has just a unique variable in the left hand side, this variable can be used as a key and the corresponding value is the left hand side of the equation given as a quoted expression `y = :x`.

## 2.2 Merging models

It is possible to combine models by merging. If we want to change the model to become a high passfilter, an alternative output equation

```math
y = -x + u
```

is defined in an anonymous Model `Model( y = :(-x + u) )`. The merging can in this case be made with the Julia `merge` function:

```julia
HighPassFilter = merge(LowPassFilter, Model( y = :(-x + u) ) )
```

In general, also recursive merging is desired and TinyModia provides a `mergeModels` function for that (see appendix 2). This function can also be invoked as a binary operator `|` (also used for merge in Python). Note, that the order of the arguments/operands are important.

Generalizing the block to have two outputs for both low and high pass filtering would be done as follows:

```julia
LowAndHighPassFilter = LowPassFilter | Model(
        outputs = :[low, high],
        y = nothing,
        low = :x,
        high = :(-x + u)
    )
```
The equation for `y` is removed by "assigning" `nothing` and two variables are defined and declared as outputs.
Model `LowAndHighPassFilter` represents the following equations:

```math
\begin{aligned}
T \cdot \frac{dx}{dt} + x &= u\\
     low &= x \\
    high &= -x + u \\
x(t_0) &= 0
\end{aligned}
```

By turning on logging of merging `setLogMerge(true)`, the translator gives the log:

```julia
Changing: outputs = [y] to outputs = [low, high]
Deleting: y
Adding: low = x
Adding: high = -x + u
```

The resulting model is pretty printed by calling `@showModel(LowAndHighPassFilter)`:

```julia
LowAndHighPassFilter = Model(
  T = 0.2,
  inputs = [u],
  outputs = [low, high],
  init = (
    x = 0,
  ),
  equations = [T * der(x) + x = u],
  low = x,
  high = -x + u,
),
```

In order to test such an input/output block, an input needs to be defined. This can be made by adding an equation for `u`. Assume we want `u` to be sinousoidial with an increasing frequency:

```julia
using Unitful

TestLowAndHighPassFilter = LowAndHighPassFilter | Model(
        w = :((time+1u"s")*u"1/s/s"),
        u = :(sin(w*time)*u"V"),
        init = Map(x=0.2u"V")
    )
```

`time` is a reserved name for the independent variable. It has unit `s` for seconds. The Julia package [Unitful](https://painterqubits.github.io/Unitful.jl/stable/) provides a means for defining units and managing unit inference. Definition of units is done with a string macro `u"..."`. In this case, the input signal was given unit Volt. The state x must then also have consistent unit, i.e. Volt. If the model equations contain systems of simultaneous equations, then approximate guess values, optionally with units, must be given `start`: `start = Map(i=0.0u"A")`.

## 2.3 Hierarchical modeling

Sofar, the composition of models have resulted in named tuples with values being numeric values or quoted expressions. Hierarchical models are obtained if the values themself are named tuples. A model with two filters can, for example, be defined as follows:

```julia
TwoFilters = (
    high = HighPassFilter,
    low = LowPassFilter,
)
```

Note, that the previous definitions of HighPassFilter and LowPassFilter was used instead of making the defintions inline.

A band pass filter is a series connection of a high pass filter and a low pass filter and is described below:

```julia
BandPassFilter = (
    inputs = :[u],
    outputs = :[y],
    high = HighPassFilter | Map(T=0.5),
    low = LowPassFilter | Map(T=2.0),
    equations = :[
        high.u = u,
        low.u = high.y,
        y = low.y]
)
```

A new input have been defined which is propagated to `high.u`. The series connection itself is obtained by the equation `low.u = high.y`. Note, that dot-notation is allowed in equations.

The input and output for the BandPassFilter when using the same input definition as for the TestLowPassFilter is shown below:

![Band Pass Filter Plot](../resources/images/BandPassFilterPlot.png)

The above examples are available in file SimpleFilters.jl.

## 2.4 Physically oriented modeling

Sofar, only signal flow modeling has been used, i.e. input/output blocks coupled with equations between outputs and inputs. For object oriented modeling more high level constructs are neccessary. Coupling is then acausal and involves potentials such as electric potential, positions, pressure, etc. and flows such as electric current, forces and torques and mass flow rate.

### 2.4.1 Connectors

Models which contain any flow variable, i.e. included in the list `flows`, are considered connectors. Connectors must have equal number of flow and potential variables, i.e. included in the list `potentials`, and have matching array sizes. Connectors may not have any equations. An example of an electrical connector with potential (in Volt) and current (in Ampere) is shown below.

```julia
Pin = Model( potentials = :[v], flows = :[i] )
```

### 2.4.2 Components

Components are declared in a similar ways as blocks. However, the interfaces between components are defined using connector instances.

An electrical resistor can be descibed as follows:

```julia
Resistor = Model(
    R = 1.0u"Ω",
    p = Pin,
    n = Pin,
    equations = :[
        0 = p.i + n.i
        v = p.v - n.v
        i = p.i
        R*i = v ]
    )
```

### 2.4.3 Inheritance

Various physical components sometimes share common properties. One mechanism to handle this is to use inheritance. In TinyModia, merging is used.

Electrical components such as resistors, capacitors and inductors are categorized as oneports which have two pins. Common properties are: constraint on currents at the pins and definitions of voltage over the component and current through the component.

```julia
OnePort = Model(
    p = Pin,
    n = Pin,
    equations = :[
        0 = p.i + n.i
        v = p.v - n.v
        i = p.i ] )
```

Having such a OnePort definition makes it convenient to define electrical component models by merging OnePort with specific parameter definitions with default values and equations:

```julia
Resistor = OnePort | Model( R = 1.0u"Ω", equation = :[ R*i = v ], )

Capacitor = OnePort | Model( C = 1.0u"F", init=Map(v=0.0u"V"), equation = :[ C*der(v) = i ] )

Inductor = OnePort | Model( L = 1.0u"H", init=Map(i=0.0u"A"), equation = :[ L*der(i) = v ] )

ConstantVoltage = OnePort | Model( V = 1.0u"V", equation = :[ v = V ] )
```
The merged `Resistor` is shown below:

```julia
Resistor = Model(
  p = (
    potentials = [v],
    flows = [i],
  ),
  n = (
    potentials = [v],
    flows = [i],
  ),
  equations = [v = p.v - n.v; 0 = p.i + n.i; i = p.i],
  R = 1.0 Ω,
  equation = [R * i = v],
),
```

### 2.4.4 Connections

Connections are described as an array of tuples listing the connectors that are connected:
```julia
    ( <connect reference 1>, <connect reference 2>, ... )
```
A connect reference has either the form 'connect instance name' or 'component instance name'.'connect instance name' with 'connect instance name' being either a connector instance, input or output variable.

Examples
```julia
    connect = :[
      (V.p, R1.p)
      (R1.n, p)
      (C1.n, V.n, R2.p)
      ...
    ]
```

For connectors, all the potentials of the connectors in the same connect tuple are set equal and the sum of all incoming flows to the model are set equal to the sum of the flows into sub-components.

### 2.4.5 Connected models

Having the above electrical component models, enables defining a filter

![Filter Circuit](../resources/images/Filter.png)

by instanciating components, setting parameters and defining connections.


```julia
Filter = (
    R = Resistor | Map(R=0.5u"Ω"),
    C = Capacitor | Map(C=2.0u"F"),
    V = ConstantVoltage | Map(V=10.0u"V"),
    connect = :[
      (V.p, R.p)
      (R.n, C.p)
      (C.n, V.n)
    ]
)
```

The connect tuples are translated to:

```julia
  V.p.v = R.p.v
  V.p.i + R.p.i = 0
  R.n.v = C.p.v
  R.n.i + C.p.i = 0
  C.n.v = V.n.v
  C.n.i + V.n.i = 0
```

### 2.4.6 Parameter propagation

Hierarchical modification of parameters is powerful but sometimes a bit inconvenient. It is also possible to propagate parameters intoduced on a high level down in the hierarchy. The following Filter model defines three parameters, `r`, `c` and `v`. The `r` parameter is used to set the resistance of the resistor R: `Map(R=:(up.r))`. A special notation `up.` is used to denote that `r` should not be searched within the resistor, but in the enclosing model, i.e. Filter.

```julia
Filter = Model(
    r = 1.0u"Ω",
    c = 1.0u"F",
    v = 1.0u"V",
    R = Resistor | Map(R=:(up.r)),
    C = Capacitor | Map(C=:(up.c)),
    V = ConstantVoltage | Map(V=:(up.v)),
    connect = :[
      (V.p, R.p)
      (R.n, C.p)
      (C.n, V.n)
    ]
)
```

Two separate filters can then be defined with:

```julia
TwoFilters = Model( f1 = Filter | Map( r = 10.0, c = 2.0), f2 = Filter )
```

### 2.4.7 Redeclarations

It is possible to reuse a particular model topology by redeclaring the model of particular components. For example, changing the filter `f1` to a voltage divider by changing C from a Capacitor to a Resistor. A predefined model `Redeclare` is used for this purpose.

```julia
VoltageDividerAndFilter = TwoFilters | Map(f1 = Map(C = Redeclare | Resistor | Map(R = 20.0)))
```

By using `Redeclare`, a new model based on a Resistor is used for `C` and the usual merge semantics with the previously defined model of `C` is not used.

The above examples are available in file FilterCircuit.jl.

# 3 Simulation

A particular model needs to be instantiated and simulated.

```julia
    using ModiaBase
    using ModiaPlot

    filter = @instantiateModel(Filter)
    simulate!(filter, stopTime=10.0)
    plot(filter, "y", figure=1)
```

The `@instantiateModel` macro takes additional arguments: 

```julia
    modelInstance = @instantiateModel(model; modelName="", modelModule=nothing, FloatType = Float64, aliasReduction=true, unitless=false,
        log=false, logModel=false, logDetails=false, logStateSelection=false, logCode=false, logExecution=false, logTiming=false)
    
Instantiates a model, i.e. performs structural and symbolic transformations and generates a function for calculation of derivatives suitable for simulation.

* `model`: model (declarations and equations)
* `FloatType`: Variable type for floating point numbers, for example: Float64, Measurements{Float64}, StaticParticles{Float64,100}, Particles{Float64,2000}
* `aliasReduction`: Perform alias elimination and remove singularities
* `unitless`: Remove units (useful while debugging models and needed for MonteCarloMeasurements)
* `log`: Log the different phases of translation
* `logModel`: Log the variables and equations of the model
* `logDetails`: Log internal data during the different phases of translation
* `logStateSelection`: Log details during state selection
* `logCode`: Log the generated code
* `logExecution`: Log the execution of the generated code (useful for finding unit bugs)
* `logTiming`: Log timing of different phases
* `return modelInstance prepared for simulation` 
```

The simulation is performed with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) using the default integrator that this package automatically selects. The simulation result is
stored inside the `filter` object and can be plotted with `ModiaPlot.plot`.

 It is also possible to specify the integrator as second argument of `simulate!`:

```julia
    using ModiaBase
    using DifferentialEquations
    using ModiaPlot

    filter = @instantiateModel(Filter)
    simulate!(filter, Tsit5(), stopTime=10.0)
    plot(filter, ["y", "x"], figure=1)
```

Integrator `DifferentialEquations.Tsit5` is an adaptive Runge-Kutta method of order 5/4.

For more information, see the documentation of [`simulate!`](@ref) and the documentation of the
[ModiaPlot](https://modiasim.github.io/ModiaPlot.jl/stable/) package.


# Appendix

## 1 Named tuples and quoted expressions

The fundamental mechanism for defining models in TinyModia are named tuples which is a list of key/value pairs enclosed in parentheses:

```julia
julia> S=(p=5, q=10)
(p = 5, q = 10)

julia> typeof(S)
NamedTuple{(:p, :q),Tuple{Int64,Int64}}
```

Named tuples are conceptually similar to dictionaries (Dict), but the constructor syntax is simpler. Note that if only one key/value pair is given, a comma must preceed the final parentheses: `(p = 5, )`. It is also possible to define named tuples using a keyword argument list, i.e. a list starting with a semi-colon: `z=(;p=5)`.

The values can also be a quoted expression, i.e. an expression enclosed in `:( )`, an array of quoted expressions encloded in `:[ ]` or just a quoted symbol, `:x`. This mechanism is used to encode equations and expressions of the model which needs to be manipulated before the model can be simulated.

Julia defines a very useful merge operation between named tuples (and dictionaries):

```julia
julia> T=(q=100, r=200)
(q = 100, r = 200)

julia> merge(S, T)
(p = 5, q = 100, r = 200)
```

If a key already exists `q` in the first named tuple, it's value is overwritten otherwise it's added, `r`. Such a merge semantics allows for unification of parameter modifications and inheritance as will be demonstrated below.

## 2 Mergemodel algorithm

The `mergeModel` algorithm is defined as follows (without logging):

```julia
function mergeModels(m1::NamedTuple, m2::NamedTuple, env=Symbol())
    mergedModels = OrderedDict{Symbol,Any}(pairs(m1)) # Convert the named tuple m1 to an OrderedDict
    for (k,v) in collect(pairs(m2))
        if typeof(v) <: NamedTuple
            if k in keys(mergedModels) && ! (:_redeclare in keys(v))
                mergedModels[k] = mergeModels(mergedModels[k], v, k)
            else
                mergedModels[k] = v
            end
        elseif v === nothing
            delete!(mergedModels, k)
        else
            mergedModels[k] = v
        end
    end
    return (; mergedModels...) # Transform OrderedDict to named tuple
end

|(m::NamedTuple, n::NamedTuple) =  mergeModels(m, n)

Redeclare = ( _redeclare = true, )
```

# Functions

```@meta
CurrentModule = Modia
```


## Instantiation

```@docs
@instantiateModel
instantiateModel
```

## Simulation

```@docs
simulate!
```


## Linearization

```@docs
linearize!
```

## Parameters/Init/Start

```@meta
CurrentModule = Modia
```

The following functions are provided to inquire values of *parameters* and of
*init/start* values before and after *evaluation*:

| Functions                         | Description                                             |
|:----------------------------------|:--------------------------------------------------------|
| [`hasParameter`](@ref)            | Return true, if a parameter/init/start name is known    |
| [`getParameter`](@ref)            | Return value of a parameter/init/start name             |
| [`getEvaluatedParameter`](@ref)   | Return value of an evaluated parameter/init/start name  |
| [`getLastValue`](@ref)            | Return last available value of a variable name          |
| [`showParameters`](@ref)          | Print the parameters and the init/start values          |
| [`showEvaluatedParameters`](@ref) | Print the evaluated parameters and init/start values    |


```@docs
hasParameter
getParameter
getEvaluatedParameter
getLastValue
showParameters
showEvaluatedParameters
```


## Results

```@meta
CurrentModule = Modia
```

The simulation result of a model `instantiatedModel` supports the functions of
[ModiaResult](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html) and
exports them, so the functions can be accessed without prefixing them with `Modia`.

The following functions are provided to access **Modia results**
(use `ìnstantiatedModel` as argument `result`, e.g. `printResultInfo(instantiatedModel)`):

| Functions                                                                                                            | Description                                               |
|:---------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------|
| [`printResultInfo` ](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaResult.printResultInfo)    | Print info of the result on stdout.                       |
| [`resultInfo`](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaResult.resultInfo)               | Return info about the result as [DataFrame](https://github.com/JuliaData/DataFrames.jl) table            |
| [`rawSignal`](https://modiasim.github.io/ModiaResult.jl/stable/AbstractInterface.html#ModiaResult.rawSignal)         | Return raw signal data given the signal name.             |
| [`getPlotSignal`](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaResult.getPlotSignal)         | Return signal data prepared for a plot package.           |
| [`signalNames`](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaResult.signalNames)             | Return all signal names.                                  |
| [`timeSignalName`](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaResult.timeSignalName)       | Return the name of the time signal.                       |
| [`hasOneTimeSignal`](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaResult.hasOneTimeSignal)   | Return true if one time signal present.                   |
| [`hasSignal`](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaResult.hasSignal)                 | Return true if a signal name is known.                    |
| [`getLastValue`](@ref)                                                                                               | Return last available value of variable name              |
| [`defaultHeading`](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaResult.defaultHeading`)      | Return default heading of a result.                       |


## Plotting

```@meta
CurrentModule = Modia
```

The simulation result of a model `instantiatedModel` supports the functions of
[ModiaResult](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html) and
exports them, so the functions can be accessed without prefixing them with `Modia`.

The following functions are provided to **define/inquire the current plot package**:

| Functions                                                                                                                      | Description                                               |
|:-------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------|
| [`@usingModiaPlot`](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaResult.@usingModiaPlot)               | Expands into `using ModiaPlot_<PlotPackageName>`          |
| [`usePlotPackage`](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaResult.usePlotPackage)                 | Define the plot package to be used.                       |
| [`usePreviousPlotPackage`](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaResult.usePreviousPlotPackage) | Define the previously defined plot package to be used.    |
| [`currentPlotPackage`](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaResult.currentPlotPackage)         | Return name defined with `usePlotPackage`                 |

The following functions are available after

1. `ENV["MODIA_PLOT"] = XXX` (e.g. in startup.jl file) or
    `usePlotPackage(XXX)` has been executed (XXX = "PyPlot", "GLMakie", "WGLMakie", "CairoMakie", "NoPlot" or "SilentNoPlot") and
2. `@usingModiaPlot` has been called,

to **plot with the currently defined plot package**
(use `ìnstantiatedModel` as argument `result`, e.g. `plot(instantiatedModel, ...)`):

| Functions                                                                                                             | Description                                               |
|:----------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------|
| [`plot`](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaPlot_PyPlot.plot)                       | Plot simulation results in multiple diagrams/figures.     |
| [`saveFigure`](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaPlot_PyPlot.saveFigure)           | Save figure in different formats on file.                 |
| [`closeFigure`](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaPlot_PyPlot.closeFigure)         | Close one figure                                          |
| [`closeAllFigures`](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaPlot_PyPlot.closeAllFigures) | Close all figures                                         |
| [`showFigure`](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#ModiaPlot_PyPlot.showFigure)           | Show figure in window (only GLMakie, WGLMakie)            |

A Modia variable `a.b.c` is identified by a String key `"a.b.c"`.
The legends/labels of the plots are automatically constructed by the
names and units of the variables. Example:

```julia
using Modia
@usingModiaPlot

instantiatedModel = @instantiateModel(...)
simulate!(instantiatedModel, ...)
plot(instantiatedModel,
     [ ("phi", "r")        ("phi", "phi2", "w");
       ("w", "w2", "phi2") "w"                ],
     heading="Matrix of plots")
```

generates the following plot:

![Matrix-of-Plots](../resources/images/matrix-of-plots.png)


## PathPlanning

```@meta
CurrentModule = Modia
```

There are some pre-defined functions to define reference paths

```@docs
PTP_path
pathEndTime
getPosition!
getPosition
getIndex
getPath
```
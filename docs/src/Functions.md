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


## Results and Plotting

The simulation result of a model `instantiatedModel` supports the abstract interface
[ModiaResult](https://modiasim.github.io/ModiaResult.jl/stable/index.html) and
exports them, so the functions can be accessed without prefixing them with `Modia`.
The following functions are provided (for details see 
[Functions of ModiaResult](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#Functions-of-ModiaResult)):

```@meta
CurrentModule = Modia
```

| Functions                        | Description                                       |
|:---------------------------------|:--------------------------------------------------|
| `@usingModiaPlot`        | Expands into `using ModiaPlot_<PlotPackageName>`          |
| `usePlotPackage`         | Define the plot package to be used.                       |
| `usePreviousPlotPackage` | Define the previously defined plot package to be used.    |
| `currentPlotPackage`     | Return name defined with `usePlotPackage`                 |
| `resultInfo`             | Return info about the result as [DataFrame](https://github.com/JuliaData/DataFrames.jl) table            |
| `printResultInfo`        | Print info of the result on stdout.                       |
| `rawSignal`              | Return raw signal data given the signal name.             |
| `getPlotSignal`          | Return signal data prepared for a plot package.           |
| `defaultHeading`         | Return default heading of a result.                       |
| `signalNames`            | Return all signal names.                                  |
| `timeSignalName`         | Return the name of the time signal.                       |
| `hasOneTimeSignal`       | Return true if one time signal present.                   |
| `hasSignal`              | Return true if a signal name is known.                    |


The following functions are available after `ENV["MODIA_PLOT"] = XXX` or
`@usingModiaPlot(XXX)` have been executed
(for details see 
[Functions of Plot Package](https://modiasim.github.io/ModiaResult.jl/stable/Functions.html#Functions-of-Plot-Package)):


| Functions          | Description                                               |
|:-------------------|:----------------------------------------------------------|
| `plot`             | Plot simulation results in multiple diagrams/figures.     |
| `saveFigure`       | Save figure in different formats on file.                 |
| `closeFigure`      | Close one figure                                          |
| `closeAllFigures`  | Close all figures                                         |
| `showFigure`       | Show figure in window (only GLMakie, WGLMakie)            |


A Modia variable `a.b.c` is identified by a String key `"a.b.c"`.
The legends/labels of the plots are automatically constructed by the
names and units of the variables. Example:

```julia
using Modia
@usingModiaPlot

instantiatedModel = @instantiatedModel(...)
simulate!(instantiatedModel, ...)
plot(instantiatedModel,
     [ ("phi", "r")        ("phi", "phi2", "w");
       ("w", "w2", "phi2") "w"                ],
     heading="Matrix of plots")
```

generates the following plot:

![Matrix-of-Plots](../resources/images/matrix-of-plots.png)


## PathPlanning

There are some pre-defined functions to define reference paths

```@docs
PTP_path
pathEndTime
getPosition!
getPosition
getIndex
getPath
```
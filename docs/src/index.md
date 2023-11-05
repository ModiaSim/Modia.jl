# Modia Documentation

[Modia](https://github.com/ModiaSim/Modia.jl) is an environment in form of a Julia package to model and simulate physical systems (electrical, mechanical, thermo-dynamical, etc.) described by differential and algebraic equations. A user defines a model on a high level with model components (like a mechanical body, an electrical resistance, or a pipe) that are physically connected together. A model component is constructed by **`expression = expression` equations** or by Julia structs/functions, such as the pre-defined [Modia3D] (https://github.com/ModiaSim/Modia3D.jl) multibody components. The defined model is symbolically processed (for example, equations might be analytically differentiated) with algorithms from package [ModiaBase.jl](https://github.com/ModiaSim/ModiaBase.jl). From the transformed model a Julia function is generated that is used to simulate the model with integrators from [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).
The basic type of the floating point variables is usually `Float64`, but can be set to any
type `FloatType <: AbstractFloat` via `@instantiateModel(..., FloatType = xxx)`, for example
it can be set to `Float32, DoubleFloat, Measurement{Float64}, StaticParticles{Float64,100}`.


## Installation

The package is registered and is installed with (Julia >= 1.7 is required):

```julia
julia> ]add Modia
```

Furthermore, one or more of the following packages should be installed in order
to be able to generate plots:

```julia
julia> ]add SignalTablesInterface_PyPlot        # if plotting with PyPlot desired

        # currently under registration
        add SignalTablesInterface_GLMakie       # if plotting with GLMakie desired
        add SignalTablesInterface_WGLMakie      # if plotting with WGLMakie desired
        add SignalTablesInterface_CairoMakie    # if plotting with CairoMakie desired
```

or call `t = getValues(instantiatedModel, "time"), y = getValues(instantiatedModel, "y")` to retrieve
the results in form of vectors and arrays and use any desired plot package for plotting, e.g., `plot(t,y)`.

Note, Modia reexports the following definitions

- `using Unitful`
- `using DifferentialEquations`
- `using SignalTables`
- and exports functions `CVODE_BDF` and `IDA` of [Sundials.jl](https://github.com/SciML/Sundials.jl).

As a result, it is usually sufficient to have `using Modia` in a model to utilize the relevant
functionalities of these packages.


## Release Notes

### Version 0.12.1

- Pull request [#170](https://github.com/ModiaSim/Modia.jl/pull/170)
  (tstops in simulation engine changed from Tuple to Vector, since otherwise incompability with new Sundials version).


### Version 0.12.0

- Improved documentation of built-in component functions.

**Non-backwards** compatible changes 

- Renamed struct `SimulationModel` to `InstantiatedModel`.
- Renamed function `get_scalar_x_segmented_value` to `copy_scalar_x_segmented_value_from_state`
- Renamed function `get_SVector3_x_segmented_value` to `copy_SVector3_x_segmented_value_from_state`
- Renamed function `get_Vector_x_segmented_value!` to `copy_Vector_x_segmented_value_from_state`
- Renamed function `add_der_x_segmented_value!` to `copy_der_x_segmented_value_to_state`
- Renamed function `add_w_segmented_value!` to `copy_w_segmented_value_to_result`


### Version 0.11.0

- Require ModiaBase 0.11.1
- Manifest.toml file removed.
- test.yml for github actions changed to use Julia 1.9.0

**Non-backwards** compatible changes 

These changes only influence models that use the new feature of built-in components.

- `_buildFunction` argument list changed (options of @instantiateModel added)



### Version 0.10.0

- Initial support of segmented simulations where the number of states can change during simulation.
  For examples, see `Modia/test/TestHeatTransfer2.jl` and models in directory `Modia3D/test/Segmented`
  (of release 0.12.0 and later). The tutorial will be updated for this feature in an upcoming version.
  

**Non-backwards** compatible changes 

These changes should usually not influence user models.

- `_buildFunction = <functionName>` changed to `_buildFunction = Par(functionName = <functionName>)` and
  changed argument list of `<functionName>`.
- `_instantiateFunction = Par(..)` changed to `_initSegmentFunction = Par(functionName = <functionName>)`
  and changed argument list of `<functionName>`.


### Version 0.9.4

- Precompile statements included (compilation of Modia package takes more time, but startup of Modia model simulations is faster).
- `@instantiateModel(..., logFile=true)`: New keyword argument `logFile` in order that log of file and line number can be
  switched off, when `@instantiateModel` is called.
- Error messages improved, when model errors result in failed evaluation of parameters.
- Log of statistics improved and included in writeSignalTable(..) of instantiatedModel.
- SignalTables.getSignalNames(..): Order of signal names improved so that the linear listing reflects the hierarchy of the names.
- writeSignalTable(..): attributes renamed to _attributes. All used simulate!(..) options included in Map experiment.
- Require SignalTables 0.4.2 (since several issues fixed with writeSignalTable(..)).
- Update to newest versions of packages.

**Bug fixes**

- DifferentialEquations 7.6.0 introduced a non-backwards compatible change with [#867](https://github.com/SciML/DifferentialEquations.jl/issues/867). Modia was corrected to cope with this change (based on [#162](https://github.com/ModiaSim/Modia.jl/pull/162)).
- `@instantiateModel(..., logCalculations=true)` skipped actual computations. This was fixed via [#161](https://github.com/ModiaSim/Modia.jl/pull/161).


### Version 0.9.3

- Requires SignalTables 0.4.0 (introduces Map-signal)
- getSignalNames(...; getVar=true, getPar=true, getMap=true): New keyword arguments to filter names.
- writeSignalTable(...) of instantiatedModel: Include attributes = Map(model=..., experiment=...).
- Some internal bug-fixes.


### Version 0.9.2

- Bug fix: integrator IDA() can be used (especially to avoid solving large linear equation systems in the model).\
  Extend some test models to use IDA().


### Version 0.9.1

- Requires SignalTables 0.3.5.

- [`@usingModiaPlot`](@ref): corrected and fixed in docu. Alternatively, @usingPlotPackage can be used,
  provided package SignalTables is present in your current environment.

- Internal: A function call in the generated code prefixed with `Modia.`.


### Version 0.9.0

- This version is slightly **non-backwards** compatible to 0.8.x. Most important, the result handling has been changed.
  Especially, package [ModiaResult.jl](https://github.com/ModiaSim/ModiaResult.jl) has been replaced by
  package [SignalTables.jl](https://github.com/ModiaSim/SignalTables.jl).
  Also the plot package interfaces SignalTablesInterface\_PyPlot, SignalTablesInterface\_GLMakie etc. have been replaced by packages
  SignalTablesInterface\_PyPlot, SignalTablesInterface\_GLMakie etc.\
  In order that plotting works again with your models, you have to add one of the new plot package interfaces, e.g.
  `]add SignalTablesInterface_PyPlot`. One benefit is, that the plot packages have now access to all attributes
  associated with a variable.

- An instantiated model (as returned from `@instantiateModel(..)`) is now a signal table according to [SignalTables.jl](https://github.com/ModiaSim/SignalTables.jl).
  This means that all functions defined for a signal table (see [function overview](https://modiasim.github.io/SignalTables.jl/stable/Functions/OverviewOfFunctions.html))
  can be applied on a simulated model. Hereby, all Var(..) and Par(..) Modia variables are seen as signals of the signal table
  (so both time varying variables, as well as parameters). See example `Modia/test/TestFirstOrder2.jl`.\
  For example, it is now possible to store simulation results (together with all parameter and start values) on file in JSON format
  with `writeSignalTable(filename, instantiatedModel)` (or in HDF5 format via [JDL](https://github.com/JuliaIO/JLD.jl)).
  You get an overview of a simulation result via `showInfo(instantiatedModel)`.

- New functions [`hasParameter`](@ref), [`getParameter`](@ref), [`getEvaluatedParameter`](@ref),
  [`showParameters`](@ref), [`showEvaluatedParameters`](@ref) to
  get parameter/init/start values by name (e.g. `getEvaluatedParameter(instantiatedModel, "a.b.c")`) or
  show all parameters.

- New functions to add states and algebraic variables from within functions that are not visible in the generated code
  (see [Variables of built-in Components](@ref) and example `Modia/test/TestLinearSystems.jl`).
  This feature is used in the next version of
  Modia3D to allow (Modia3D) model changes after code generation and to get more light weight code.

- simulate!(..): Maximum number of iterations is switched off (DifferentialEquations.jl option set to: maxiters = Int(typemax(Int32)) ≈ 2e9).

- Docu improved.


**Bug fixes**

1. A hierarchical model name with a derivative operator, say `der(a.b.c)`, has now the correct name `a.b.der(c)`
   in the result. For example, the plot command needs to be changed to `plot(..., "a.b.der(c)")` instead of the previous
   command `plot(..., "der(a.b.c)")`.
2. The initial state vector was not always correctly filled with start/init values of the model (is now fixed).
3. `signalNames(instantiatedModel)` did sometimes not show the correct signal names available in the result (is now fixed).
   `signalNames` is deprecated. Use instead [getSignalNames](https://modiasim.github.io/SignalTables.jl/stable/Functions/SignalTables.html#SignalTables.getSignalNames).


**Non-backwards compatible changes**

- Bug fix 1 can lead for some models to warnings and the selected variable is no longer plotted (-> the model needs to be changed).

- Bug fix 2 can lead for some models to a different result (without notice).

- The result data structure is now constructed with `deepcopy(..)` of every involved result variable.
  Previously, for some result variables just the variable reference was stored.
  The effect is that if previously a complex internal data structure was incorporated into the result data structure,
  then it was present just once. Now, a deepcopy of the data structure is stored at every time instant.
  Note, a variable `v` (especially, a complex internal data structure) is not stored in the result if defined as
  `v = Var(hideResult=true)`.
  In some rare cases, `deepcopy(..)` gives an error (if module variables are, for whatever reason, tried to be copied).
  Such variables `v` need to be declared with `v = Var(hideResult=true)`, in order that this error does not appear
  (and these variables are then not stored in the result).

- Function `rawSignal(instantiatedModel, name)`  is no longer supported.
  Use [getValues](https://modiasim.github.io/SignalTables.jl/stable/Functions/SignalTables.html#SignalTables.getValues)
  or [getSignal](https://modiasim.github.io/SignalTables.jl/stable/Functions/SignalTables.html#SignalTables.getSignal) instead.

- Function `getPlotSignal(instantiatedModel, name)` is no longer supported.
  Use [getFlattenedSignal](https://modiasim.github.io/SignalTables.jl/stable/Functions/SignalTables.html#SignalTables.getFlattenedSignal) instead.

- Function `getPath(path, ...)` does no longer return a dictionary but a [SignalTable](https://modiasim.github.io/SignalTables.jl/stable/Functions/SignalTables.html#SignalTables.SignalTable).


### Version 0.8.4

- Fix issue with DiffEqBase, version 6.91.6 and later.


### Version 0.8.3

- Bug fix: Parameters that are Numbers, but not AbstractFloats, and have no unit defined,
  e.g. a Bool or an Int parameter, are no longer converted to FloatType in the generated Code.


### Version 0.8.2

- New exported functions
  - modelToJSON(model; expressionsAsStrings=true)
  - JSONToModel(json)
  - writeModel(filename, model; log=true)
  - readModel(filename; log=true)
  writeModel saves the model in JSON format on file and readModel reads the model from a JSON file.

- `Modia/examples/ServoSystem.jl` enhanced, so that the hierarchical model with units is
  first stored in JSON format on file, then the model is read from file, a parameter is modified,
  and then the model is passed to @instantiateModel(..).

- `@instantiateModel(...)`:
  - `@instantiateModel(..., logCode=true, ...)` provides now correct unit type casts for scalars.
  - `@instantiateModel(..., saveCodeOnFile=fileName, ...)` stores the generated code on file `fileName`.
  - Automatically use `@instantiatedModel(..., unitless=true, ..)`, if `FloatType = MonteCarloMeasurements.XXX`,
    because there are easily cases where this fails, if units are present.

- `@showModel model`: Nicer pretty print if model is hierarchical.

- New function `Modia.unitAsString(unitOfQuantity)`,
  see Unitful issue 412 (https://github.com/PainterQubits/Unitful.jl/issues/412).

- Remove empty hierarchies in model parameters (seen with `simulate!(..., logParameters=true)`).

- Memory allocation reduced if states or tearing variables are SVectors.

- Improved casting and checking of types in the generated code
  (see new test model Modia/test/TestUnitAsString.jl).

- Moved ModiaBase.Symbolic.makeDerVar from ModiaBase to new file `Modia/src/Symbolic.jl` (because makeDerVar needs FloatType for
  generating type-stable code and FloatType is available in Modia but not in ModiaBase).

- Github actions workflow added for automatic tests on Linux/Windows/MacOS, for pull requests on main.

**Bug fixes**

- Fixed issue with unit on macOS (exponents had been displayed as Unicode superscripts when converting
  the unit to a string, leading to errors in the further processing).

- Hide result only if `Var(hideResult=true)` (previously, hideResult=false was treated as true).

- `Modia/models/Rotational.jl`: Change some Int to Float64 values, because errors occured in some situations.


### Version 0.8.1

- Missing file Modia/test/TestLinearEquations.jl added.

### Version 0.8.0

**Non-backwards** compatible changes

The Modia packages are slightly restructured to allow more efficient operations.
Previously, Modia was planned to include all the functionality with all model libraries.
This is now changed and Modia includes now equation-oriented modeling and basic model libraries.
Further model libraries, such as Modia3D (and other model libraries in the future) must be
explicitly imported and are no longer automatically imported by Modia.
To simplify the structuring, ModiaLang is merged into Modia
and some functionality for the code generation is moved from ModiaBase to Modia.
Overall, the benefit is that loading and compilation times are reduced, if Modia3D is not needed.
Furthermore, the generated code contains only references to Modia functionality and no longer to ModiaBase.
Details of the changes:

- ModiaLang#main 0.11.3 and ModiaLang#development merged into Modia 0.7.0 resulting
  in the new Modia version 0.8.0 (hereby history of both ModiaLang and of Modia is preserved).

- Modia3D is removed from Modia (so when a model is using Modia3D, the package must be explicitly imported
  and is no longer automatically imported from Modia).

- Require ModiaBase 0.10 (where EquationAndStateInfo.jl and StateSelection.jl are removed and
  added to Modia, in order that only references to Modia are in the generated code and no longer
  references to ModiaBase).


## Old Release Notes (until 28.2.2022)

### Release Notes of Modia (until 28.2.2022)

#### Version 0.7.0

**Non-backwards** compatible changes (basically, these changes are, erronously, in 0.6.1):

- Equations can only be defined with key `equations` and no other key
  (still, expressions can be associated with one variable, such as `b = Var(:(2*a))`).
  In versions 0.6.0 and before, equations could be associated with any key.

- The merge operator `|` appends the expression vectors of `equations`, so
  `m1 | m2` basically appends the vector of `m2.equations` to the vector of `m1.equations`.
  In versions 0.6.0 and before, the merge operator did not handle `equations` specially,
  and therefore `m1 | m2` replaced `m1.equations` by `m2.equations`.

- Parameter values in the code are now type cast to the type of the parameter value from the
  `@instantiatedModel(..)` call. The benefit is that access of parameter values in the code is type stable
  and operations with the parameter value are more efficient and at run-time no memory is allocated.
  Existing models can no longer be simulated, if parameter values provided via `simulate!(.., merge=xx)` are not
  type compatible to their definition. For example, an error is thrown if the @instantedModel(..) uses a Float64 value and the
  `simulate!(.., merge=xx)` uses a `Measurement{Float64}` value for the same parameter

- Operator `buildModia3D(..)` as used in Modia3D models is removed. Instead, the new constructor
  `Model3D(..)` must be used at the top level of a Modia3D definition. It is now possible to define
  several, independent multibody systems (currently, only one of them can have animation and animation export).

- `Var(init=[...])` or `Var(start=[..])` of FreeMotion joints must be defined as
  `Var(init=SVector{3,Float64}(..))` or `Var(start=SVector{3,Float64}(..))`.
  Otherwise, errors occur during compilation.


Other changes

- Documentation (especially tutorial) adapted to the new version.

- Examples and test models (Modia/examples, Modia/tests) adapted to the new version, especially
  to the non-backwards compatible changes.

- For further changes of equation-based models, see the release notes of [ModiaLang 0.11.0](https://github.com/ModiaSim/ModiaLang.jl/releases/tag/v0.11.0).

- For further changes of Modia3D models, see the release notes of [Modia3D 0.9.0](https://github.com/ModiaSim/Modia3D.jl/releases/tag/v0.9.0).


#### Version 0.6.1

This version was erronously released as 0.6.1. Since it contains non-backwards compatible
changes with respect to 0.6.0, this is wrong and should have been released as version 0.7.0.

- See release notes of [ModiaLang](https://github.com/ModiaSim/ModiaLang.jl/releases/tag/v0.11.0) and
  of [Modia3D](https://github.com/ModiaSim/Modia3D.jl/releases/tag/v0.9.0).

- Project.toml and Manifest.toml updated due to new versions of Modia3D and ModiaLang

- docu: fix some typing and formatting


#### Version 0.6.0

- Modia is restricted to Julia 1.7
- cyclic dependencies with Modia3D package are removed


#### Version 0.5.2

- Fully reexporting Modia3D and removing duplicate ModiaInterface (see [Modia3D release notes 0.6.0](https://github.com/ModiaSim/Modia3D.jl/releases/tag/v0.6.0)).

#### Version 0.5.1

- Using and reexporting ModiaLang 0.8.3 (see [release notes 0.8.3 and 0.8.2](https://github.com/ModiaSim/ModiaLang.jl/releases)).
- Using and partially reexporting Modia3D 0.5.1 (see [release notes 0.5.1](https://github.com/ModiaSim/Modia3D.jl/releases/tag/v0.5.1)).



#### Version 0.5.0

- Using and reexporting ModiaLang 0.8.1 (see [release notes](https://modiasim.github.io/ModiaLang.jl/stable/)).
- Using and partially reexporting Modia3D 0.5.0 (see [release notes](https://modiasim.github.io/Modia3D.jl/stable/#Release-Notes)).
- New plot package interface via [ModiaResult](https://github.com/ModiaSim/ModiaResult.jl). Additional support for PyPlot, WGLMakie, CairoMakie (besides GLMakie).


#### Version 0.4.0

- Initial version of new Modia design.


### Release Notes of ModiaLang (until 28.2.2022)

#### Version 0.11.3

- @instantiateModel(..): `Var(hideResult=true)` is no longer ignored if present in a sub-component.

- simulate!(..): Unnecessary evaluation of the parameters dictionary is avoided
  (if merge = missing, nothing or has no elements).


#### Version 0.11.2

- Minor (efficiency) improvement if states are SVectors.

- Require ModiaBase 0.9.2 (to get rid of performance issues in Modia3D).

- Replace ustrip(..) with ustrip.(..) at some places to get rid of warnings.


#### Version 0.11.1

- Update of Manifest.toml file
- Require ModiaBase 0.9.1 (with updated Manifest.toml file)


#### Version 0.11.0

Non-backwards compatible changes

- Equations can only be defined with key `equations` and no other key.

- Parameter values in the code are now type cast to the type of the parameter value from the
  `@instantiatedModel(..)` call. The benefit is that access of parameter values in the code is type stable
  and operations with the parameter value are more efficient and at run-time no memory is allocated.
  Existing models can no longer be simulated, if parameter values provided via `simulate!(.., merge=xx)` are not
  type compatible to their definition. For example, an error is thrown if the @instantedModel(..) uses a Float64 value and the
  `simulate!(.., merge=xx)` uses a `Measurement{Float64}` value for the same parameter

Other changes

- Hierarchical names in function calls supported (e.g. `a.b.c.fc(..)`).

- Functions can return multiple values, e.g. `(tau1,tau2) = generalizedForces(derw1, derw2)`.

- Support for StaticArrays variables (the StaticArrays feature is kept in the generated AST).
  For an example, see `ModiaLang/test/TestArrays.jl`.

- Support for Array variables (especially of state and tearing variables)
  where the dimension can change after `@instantiateModel(..)`.
  For examples, see `ModiaLang/test/TestArrays.jl` and `TestMultiReturningFunction10.jl`.

- New keyword `Var(hideResult=true)` removes variable from the result (has no effect on states, derivative of states and parameters).
  For an example, see `ModiaLang/test/TestMultiReturningFunction10.jl`

- New feature of @instantiatedModel(..): If a Model(..) has key `:_buildFunction`, call this function to merge additional code to the model.
  For details see the docu of function buildSubModels! in ModiaLang.jl.
  For examples, see `ModiaLang/test/TestMultiReturningFunction10.jl` and
  constructor `Model3D(..)` in `Modia3D/src/ModiaInterface/model3D.jl` and `Modia3D/src/ModiaInterface/buildModia3D.jl`.

- Generalized connection semantics.

- Functions converting model to/from JSON: `modelToJSON(model)`, `JSONtoModel(json_string)`

- `simulate!(..):`
   - New option `logProgress=false` in function `simulate!(..)` to print current simulation time every 5s (cpu-time).
   - If tolerance is too small, a warning is prented and it is automatically enlarged to a meaningful value
     (e.g. tolerance = 1e-8 is not useful if `FloatType=Float32`)
   - Logging improved: If log=true or logTiming=true, then timing, memory allocation and compilation time is
     reported for initialization (ths includes compilation of the generated getDerivatives(..) function).
     The remaining log shows cpu-time and memory allocation **without** initialization
     (and without the resources needed to compile getDerivatives(..)).
   - Prefix messages of the timers with "ModiaLang" or "DifferentialEquations" to more clearly see
     the origin of a message in the timer log.

- Large speedup of symbolic transformation, if function depends on many input (and output) arguments
  (includes new operator `implicitDependency(..)`).

- Included DAE-Mode in solution of linear equation system (if DAE integrator is used and all unknowns of a linear
  equation system are part of the DAE states, solve the linear equation system during continuous integration
  via DAE solver (= usually large simulation speed-up, for larger linear equation systems)

Bug fixes

- If unitless=true, units in instantiatedModel.evaluatedParameters are removed.

- The unit macro is kept in the generated code and is no longer expanded. For example, `u"N"`, is kept in the code that is
  displayed with `logCode=true` (previously, this was expanded and the unit was displayed in the code as `N` which is not correct Julia code).

- Function `ModiaLang.firstInitialOfAllSegments(..)` now correctly returns true for the first call of the getDerivatives function during the simulation.


#### Version 0.10.2

- Minor (efficiency) improvement if states are SVectors.
- Require ModiaBase 0.9.2 (to get rid of performance issues in Modia3D).
- Replace ustrip(..) with ustrip.(..) at some places to get rid of warnings.


#### Version 0.10.1

- Update of Manifest.toml file
- Require ModiaBase 0.9.1 (with updated Manifest.toml file).


#### Version 0.10.0

- Require DifferentialEquations.jl version 7.
- Cleanup of using/export
- Cleanup of Project.toml/Manifest.toml.´
- @reexport using Unitful
- @reexport using DifferentialEquations
- Cleanup of test files (besides ModiaLang, no other package needed in the environment to run the tests).
- Change `InstantiatedModel{FloatType,ParType,EvaluatedParType,TimeType}` to `InstantiatedModel{FloatType,TimeType}`


#### Version 0.9.1

- New function plotPath to plot a PTP_path
- Replace ustrip(..) with ustrip.(..) at some places to get rid of warnings.
- Include time in error message, if simulation failed


#### Version 0.9.0

- Require Julia 1.7
- Upgrade Manifest.toml to version 2.0
- Update Project.toml/Manifest.toml



#### Version 0.8.7

- Packages used in test models, prefixed with ModiaLang. to avoid missing package errors.
- Deactivating test with DoubleFloats, since not in Project.toml
- Version/date updated



#### Version 0.8.6

- Require ModiaResult, version 0.3.9
- Project.toml/Manifest.toml updated


#### Version 0.8.5

- simulate!(..):
  - Trigger an error, if simulation is not successful (retcode is neither :Default nor :Success nor :Terminate)
  - Use RightRootFind for zero crossings (improves state events based on new DifferentialEquations option)
  - New keyword argument requiredFinalStates_atol=0.0.
  - Improve docu (e.g. add return argument solution).
  - Show correct integrator name QBDF in simulation log (instead of QNDF)
  - Raise an error, if (relative) tolerance is too small for FloatType
  - Use FloatType for zero crossing hysteresis, instead of Float64
  - If log=true print info about end of initialization.

- Support of MonteCarloMeasurements with units + new test model TestLinearEquationSystemWithUnitsAndMonteCarlo.jl

- Fixing and activating the deactivated test TestTwoInertiasAndIdealGearWithUnitsAndMonteCarlo.jl.


#### Version 0.8.4

- FloatType is included in the name space of Core.eval when evaluating parameters.

- Version and Date updated

- Included Version in printout of runtests.jl and runtests_withPlot.jl

- Print difference of finalStates and requiredFinalStates in case they do not match with the given tolerance.


#### Version 0.8.3

- Project.toml, Manifest.toml updated: Require newest version 0.7.7 of ModiaBase (containing a bug fix)

- Minor correction of simulate!(log=true) output


#### Version 0.8.2

- Issue with tearing fixed: Variables are only explicitly solved, if linear factor is a non-zero literal number
  (previously a division by zero could occur, if the linear factor became zero during simulation).

- Issue with unit of tearing variable fixed, if it is a derivative of a variable
  (previously, the generated code for unitless=false was wrong, if the tearing variable was
   a derivative, since the unit was not taken into account).

- simulate!(..):
  - Support DAE integrators, especially IDA() from Sundials.
  - New keyword `useRecursiveFactorizationUptoSize=0`: Linear equation systems A*v=b are solved with
    [RecursiveFactorization.jl](https://github.com/YingboMa/RecursiveFactorization.jl) instead of
    the default `lu!(..)` and `ldiv!(..)`, if
    `length(v) <= useRecursiveFactorizationUptoSize`.
    According to `RecursiveFactorization.jl` docu, it is faster as `lu!(..)` with OpenBLAS,
    for `length(v) <= 500` (typically, more as a factor of two).
    Since there had been some cases where `lu!(..)!` was successful,
    but `RecursiveFactorization.jl` failed due to a singular system, the default is to use `lu!(..)!`.
  - If log=true, sizes of linear equation systems are listed, as well as whether
    RecursiveFactorization.jl is used for the respective system.

- Test for RecursiveFactorization.jl added in TestTwoInertiasAndIdealGear.jl

- Some test models corrected (since leading to errors with the above changes).

- Updated Project.toml and Manifest.toml with newest versions of packages
  (including MonteCarloMeasurements, version >= 1)
  and improved Project.toml file to reduce issues with package constraints


#### Version 0.8.1

- Added a minimal documentation, including release notes.
- No message anymore, when ModiaLang is started.
- Fixed bug that `using ModiaResult` is needed, when calling `@usingModiaPlot`.


#### Version 0.8.0

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


#### Version 0.7.3

- Evaluation and propagation of parameter expressions (also in simulate!(..., merge=Map(...))).
  Propagation of start/init values of states is not yet supported.

- State events supported.


#### Version 0.7.2

- Missing dependency of Test package added.

#### Version 0.7.1

- Variable constructor `Var(...)` introduced. For example:
  `v = input | Var(init = 1.2u"m")`.

- Functions are called in the scope where macro `@instantiateModel` is called.

- New arguments of function `simulate!`:
  - Parameter and init/start values can be changed with argument `merge`.
  - A simulation can be checked with argument `requiredFinalStates`.
  - Argument `logParameters` lists the parameter and init/start values used for the simulation.
  - Argument `logStates` lists the states, init, and nominal values used for the simulation.

- `end` in array ranges is supported, for example `v[2:end]`.

- New (small) model library `Modia/models/HeatTransfer.jl`.

- Modia Tutorial improved.

- Functions docu improved.

#### Version 0.7.0

- Initial version, based on code developed for Modia 0.6 and ModiaMath 0.6.


## Main developers

- [Hilding Elmqvist](mailto:Hilding.Elmqvist@Mogram.net), [Mogram](http://www.mogram.net/).

- [Martin Otter](https://rmc.dlr.de/sr/en/staff/martin.otter/),
  [DLR - Institute of System Dynamics and Control](https://www.dlr.de/sr/en).

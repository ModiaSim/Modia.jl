# Internal

This chapter documents internal functions that are typically only
for use of the developers of a model library or of Modia.

## Variables of built-in Components

```@meta
CurrentModule = Modia
```

The following functions are provided to define and access new variables
in built-in components (seee for example model `InsulatedRod2` in `Modia/models/HeatTransfer.jl`).

| Functions                                               | Description                                                                       |
|:--------------------------------------------------------|:----------------------------------------------------------------------------------|
| [`new_x_segmented_variable!`](@ref)                     | Generate new state variable (`x_segmented` and `der_x_segmented` variables)       |
| [`new_w_segmented_variable!`](@ref)                     | Generate new local variable (`w_segmented` variable)                              | 
| [`new_alias_segmented_variable!`](@ref)                 | Generate new alias variable                                                       |
| [`new_z_segmented_variable!`](@ref)                     | Generate new zero crossing variables (`z_segmented` variables)                    |
| [`get_x_startIndex_from_x_segmented_startIndex`](@ref)  | Return start index of `x_segmented` variable with respect to state vector `x`     |
| [`copy_scalar_x_segmented_value_from_state`](@ref)      | Return value of scalar `x_segmented` variable from state vector `x`               |
| [`copy_SVector3_x_segmented_value_from_state`](@ref)    | Return value of `SVector{3,FloatType}` x_segmented variable from state vector `x` |
| [`copy_Vector_x_segmented_value_from_state`](@ref)      | Return value of `Vector{FloatType}` x_segmented variable from state vector `x`    |
| [`copy_der_x_segmented_value_to_state`](@ref)           | Copy value of `der_x_segmented` variable to state derivative vector `der(x)`      |
| [`copy_w_segmented_value_to_result`](@ref)              | Copy value of local variable (`w-segmented`) to result                            |


```@docs
new_x_segmented_variable!
new_w_segmented_variable!
new_alias_segmented_variable!
new_z_segmented_variable!
get_x_startIndex_from_x_segmented_startIndex
copy_scalar_x_segmented_value_from_state
copy_SVector3_x_segmented_value_from_state
copy_Vector_x_segmented_value_from_state
copy_der_x_segmented_value_to_state
copy_w_segmented_value_to_result
```

## Inquiries in built-in Components

The following functions are provided to inquire properties 
in built-in components at the current state of the simulation
(see for example model `InsulatedRod2` in `Modia/models/HeatTransfer.jl`).

```@docs
isInitial
isFirstInitialOfAllSegments
isTerminal
isTerminalOfAllSegments
isEvent
isFirstEventIteration
isFirstEventIterationDirectlyAfterInitial
isFullRestart
isAfterSimulationStart
isZeroCrossing
storeResults
getTime
```

## Code Generation

This section lists internal functions to **generate Julia code** of the
transformed equations.

```@meta
CurrentModule = Modia
```

```@docs
InstantiatedModel
generate_getDerivatives!
init!
outputs!
terminate!
derivatives!
DAEresidualsForODE!
affectEvent!
zeroCrossings!
affectStateEvent!
timeEventCondition!
affectTimeEvent!
addToResult!
getFloatType
measurementToString
```














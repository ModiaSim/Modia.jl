# Internal

This chapter documents internal functions that are typically only
for use of the developers of package Modia.

## Code Generation

This section provides functions to **generate Julia code** of the
transformed equations.

```@meta
CurrentModule = Modia
```

```@docs
SimulationModel
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

## Inquiries in Model

The functions in this section can be called in the model code or in 
functions that are called from the model code.

```@docs
isInitial
isTerminal
isEvent
isFirstEventIteration
isFirstEventIterationDirectlyAfterInitial
isAfterSimulationStart
isZeroCrossing
storeResults
getTime
```

## Variable definitions in functions

The following functions can be used to define states and algebraic variables inside functions:

```@docs
new_x_segmented_variable!
new_w_segmented_variable!
new_alias_segmented_variable!
new_z_segmented_variable!
get_x_startIndex_from_x_segmented_startIndex
get_scalar_x_segmented_value
get_SVector3_x_segmented_value
get_Vector_x_segmented_value!
add_der_x_segmented_value!
add_w_segmented_value!
```















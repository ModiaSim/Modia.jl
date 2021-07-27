# Internal

This chapter documents internal functions that are typically only
for use of the developers of package Modia.

## Code Generation

This section provides functions to **generate Julia code** of the
transformed equations.

```@meta
CurrentModule = ModiaLang
```

```@docs
SimulationModel
generate_getDerivatives!
init!
outputs!
addToResult!
getFloatType
baseType
measurementToString
```

## Inquiries in Model

```@docs
isInitial
isTerminal
isEvent
isFirstEventIteration
isFirstEventIterationDirectlyAfterInitial
isAfterSimulationStart
isZeroCrossing
storeResults
```

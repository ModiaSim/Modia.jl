# License for this file: MIT (expat)
# Copyright 2020-2021, DLR Institute of System Dynamics and Control
# Author: Martin Otter, DLR-SR
#
# Provide information about the structure of the ODE or DAE equations.

import DataFrames
import OrderedCollections
using  LinearAlgebra
using  StaticArrays
import TimerOutputs
import RecursiveFactorization


export StateCategory, XD, XALG, XLAMBDA, XMUE
export ResidualCategory, FD, FC_ALG, FC_LOW_HIGH, FC_MUE

export LinearEquations

export EquationInfo, StateElementInfo, get_stateNames, get_x_table


"""
    isFixedLengthStartOrInit(startOrInit::Any, name::AbstractString)

Returns true, if `startOrInit` (so start value, init value, `nothing`) characterizes a value
with fixed length, so `length(v_startOrInit)` is fixed after compilation (= either a number or a StaticArray).
Note, if a start-value is not defined (`startOrInit = nothing`) a scalar default is used.

An error is triggered, if `startOrInit` is neither a `Number` nor an `AbstractVector` nor `nothing`.
In this case argument `name` is used in the error message (currently, EquationAndStateInfo is restricted
to numbers and vectors. In a subsequent version, two and more-dimensional arrays might be supported).

# Examples
```
isFixedLengthStartOrInit(1.0)                 # = true
isFixedLengthStartOrInit(SVector{3}(1,2,3))   # = true
isFixedLengthStartOrInit([1,2,3])             # = false (size can be changed after compilation)
isFixedLengthStartOrInit(nothing)             # = true (is assumed to be a scalar with value zero)
isFixedLengthStartOrInit(missing)             # error
```
"""
isFixedLengthStartOrInit(startOrInit, name::AbstractString) = begin
    if !(startOrInit isa Nothing || startOrInit isa Number || startOrInit isa AbstractVector)
        error("Start or init value $startOrInit of variable $name\nis neither a subtype of Number, subtype of AbstractVector, nor is nothing")
    end
    !(startOrInit isa AbstractVector && !(startOrInit isa StaticArray))
end


"""
    @enum StateCategory

Type of a state.

| value     | description                                                                       |
|:----------|:----------------------------------------------------------------------------------|
| `XD`      | Differential variable (derivative of variable appears in equations)               |
| `XALG`    | Algebraic variable that is included in x                                          |
| `XLAMBDA` | Algebraic variable that is included in der_x                                      |
| `XMUE`    | Algebraic variable used for stabilization and included in der_x (exact value = 0) |
"""
@enum StateCategory XD=1 XALG XLAMBDA XMUE


"""
    @enum ResidualCategory

Type of a residual (see also [`StateCategory`](@ref) for the type of a state).

| value         | description                                                              |
|:--------------|:-------------------------------------------------------------------------|
| `FD`          | Differential equation (depends on XD, XLAMBDA, XMUE variables)           |
| `FC_ALG`      | Constraint equation that is purely algebraic (depends on XALG variables) |
| `FC_LOW_HIGH` | Constraint equation on lowest and on highest derivative level            |
| `FC_MUE`      | Constraint equation that is associated with a stabilizer XMUE            |
"""
@enum ResidualCategory FD=1 FC_ALG FC_LOW_HIGH FC_MUE



const niter_max = 20  # Maximum number of fixed-point iterations to solve A*x = b

"""
    leq = LinearEquations{FloatType}(x_names::Vector{String}, x_vec_julia_names::AbstractVector,
                                     x_lengths::Vector{Int}, nx_fixedLength::Int, A_is_constant::Bool;
                                     useRecursiveFactorizationUptoSize = 0)

Define linear equation system "A*x=b" with `x::Vector{FloatType}`.

- `x` is constructed from a set of scalar or vector variables i with names `x_names[i]`
   and length `x_length[i]` (that is `length(x) = sum(x_lengths)`).
  x_names[1:nx_fixedLength] are elements with fixed lengths (dimensions do not change after compilation).
  x_names[nx_fixedLength+1:end] are vector-valued elements where the dimensions may changer after compilation.

- `x_vec_julia_names` are the Julia names of the vector-valued elements.

- If `A_is_constant = true` then `A` is a matrix that is constant after initialization

- If length(x) <= useRecursiveFactorizationUptoSize, then linear equation systems will be solved with
  `RecursiveFactorization.jl` instead of the default `lu!(..)` and `ldiv!(..)`.

For details how to use this constructor, see [`LinearEquationsIteration!`](@ref).
"""
mutable struct LinearEquations{FloatType <: Real}
    odeMode::Bool                     # Set from the calling function after LinearEquations was instantiated (default: true)
                                      # = true (standard mode): Compute "x" from equation "residuals = A*x - b"
                                      # = false (DAE solver): During events (including initialization)
                                      #   compute "x" as for odeMode=true. Outside of events:
                                      #   (1) "x" is set from the outside (= der(x_dae) provided by DAE solver)
                                      #   (2) Compute "residuals" from "residuals := A*x - b"
                                      #   (3) From the outside copy "residuals" into "residuals_dae" of the DAE solver.

    A_is_constant::Bool               # = true, if A-matrix is constant
    x_names::Vector{String}           # Names of the x-variables
    x_vec_julia_names::AbstractVector # Julia names of the vector-valued x-variables
                                      # (needed to generate code for
                                      #     if _m.storeResult
                                      #        x_vec_julia_names[1] = deepcopy(x_vec_julia_names[1])
                                      #        ...)
    x_lengths::Vector{Int}            # Lengths of the x-variables (sum(x_lengths) = length(x))
    nx_fixedLength::Int               # x_names[1:nx_fixedLength] are elements with fixed lengths (after compilation)
                                      # x_names[nx_fixedLength+1:end] are vector elements where the size may change after compilation
    x_vec::Vector{Vector{FloatType}}  # x_vec[1] is the value vector of x_names[nx_fixedLength+1]
                                      # x_vec[2] is the value vector of x_names[nx_fixedLength+2]
                                      # <..>


    A::Matrix{FloatType}
    b::Vector{FloatType}
    x::Vector{FloatType}              # Values of iteration variables
    pivots::Vector{Int}               # Pivot vector if recursiveFactorization = true
    residuals::Vector{FloatType}      # Values of the residuals FloatType vector; length(residuals) = length(x)

    # Iteration status of for-loop
    mode::Int       # Operational mode (see function LinearEquationsIteration!)
    niter::Int      # Current number of iterations in the fix point iteration scheme
    niter_max::Int  # Maximum number of iterations
    success::Bool   # = true, if solution of A*x = b is successfully computed
                    # = false, if solution is not computed; continue with fix point iteration

    # For warning message if niter > niter_max
    inconsistentPositive::Vector{String}
    inconsistentNegative::Vector{String}

    # Constructed during initialization
    useRecursiveFactorizationUptoSize::Int
    useRecursiveFactorization::Bool             # = true, if RecursiveFactorization.jl shall be used to solve the linear equation system

    luA::LU{FloatType,Array{FloatType,2}}       # lu-Decomposition of A

    function LinearEquations{FloatType}(x_names::Vector{String}, x_vec_julia_names::AbstractVector,
                                        x_lengths::Vector{Int}, nx_fixedLength::Int, A_is_constant::Bool;
                                        useRecursiveFactorizationUptoSize::Int = 0) where {FloatType <: Real}
        @assert(length(x_names) > 0)
        @assert(length(x_names) == length(x_lengths))
        nx = sum(x_lengths)
        @assert(nx >= 0)
        useRecursiveFactorization = nx <= useRecursiveFactorizationUptoSize

        # Allocate value storage for vector elements
        nx_vec = length(x_names) - nx_fixedLength
        x_vec  = fill(FloatType[], nx_vec)
        j      = nx_fixedLength
        for i = 1:nx_vec
            j += 1
            x_vec[i] = zeros(FloatType, x_lengths[j])
        end

        new(true, A_is_constant, x_names, Any[], x_lengths, nx_fixedLength, x_vec,
            zeros(FloatType,nx,nx), zeros(FloatType,nx), zeros(FloatType,nx), fill(0,nx), zeros(FloatType,nx),
            -2, 0, niter_max, false, String[], String[],
            useRecursiveFactorizationUptoSize, useRecursiveFactorization)
    end
end
LinearEquations(args...) = LinearEquations{Float64}(args...)


"""
    b = copy_x_into_x_vec!(leq)

If vector valued elements of x-vector, copy vector valued elements of leq.x
into vector leq.x_vec. The function returns true
"""
function copy_x_into_x_vec!(leq::LinearEquations{FloatType})::Bool where {FloatType}
    i_vec = 0
    i_x   = leq.nx_fixedLength
    for i = leq.nx_fixedLength+1:length(leq.x_names)
        i_vec += 1
        x_vec = leq.x_vec[i_vec]
        for j = 1:leq.x_lengths[i]
            x_vec[j] = leq.x[i_x+j]
        end
        i_x += leq.x_lengths[i]
    end
    return true
end



"""
     iterating = LinearEquationsIteration!(leq, isInitial, [solve, isStoreResult,] time, timer)

This function solves a linear equation system in residual form "residual = A*x - b"
by iterating with a while loop over this system (arguments `solve, isStoreResult` are optional
and have only an effect if leq.odeMode=false, that is if the linear equation system is
solved in DAE mode):

```julia
function getDerivatives!(_der_x, _x, _m, _time)::Nothing
    _leq::Union{Nothing,LinearEquations{FloatType}} = nothing
    ...
    _leq      = _m.linearEquations[<nr>]   # leq::LinearEquations{FloatType}
    _leq.mode = -3  # initializes the iteration
    while LinearEquationsIteration!(_leq, _m.isInitial, _m.solve_leq, _m.storeResult, _m.time, _m.timer)
        x1 = _leq.x[1]
        x2 = _leq.x[2]
        x3 = _leq.x_vec[1]
        ...
        v_solved1 = f(x1, x2, ..., positive(x1,.., _leq))
        v_solved2 = f(x1, x2, ..., v_solved1)
        ...
        appendResidual!(_leq.residuals, < getResiduals(x1,x2,...) > )
    end
    ...
end
```

Note, A and b can be functions of linear event operators, such as `positive(c_i'*x - d_i)`.
In this case this system is solved by a fixed point iteration scheme,
that is `A*x = b` is solved until potentially present
`positive(c_i'*x - d_i)` calls are consistent to x
(`c_i'*x - d_i` and positive(..) must have the same sign).
Iteration is terminated after 20 iterations if no convergence is found.
In this case a warning message is triggered and simulation is continued
(it might be that simulation is still successful, even if `x` and
`positive(..)` are temporarily not consistent to each other).

The current values of A,x,b,residuals are stored in leq.
If A is fixed after initialization (leq.A_is_constant = true), then A is computed
only once at initialization, the LU decomposition of A is stored in leq
and used in subsequent calls to solve the equation system.


# Input arguments

- `leq::LinearEquations{FloatType}`: Instance of `LinearEquations`.
- `isInitial::Bool`: = true: Called during initialization.
- `solve::Bool`: = true: leq.x is computed by LinearEquationsIteration!.
                 = false: leq.x has been set by calling environment
                          (typically when called from DAE integrator).
                          Note, at events and at terminate, solve must be true).
- `isStoreResult::Bool`: = true: Called at a communication point (store results).
- `time`: Actual time instant (only used for warning/error messages).
- `timer::TimerOutputs.TimerOutput`: Timer used to measure the solution of the linear equation system.


# Output arguments

- iterating::Bool: = true : Continue iteration of while-loop.
                   = false: Immediately break while-loop after this function was terminated.


# Enhancing efficiency

Variable `leq.mode::Int` can be accessed (read-only) in the body of the while-loop to enhance efficiency.
For example, if residuals are computed in a function and the function evaluation is expensive.

```
leq.mode = -2: Residuals might be computed, but they are not used.
leq.mode = -1: Compute "residuals .= A*x - b"
leq.mode =  0: Compute "residuals .= A*0 - b"
leq.mode >  0: Compute "residuals .= A*e_i - b"   # e_i = unit vector with i = leq.mode
```


# Hidden argument `leq.mode::Int` on input

```julia
leq.mode = -3  # LinearEquationsIteration! is called the first time
               if leq.odeMode || solve
                  # ODE mode or DAE mode at an event (solve "x" from equation "residuals = A*x - b")
                  # Initialize fixed point iteration or continue fixed point iteration (if in DAE mode)
                  leq.niter     = 0      # Number of event iterations
                  leq.success   = false  # Event iteration was not yet successful
                  leq.x        .= 0
                  leq.mode      = 0      # Compute "residuals .= A*0 - b"
               else # solve = false
                  # DAE mode (but not at an event)
                  leq.niter     = 0
                  leq.success   = true
                  leq.x_is_zero = false
                  if isStoreResult
                     leq.mode = -2  # Residuals might be computed, but they are not used.
                  else
                     leq.mode = -1  # Compute "residuals .= A*x - b"
                  end
               end
               return true  # Continue while-loop

leq.mode = -2  # Terminate while-loop or initialize next event iteration
               if leq.success           # Set from leq.mode=-2 or leq.mode>=0 or from event operators, such as positive!(...)
                  return false          # Terminate while-loop
               elseif leq.niter > 20
                  <warning>             # Event iteration failed (variables are not consistent)
                  return false          # Terminate while-loop
               else
                  # Initialize next event iteration
                  leq.x        .= 0
                  leq.x_is_zero = true
                  leq.mode      = 0     # Compute "residuals .= A*0 - b"
                  return true           # Continue while-loop
               end

leq.mode = -1: @assert(!leq.odeMode && !solve)
               # DAE mode, but not at an event
               return false      # Terminate while-loop

leq.mode =  0: @assert(leq.odeMode || solve)
               # ODE mode or DAE mode at an event (solve "x" from equation "residuals = A*x - b")
               leq.b    = -leq.residuals
               leq.x[1] = 1.0
               leq.mode = 1      # Compute "residuals := A*e_1 - b"
               return true       # Continue while-loop

leq.mode >  0: @assert(leq.odeMode || solve)
               j          = leq.mode
               leq.A[:,j] = leq.residuals + leq.b
               leq.x[j]   = 0.0
               if j < nx
                   leq.x[j+1] = 1.0
                   leq.mode   = j+1      # Compute A[:,j+1] after the next iteration
                   return true           # Continue while-loop
               elseif j == nx
                   leq.x = solve(leq.A,leq.b) # Solve linear equation system A*x = b for x.
                   leq.mode    = -1      # Compute all variables IN the next iteration as function of x
                   leq.niter  += 1       # Increment number of iterations to solve A*x = b
                   leq.success = true    # Terminate for-loop at the beginning of the next iteration,
                                         # provided no positive(..) call changes its value.
                                         # (leq.success is set to false in positive(..), if the return value changes).
                   return true           # Continue while-loop
               end
```
"""
LinearEquationsIteration!(leq::LinearEquations, isInitial, time, timer) = LinearEquationsIteration!(leq, isInitial, true, false, time, timer)
function LinearEquationsIteration!(leq::LinearEquations{FloatType}, isInitial::Bool, solve::Bool,
                                   isStoreResult::Bool, time, timer)::Bool where {FloatType}
    mode = leq.mode
    nx   = length(leq.x)

    if mode == -3
        # LinearEquationsIteration! is called the first time in the current model evaluation
        leq.niter = 0      # Number of event iterations
        empty!(leq.inconsistentPositive)
        empty!(leq.inconsistentNegative)

        if leq.odeMode || solve
           # ODE mode or DAE mode at an event (solve "x" from equation "residuals = A*x - b")
           # Initialize fixed point iteration or continue fixed point iteration (if in DAE mode)
           leq.success  = false  # Event iteration was not yet successful
           leq.x       .= 0
           leq.mode     = 0      # Initialize leq and compute "residuals .= A*0 - b"
        else
           # DAE mode (but not at an event)
           leq.success   = true
           if isStoreResult
              leq.mode = -2   # Residuals might be computed, but they are not used.
           else
              leq.mode = -1   # Initialize leq and compute "residuals .= A*x - b"
           end
        end
        empty!(leq.residuals)
        return copy_x_into_x_vec!(leq)  # Continue while-loop

    elseif mode == -2
        # Terminate while-loop or initialize next event iteration
        if leq.success   # Set from leq.mode=-2 or leq.mode>=0 or from event operators, such as positive!(...)
            return false # Terminate while-loop
        elseif leq.niter > niter_max
            str = ""
            if length(leq.inconsistentPositive) > 0
                str = str * "positive(expr) is inconsistent for expr = $(leq.inconsistentPositive)."
            end
            if length(leq.inconsistentNegative) > 0
                if length(str) > 0
                    str = str * "\n"
                end
                str = str * "negative(expr) is inconsistent for expr = $(leq.inconsistentNegative)."
            end
            @warn "At time = $time, no consistent solution found for mixed linear equation system.\n" *
                  "Simulation is continued although some variables might not be correct at this time instant.\n$str"
            return false       # Terminate while-loop
        end
        leq.x        .= 0
        leq.mode      = 0      # Compute "residuals .= A*0 - b"
        empty!(leq.residuals)
        return copy_x_into_x_vec!(leq) # Continue while-loop

    elseif mode < -3 || mode > nx
        @goto ERROR
    end

    x = leq.x
    A = leq.A
    b = leq.b
    residuals = leq.residuals

    if length(residuals) != length(x)
        error("Function LinearEquationsIteration! wrongly used:\n",
              "length(leq.residuals) = ", length(leq.residuals), ", length(leq.x) = ", length(leq.x))
    end

    if mode == -1
        @assert(!leq.odeMode && !solve)
        return false   # Terminate while-loop (leq.residuals must be copied into DAE residuals)

    elseif !leq.A_is_constant || isInitial  # A is not constant or A is constant and isInitial = true
        if mode == 0
            # residuals = A*x - b -> b = -residuals)
            for i = 1:nx
                b[i] = -residuals[i]
            end
            leq.mode = 1
            x[1] = max(FloatType(1e-3), abs(b[1]))     #convert(FloatType, 1)
            empty!(leq.residuals)
            return copy_x_into_x_vec!(leq)
        end

        # residuals = A*e_j*x_j - b -> A[:,j] = (residuals + b)/x[j]
        j = mode
        for i = 1:nx
            A[i,j] = (residuals[i] + b[i])/x[j]
        end
        x[j] = 0

        if j < nx
            leq.mode += 1
            x[leq.mode] = max(FloatType(1e-3), abs(b[leq.mode]))    # convert(FloatType, 1)
            empty!(leq.residuals)
            return copy_x_into_x_vec!(leq)
        end

        # Solve linear equation system
        if nx == 1
            x[1] = b[1]/A[1,1]
            if !isfinite(x[1])
                error("Linear scalar equation system is singular resulting in: ", leq.x_names[1], " = ", x[1])
            end
        else
            x .= b
            if leq.useRecursiveFactorization
                Modia.TimerOutputs.@timeit timer "Modia LinearEquationsIteration! (solve A*x=b Rec.Fac.)" begin
                    leq.luA = RecursiveFactorization.lu!(A, leq.pivots)
                    ldiv!(leq.luA, x)
                end
            else
                Modia.TimerOutputs.@timeit timer "Modia LinearEquationsIteration! (solve A*x=b)" begin
                    leq.luA = lu!(A)
                    ldiv!(leq.luA, x)
                end
            end
        end

    elseif leq.A_is_constant && !isInitial  # isInitial=false, LU decomposition of A is available in leq.luA
        for i = 1:nx
            x[i] = -residuals[i]
        end

        # Solve linear equation system
        if nx == 1
            x[1] = x[1]/A[1,1]
            if !isfinite(x[1])
                error("Linear scalar equation system is singular resulting in: ", leq.x_names[1], " = ", x[1])
            end
        else
            ldiv!(leq.luA, x)
        end

    else
        @goto ERROR
    end

    # Linear equation system solved
    leq.mode    = -2   # Compute all variables IN the next iteration as function of x
    leq.niter  += 1    # Increment number of iterations to solve A*x = b
    leq.success = true # Terminate for-loop at the beginning of the next iteration,
                       # provided no positive(..) call changes its value.
                       # (leq.success is set to false in positive(..), if the return value changes).

    empty!(leq.residuals)
    return copy_x_into_x_vec!(leq)

    @label ERROR
    error("Should not occur (Bug in file Modia/src/EquationAndStateInfo.jl):\n",
          "   time           = $time\n",
          "   isInitial      = $isInitial,\n",
          "   solve          = $solve,\n",
          "   isStoreResult  = $isStoreREsult,\n",
          "   leq.odeMode    = $(leq.odeMode),\n",
          "   leq.mode       = $(leq.mode),\n",
          "   leq.x_names    = $(leq.x_names),\n",
          "   leq.x_lengths  = $(leq.x_lengths),\n")
end


"""
    xe_info = StateElementInfo(x_name, x_name_julia, der_x_name, der_x_name_julia,
                               stateCategory, unit, startOrInit, fixed, nominal, unbounded;
                               startIndex=-1, x_segmented_startIndex=-1)

Return an instance of the mutable struct `StateElementInfo` that defines the information
for one element of the state vector.

# Arguments
- x_name: Name of x-element or "" if no name (if stateCatebory = XLAMBDA or XMUE)
- x_name_julia: Julia name of x-element in getDerivatives!/getResiduals! function
  or "" if not needed (since no code generation).
- der_x_name: Name of der_x-element or "" if either `der(x_name)` or if no name
  (if stateCategory = XALG).
- der_x_name_julia: Julia name of der_x-element in getDerivatives!/getResiduals! function
  or "" if not needed (since no code generation).
- stateCategory::StateCategory: Category of the state
- unit: unit of XD, XALG (x_name) or XLAMBDA, XMUE (der_x_name) or "" if not yet known)
- startOrInit: Start or init value (scalar or vector)
- fixed: false (= guess value) or true (= not changed by initialization).
         Only relevant for ode=false, otherwise ignored.
- nominal: Nominal value (NaN if determined via startOrInit value)
- unbounded: false or true
- startIndex: = -1, if not yet known. startIndex with respect to x.
- x_segmented_startIndex: =-1, if not yet known. startIndex with respect to x_segmented.
"""
mutable struct StateElementInfo
    x_name::String                # Modia name of x-element or "" if no name
    x_name_julia                  # Julia name of x-element in getDerivatives! function
                                  # or "", if not needed (since no code generation).
    der_x_name::String            # Modia name of der_x-element or "" if either "der(x_name)" or if no name,
    der_x_name_julia              # Julia name of der_x-element in getDerivatives! function
                                  # or not needed (since no code generation)
    stateCategory::StateCategory  # category of the state
    unit::String                  # unit of x-element as string (or "" if not yet known)
    startOrInit::Any              # start or init value or nothing depending on fixed
                                  # (init   : if value and fixed=true
                                  #  start  : if value and fixed=false
                                  #  nothing: neither start nor init value; fixed=false)
    fixed::Bool                   #
    nominal::Float64              # nominal value (NaN if determined via start value)
    unbounded::Bool               # false or true
    fixedLength::Bool             # = true, if length of value is fixed at compile time (= scalar or StaticArray)
                                  # = false, if length of value is determined before initialization
    scalar::Bool                  # = true, if scalar
                                  # = false, if vector
    length::Int                   # length of x-element or -1 if not yet known
    startIndex::Int               # start index of state with respect to x-vector or -1 if not yet known
    x_segmented_startIndex::Int   # start index of segmented state with respect to x_segmented vector
                                  # or -1, if it is no segmented state (for a segmented state, x_segmented_startIndex
                                  # is consistently set when it is added via newHiddenState(..)).
end



# Constructor for code-generation
StateElementInfo(x_name, x_name_julia, der_x_name, der_x_name_julia,
                 stateCategory, unit, startOrInit, fixed, nominal, unbounded;
                 startIndex = -1, x_segmented_startIndex=-1) = StateElementInfo(
                 x_name, x_name_julia, der_x_name, der_x_name_julia,
                 stateCategory, unit, deepcopy(startOrInit), fixed, nominal, unbounded,
                 isFixedLengthStartOrInit(startOrInit, x_name), !(startOrInit isa AbstractArray),
                 startOrInit isa Nothing ? 1 : length(startOrInit), startIndex, x_segmented_startIndex)

function Base.show(io::IO, xe_info::StateElementInfo)
    print(io, "Modia.StateElementInfo(")
    show( io, xe_info.x_name)
	#print(io, ",")
    #show( io, xe_info.x_name_julia)
	print(io, ",")
    show( io, xe_info.der_x_name)
	#print(io, ",")
    #show( io, xe_info.der_x_name_julia)
    print(io, ",Modia.", xe_info.stateCategory)
	print(io, ",")
    show( io, xe_info.unit)
    print(io, ",", xe_info.startOrInit)
    print(io, ",", xe_info.fixed)
    print(io, ",", xe_info.nominal)
    print(io, ",", xe_info.unbounded)
    print(io, ",", xe_info.fixedLength)
    print(io, ",", xe_info.scalar)
    print(io, ",", xe_info.length)
    print(io, ",", xe_info.startIndex)
    print(io, ",", xe_info.x_segmented_startIndex)
    print(io, ")")
    return nothing
end


"""
    x_table = get_x_table(x_info::Vector{StateElementInfo})

Return the state element info as DataFrames.DataFrame table, for example
to print it as:

```
show(x_table, allrows=true, allcols=true, summary=false, eltypes=false)
```
"""
function get_x_table(x_info::Vector{StateElementInfo})
    x_table = DataFrames.DataFrame(name=String[], unit=String[], startOrInit=[], fixed=Bool[], nominal=Float64[], unbounded=Bool[])

    for xe_info in x_info
        push!(x_table, (xe_info.x_name, xe_info.unit, xe_info.startOrInit, xe_info.fixed, xe_info.nominal, xe_info.unbounded))
    end

    return x_table
end


@enum EquationInfoStatus EquationInfo_Instantiated EquationInfo_Initialized_Before_All_States_Are_Known EquationInfo_After_All_States_Are_Known

"""
    eqInfo = EquationInfo()

Return initial instance `eqInfo` that defines the information for the model equations
(especially, states and linear equation systems).
"""
mutable struct EquationInfo
    status::EquationInfoStatus                     # Defines which information is consistent in EquationInfo
    ode::Bool                                      # = true  if ODE-interface (`getDerivatives!`),
                                                   # = false if DAE-interface (`getResiduals!`).
    nz::Int                                        # Number of crossing functions
    x_info::Vector{StateElementInfo}               # Vector of StateElementInfo elements provding info for every x-element
    residualCategories::Vector{ResidualCategory}   # If ode=true, length(residualCategories) = 0
                                                   # If ode=false, residualCategories[j] is the ResidualCategory of residual[j]
    linearEquations::Vector{Tuple{Vector{String},AbstractVector,Vector{Int},Int,Bool}}
                                                   # linearEquations[i] defines a
                                                   # Modia.LinearEquations system, where the first tuple value
                                                   # is a vector of the names of the unknowns,
                                                   # the second tuple value is a vector of the Julia names of the vector-valued elements,
                                                   # the third tuple value is a vector with the lengths of the unknowns,
                                                   # the fourth tuple value is the number of residuals and the fifth tuple value
                                                   # defines whether the coefficient matrix A
                                                   # has only constant entries (=true) or not (=false).
    vSolvedWithFixedTrue::Vector{String}           # Vector of variables that are computed
                                                   # from other variables and have `fixed=true`. During initialization
                                                   # it is checked whether the calculated values and the start values of
                                                   # these variables agree. If this is not the case, an error is triggered.
    nx::Int                                        # = length(x) or -1 if not yet known
                                                   # This variable is updated once all states are known.
    nxInvariant::Int                               # = number of invariant x-elements (so x[1:nxInvariant] are invariant states) or -1 if not yet known
                                                   # This variable is updated once all states are known.
    nxSegmented::Int                               # = number of segmented x-elements (x[nxInvariant+1:nxInvariant+nxSegmented]).
                                                   # This variable is always updated consistently via function new_x_segmented_variable!(..)
                                                   # (nxSegmented=0, if there are no segmented states yet).
    nx_info_fixedLength::Int                       # x_info[1:nx_info_fixedLength] are states with fixed length (does not change after compilation) or -1 if not yet known
    nx_info_invariant::Int                         # x_info[1:nx_info_invariant] are states that are visible in getDerivatives!(..) or -1 if not yet known
                                                   # x_info[nx_info_invariant+1:end] are states defined in functions that are not visible in getDerivatives!(..)
    #x_infoByIndex::Vector{Int}                    # i = x_infoByIndex[j] -> x_info[i]
    #                                              # or empty vector, if not yet known.
    x_dict::OrderedCollections.OrderedDict{String,Int}           # x_dict[x_name] returns the index of x_name with respect to x_info
    der_x_dict::OrderedCollections.OrderedDict{String,Int}       # der_x_dict[der_x_name]  returns the index of der_x_name with respect to x_info
    defaultParameterAndStartValues::Union{AbstractDict,Nothing}  # Dictionary of default parameter and default start values.

    function EquationInfo()
        ode                   = true
        nz                    = 0
        x_info                = StateElementInfo[]
        residualCategories    = ResidualCategory[]
        linearEquations       = Tuple{Vector{String},AbstractVector,Vector{Int},Int,Bool}[]
        vSolvedWithFixedTrue  = String[]
        nx                    = -1
        nxInvariant           = -1
        nxSegmented           =  0
        nx_info_fixedLength   = -1
        nx_info_invariant     = -1
        x_dict                = OrderedCollections.OrderedDict{String,Int}()
        der_x_dict            = OrderedCollections.OrderedDict{String,Int}()
        defaultParameterAndStartValues = nothing
        new(EquationInfo_Instantiated, ode, nz, x_info, residualCategories, linearEquations, vSolvedWithFixedTrue,
            nx, nxInvariant, nxSegmented, nx_info_fixedLength, nx_info_invariant, x_dict, der_x_dict,
            defaultParameterAndStartValues)
    end
end


"""
    initEquationInfo!(eqInfo::EquationInfo, nx_info_fixedLength)

Initialize `eqInfo` before code generation.
`nx_info_fixedLength` are the number of states with fixed length (do not change after compilation).

The following variables get a consistent value:

```
eqInfo.x_dict[x_name]         = ...,
eqInfo.der_x_dict[der_x_name] = ...,
eqInfo.nx, eqInfo.nxInvariant,
eqInfo.nx_info_fixedLength, eqInfo.nx_info_invariant,
eqInfo.x_info[:].startIndex.
```

This function must be called before segmented states are set (eqInfo.nxSegmented=0).
"""
function initEquationInfo!(eqInfo::EquationInfo, nx_info_fixedLength::Int)::Nothing
    @assert(eqInfo.status == EquationInfo_Instantiated)
    @assert(eqInfo.nxSegmented == 0)
    x_dict     = eqInfo.x_dict
    der_x_dict = eqInfo.der_x_dict
    startIndex = 1
    for (i, xi_info) in enumerate(eqInfo.x_info)
        x_dict[xi_info.x_name] = i
        der_x_dict[xi_info.der_x_name] = i
        xi_info.startIndex = startIndex
        startIndex += xi_info.length
    end
    eqInfo.nx           = startIndex - 1
    eqInfo.nxInvariant  = eqInfo.nx
    eqInfo.nx_info_fixedLength = nx_info_fixedLength
    eqInfo.nx_info_invariant = length(eqInfo.x_info)

    eqInfo.status = EquationInfo_Initialized_Before_All_States_Are_Known
    return nothing
end


"""
    removeSegmentedStates!(eqInfo::EquationInfo)

Remove all segmented states from `eqInfo`.
"""
function removeSegmentedStates!(eqInfo::EquationInfo)::Nothing
    @assert(eqInfo.status == EquationInfo_After_All_States_Are_Known)
    if eqInfo.nx > eqInfo.nxInvariant
        for i = eqInfo.nx_info_invariant+1:length(eqInfo.x_info)
            xi_info = eqInfo.x_info[i]
            delete!(eqInfo.x_dict    , xi_info.x_name)
            delete!(eqInfo.der_x_dict, xi_info.der_x_name)
        end
        resize!(eqInfo.x_info, eqInfo.nx_info_invariant)
        eqInfo.nx = eqInfo.nxInvariant
    end
    eqInfo.nxSegmented = 0
    eqInfo.status = EquationInfo_Initialized_Before_All_States_Are_Known
    return nothing
end


"""
    x_start = initialStateVector!(eqInfo::EquationInfo, FloatType, isFirstSegment, x_terminate)::Vector{FloatType}

The function updates `eqInfo` (e.g. sets eqInfo.nx, eqInfo.nxInvariant) and returns the initial state vector x_start.

This function must be called, after all states are known (after calling propagateEvaluateAndInstantiate!(..)).
"""
function initialStateVector!(eqInfo::EquationInfo, FloatType::Type, isFirstSegment::Bool, x_terminate)::Vector{FloatType}
    @assert(eqInfo.status == EquationInfo_Initialized_Before_All_States_Are_Known)
    nx_info_fixedLength = eqInfo.nx_info_fixedLength
    x_info = eqInfo.x_info

    if nx_info_fixedLength == 0
        startIndex = 1
    else
        xi_info = x_info[nx_info_fixedLength]
        startIndex = xi_info.startIndex + xi_info.length
    end

    # If startOrInit is not defined, use a default value of zero.
    for xi_info in eqInfo.x_info
        if isnothing(xi_info.startOrInit)
            @info "State $(xi_info.x_name) has no start or init value defined. Using start value = 0.0."
            xi_info.startOrInit = FloatType(0)
            xi_info.scalar = true
        end
    end

    # Set startIndex for invariant states where the size was not fixed before code generation
    for i = nx_info_fixedLength+1:eqInfo.nx_info_invariant
        xi_info = x_info[i]
        xi_info.length     = length(xi_info.startOrInit)
        xi_info.startIndex = startIndex
        startIndex        += xi_info.length
    end
    eqInfo.nxInvariant = startIndex - 1

    # Set startIndex for segmented states
    for i = eqInfo.nx_info_invariant+1:length(x_info)
        xi_info = x_info[i]
        xi_info.length     = length(xi_info.startOrInit)
        xi_info.startIndex = startIndex
        startIndex        += xi_info.length
    end
    eqInfo.nx = startIndex - 1
    @assert(eqInfo.nx == eqInfo.nxInvariant + eqInfo.nxSegmented)

    # Construct x_start
    x_start = zeros(FloatType, eqInfo.nx)
    if isFirstSegment
        startIndex = 1
        for xe_info in x_info
            if xe_info.scalar
                @assert(length(xe_info.startOrInit) == 1)
                x_start[startIndex] = FloatType(ustrip(xe_info.startOrInit))
                startIndex += 1
            else
                xe_start = Vector{FloatType}(ustrip(xe_info.startOrInit))
                @assert(length(xe_start) == xe_info.length)
                copyto!(x_start, startIndex, xe_start, 1, length(xe_start))
                startIndex += length(xe_start)
            end
        end
    else
        for i in 1:eqInfo.nxInvariant
            x_start[i] = x_terminate[i]
        end
        startIndex = eqInfo.nxInvariant+1
        for i = eqInfo.nx_info_invariant+1:length(x_info)
            xe_info = x_info[i]
            if xe_info.scalar
                @assert(length(xe_info.startOrInit) == 1)
                x_start[startIndex] = FloatType(ustrip(xe_info.startOrInit))
                startIndex += 1
            else
                xe_start = Vector{FloatType}(ustrip(xe_info.startOrInit))
                @assert(length(xe_start) == xe_info.length)
                copyto!(x_start, startIndex, xe_start, 1, length(xe_start))
                startIndex += length(xe_start)
            end
        end
    end

    @assert(eqInfo.nx == startIndex - 1)
    eqInfo.status = EquationInfo_After_All_States_Are_Known

    # Final check
    for (i, xi_info) = enumerate(eqInfo.x_info)
        @assert(xi_info.startIndex > 0)
        if i <= eqInfo.nx_info_invariant
            @assert(xi_info.x_segmented_startIndex == -1)
        else
            @assert(xi_info.x_segmented_startIndex > 0)
        end
    end
    return x_start
end


"""
    nvec = get_equationSizes(equationInfo)

Return the sizes of the linear equations as `nvec::Vector{Int}`.
(e.g. `nvec=[10,3]` means that there are two linear equation systems with
size 10 and size 3).
"""
get_equationSizes(eqInfo::EquationInfo) = Int[leq[3] for leq in eqInfo.linearEquations]


function Base.show(io::IO, eqInfo::EquationInfo; indent=4)
    indentation  = repeat(" ", indent)
    indentation2 = repeat(" ", 2*indent)
    indentation3 = repeat(" ", 3*indent)
    println(io, "Modia.EquationInfo(")
    println(io, indentation2, "ode = ", eqInfo.ode, ",")
    if eqInfo.nz > 0
        println(io, indentation2, "nz = ", eqInfo.nz, ",")
    end
    print(io, indentation2, "x_info = Modia.StateElementInfo[")
    for (i, xe_info) in enumerate(eqInfo.x_info)
        if i == 1
            print(io, "\n", indentation3)
        else
            print(io, ",\n", indentation3)
        end
        show(io, xe_info)
    end
    print(io,"]")

    if length(eqInfo.residualCategories) > 0
        print(io, ",\n", indentation2, "residualCategories = [")
        for (i, rcat) in enumerate(eqInfo.residualCategories)
            if i > 1
                print(io, ",")
            end
            show(io, rcat)
        end
        print(io,"]")
    end

    leqs = eqInfo.linearEquations
    if length(leqs) > 0
        println(io, ",\n", indentation2, "linearEquations = [")
        for (i, leq) in enumerate(leqs)
            print(io, indentation3, "(", leq[1], ",\n",
                      indentation3, " ", leq[2], ",\n",
                      indentation3, " ", leq[3], ", ", leq[4], ")")
            if i < length(leqs)
                print(io, ",\n")
            end
        end
        print(io, "]")
    end


    if length(eqInfo.vSolvedWithFixedTrue) > 0
        print(io, ",\n", indentation2, "vSolvedWithFixedTrue = ")
        show(io, eqInfo.vSolvedWithFixedTrue)
    end

    if eqInfo.nx >= 0
        nx = eqInfo.nx
        nxInvariant = eqInfo.nxInvariant
        nxSegmented  = eqInfo.nxSegmented
        nx_info_fixedLength = eqInfo.nx_info_fixedLength
        nx_info_invariant = eqInfo.nx_info_invariant
        print(io, ",\n", indentation2, "nx = $nx, nxInvariant = $nxInvariant, nxSegmented = $nxSegmented")
        print(io, ",\n", indentation2, "nx_info_fixedLength = $nx_info_fixedLength, nx_info_invariant = $nx_info_invariant")
    end

    if !isnothing(eqInfo.defaultParameterAndStartValues)
        print(io, ",\n", indentation2, "defaultParameterAndStartValues = ")
        show(io, eqInfo.defaultParameterAndStartValues, indent=12, finalLineBreak=false)
    end

    println(io, "\n", indentation, ")")
    return nothing
end



"""
    names = get_stateNames(equationInfo::EquationInfo)

Return the names of the states defined in `equationInfo` as a Vector of strings.
"""
get_stateNames(eqInfo::EquationInfo) = String[xi_info.x_name for xi_info in eqInfo.x_info]


"""
    names = get_xNames(eqInfo::EquationInfo)

Return the names of all elements of the x-vector as a vector of strings.
"""
function get_xNames(eqInfo::EquationInfo)::Vector{String}
    xNames = Vector{String}(undef, eqInfo.nx)
    for xe_info in eqInfo.x_info
        if xe_info.length == 1
            xNames[xe_info.startIndex] = xe_info.x_name
        else
            for i = 1:xe_info.length
                xNames[xe_info.startIndex+i-1] = xe_info.x_name*"["*string(i)*"]"
            end
        end
    end
    return xNames
end

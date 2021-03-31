# License for this file: MIT (expat)
# Copyright 2020-2021, DLR Institute of System Dynamics and Control


using  ModiaBase
using  Unitful
using  Measurements
import MonteCarloMeasurements
using  DataStructures: OrderedDict, OrderedSet
using  DataFrames

export SimulationModel, measurementToString, get_lastValue


"""
    baseType(T)

Return the base type of a type T.

# Examples
```
baseType(Float32)                # Float32
baseType(Measurement{Float64})   # Float64
```
"""
baseType(::Type{T})                                           where {T}                  = T
baseType(::Type{Measurements.Measurement{T}})                 where {T<:AbstractFloat}   = T
baseType(::Type{MonteCarloMeasurements.Particles{T,N}})       where {T<:AbstractFloat,N} = T
baseType(::Type{MonteCarloMeasurements.StaticParticles{T,N}}) where {T<:AbstractFloat,N} = T


"""
    str = measurementToString(v)

Return variable `v::Measurements.Measurement{FloatType}` or a vector of such variables
in form of a string will the full number of significant digits.
"""
measurementToString(v::Measurements.Measurement{FloatType}) where {FloatType} =
       string(Measurements.value(v)) * " ± " * string(Measurements.uncertainty(v))

function measurementToString(v::Vector{Measurements.Measurement{FloatType}}) where {FloatType}
    str = string(typeof(v[1])) * "["
    for (i,vi) in enumerate(v)
        if i > 1
            str = str * ", "
        end
        str = str * measurementToString(vi)
    end
    str = str * "]"
    return str
end



"""
    simulationModel = SimulationModel{FloatType}(
            modelModule, modelName, getDerivatives!, equationInfo, x_startValues,
            parameters, variableNames;
            modelModule = nothing,
            vSolvedWithInitValuesAndUnit::OrderedDict{String,Any}(),
            vEliminated::Vector{Int}=Int[],
            vProperty::Vector{Int}=Int[],
            var_name::Function = v->nothing)


# Arguments

- `modelModule`: Module in which `@instantiateModel` is invoked (it is used for `Core.eval(modelModule, ...)`),
  that is evaluation of expressions in the environment of the user.
- `modelName::String`: Name of the model
- `getDerivatives::Function`: Function that is used to evaluate the model equations,
  typically generated with [`TinyModia.generate_getDerivatives!`].
- `equationInfo::ModiaBase.EquationInfo`: Information about the states and the equations.
- `x_startValues`:: Deprecated (is no longer used).
- `parameters`: A hierarchical NamedTuple of (key, value) pairs defining the parameter and init/start values.
- variableNames: A vector of variable names. A name can be a Symbol or a String.
"""
mutable struct SimulationModel{FloatType}
    modelModule
    modelName::String
    getDerivatives!::Function
    equationInfo::ModiaBase.EquationInfo
    linearEquations::Vector{ModiaBase.LinearEquations{FloatType}}
    variables::OrderedDict{String,Int}      # Dictionary of variables and their result indices (negated alias has negativ index)
    zeroVariables::OrderedSet{String}       # Set of variables that are identically to zero
    vSolvedWithInitValuesAndUnit::OrderedDict{String,Any}   # Dictionary of (names, init values with units) for all explicitly solved variables with init-values defined
    p::AbstractVector                       # Parameter and init/start values 
                                            # (parameter values are used in getDerivatives!,
                                            #  init/start values are extracted and stored in x_start)                                     
    separateObjects::OrderedDict{Int,Any}   # Dictionary of separate objects    
    storeResult::Bool
    isInitial::Bool
    isTerminal::Bool
    time::FloatType
    nGetDerivatives::Int                    # Number of getDerivatives! calls
    x_start::Vector{FloatType}
    der_x::Vector{FloatType}
    result::Vector{Tuple}
    algorithmType::DataType                  # Type of integration algorithm (used in default-heading of plot)

    function SimulationModel{FloatType}(modelModule, modelName, getDerivatives!, equationInfo, x_startValues,
                                        parameters, variableNames;
                                        vSolvedWithInitValuesAndUnit::AbstractDict = OrderedDict{String,Any}(),
                                        vEliminated::Vector{Int} = Int[],
                                        vProperty::Vector{Int}   = Int[],
                                        var_name::Function       = v -> nothing) where {FloatType}

        # Construct result dictionaries
        variables = OrderedDict{String,Int}()
        zeroVariables = OrderedSet{String}()
        vSolvedWithInitValuesAndUnit2 = OrderedDict{String,Any}( [(string(key),vSolvedWithInitValuesAndUnit[key]) for key in keys(vSolvedWithInitValuesAndUnit)] )
            # Store variables
            for (i, name) in enumerate(variableNames)
                variables[string(name)] = i
            end

            # Store eliminated variables
            for v in vEliminated
                name = var_name(v)
                if ModiaBase.isZero(vProperty, v)
                    push!(zeroVariables, name)
                elseif ModiaBase.isAlias(vProperty, v)
                    name2 = var_name( ModiaBase.alias(vProperty, v) )
                    variables[name] = variables[name2]
                else # negated alias
                    name2 = var_name( ModiaBase.negAlias(vProperty, v) )
                    variables[name] = -variables[name2]
                end
            end

        # Construct parameter values that are copied into the code
        parameterValues = [eval(p) for p in values(parameters)]

        # Construct data structure for linear equations
        linearEquations = ModiaBase.LinearEquations{FloatType}[]
        for leq in equationInfo.linearEquations
            push!(linearEquations, ModiaBase.LinearEquations{FloatType}(leq...))
        end

        # Construct dictionary for separate objects
        separateObjects = OrderedDict{Int,Any}()
                
        # Initialize execution flags
        isInitial       = true
        isTerminal      = false
        storeResult     = false
        nGetDerivatives = 0


        new(modelModule, modelName, getDerivatives!, equationInfo, linearEquations, variables, zeroVariables,
            vSolvedWithInitValuesAndUnit2, parameterValues,
            separateObjects, storeResult, isInitial, isTerminal, convert(FloatType, 0), nGetDerivatives, 
            zeros(FloatType,0), zeros(FloatType,0), Tuple[])
    end
end


"""
    floatType = getFloatType(simulationModel::SimulationModel)

Return the floating point type with which `simulationModel` is parameterized
(for example returns: `Float64, Float32, DoubleFloat, Measurements.Measurement{Float64}`).
"""
getFloatType(m::SimulationModel{FloatType}) where {FloatType} = FloatType


"""
    hasParticles(value)
    
Return true, if `value` is of type `MonteCarloMeasurements.StaticParticles` or
`MonteCarloMeasurements.Particles`.
"""
hasParticles(value) = typeof(value) <: MonteCarloMeasurements.StaticParticles ||
                      typeof(value) <: MonteCarloMeasurements.Particles
                          
                          
"""
    get_value(obj::NamedTuple, name::String)
    
Return the value identified by `name` from the potentially hierarchically 
`NamedTuple obj`. If `name` is not in `obj`, the function returns `missing`.

# Examples
```julia
s1 = (a = 1, b = 2, c = 3)
s2 = (v1 = s1, v2 = (d = 4, e = 5))
s3 = (v3 = s2, v4 = s1)

@show get_value(s3, "v3.v1.b")   # returns 2
@show get_value(s3, "v3.v2.e")   # returns 5
@show get_value(s3, "v3.v1.e")   # returns missing
```
"""
function get_value(obj::NamedTuple, name::String)
    if length(name) == 0 || length(obj) == 0
        return missing
    end
    j = findnext('.', name, 1)
    if isnothing(j)
        key = Symbol(name)
        return haskey(obj,key) ? obj[key] : missing
    elseif j == 1
        return missing
    else
        key = Symbol(name[1:j-1])
        if haskey(obj,key) && typeof(obj[key]) <: NamedTuple && length(name) > j
            get_value(obj[key], name[j+1:end])
        else
            return missing
        end
    end
end


"""
    appendName(path::String, name::Symbol)
   
Return `path` appended with `.` and `string(name)`.
"""
appendName(path, key) = path == "" ? string(key) : path * "." * string(key)


"""
    get_names(obj::NamedTuple)
    
Return `Vector{String}` containing all the names present in `obj`

# Examples
```julia
s1 = (a = 1, b = 2, c = 3)
s2 = (v1 = s1, v2 = (d = 4, e = 5))
s3 = (v3 = s2, v4 = s1)

@show get_names(s3)
```
"""
function get_names(obj::NamedTuple)
    names = String[]
    get_names!(obj, names, "")
    return names
end
function get_names!(obj::NamedTuple, names::Vector{String}, path::String)::Nothing
    for (key,value) in zip(keys(obj), obj)
        name = appendName(path, key)
        if typeof(value) <: NamedTuple
            get_names!(value, names, name)
        else
            push!(names, name)
        end
    end
    return nothing
end


"""
    get_lastValue(model::SimulationModel, name::String; unit=true)

Return the last stored value of variable `name` from `model`.
If `unit=true` return the value with its unit, otherwise with stripped unit.

If `name` is not known or no result values yet available, an info message is printed
and the function returns `nothing`.
"""
function get_lastValue(m::SimulationModel, name::String; unit=true)
    if haskey(m.variables, name)
        if length(m.result) == 0
            @info "get_lastValue(model,\"$name\"): No results yet available."
            return nothing
        end

        resIndex = m.variables[name]
        negAlias = false
        if resIndex < 0
            resIndex = -resIndex
            negAlias = true
        end
        value = m.result[end][resIndex]
        if negAlias
            value = -value
        end
    elseif name in m.zeroVariables
        # Unit missing (needs to be fixed)
        value = 0.0
    else
        value = get_value(m.p[1], name)
        if ismissing(value)
            @info "get_lastValue: $name is not known and is ignored."
            return nothing;
        end
    end

    return unit ? value : ustrip(value)
end


get_xe(x, xe_info) = xe_info.length == 1 ? x[xe_info.startIndex] : x[xe_info.startIndex:xe_info.startIndex + xe_info.length-1]

function set_xe!(x, xe_info, value)::Nothing
    if xe_info.length == 1 
        x[xe_info.startIndex] = value
    else
        x[xe_info.startIndex:xe_info.startIndex + xe_info.length-1] = value
    end
    return nothing
end


"""
    success = init!(simulationModel, startTime, tolerance, merge, log, logParameters, logStates)


Initialize `simulationModel::SimulationModel` at `startTime`. In particular:

- Empty result data structure.

- Merge parameter and init/start values into simulationModel.

- Construct x_start.

- Call simulationModel.getDerivatives! once with isInitial = true to
  compute and store all variables in the result data structure at `startTime`
  and initialize simulationModel.linearEquations.

- Check whether explicitly solved variables that have init-values defined,
  have the required value after initialization (-> otherwise error).
  
If initialization is successful return true, otherwise false.
"""
function init!(m::SimulationModel, startTime, tolerance, merge, 
               log::Bool, logParameters::Bool, logStates::Bool)::Bool
    empty!(m.result)

	# Apply updates from merge Map and log parameters
    if isnothing(merge)
        if logParameters
            parameters = m.p[1]
            @showModel parameters
        end
    else
        m.p = [recursiveMerge(m.p[1], merge)]
        if logParameters
            updatedParameters = m.p[1]
            @showModel updatedParameters
        end
    end
	
    # Re-initialize dictionary of separate objects
    empty!(m.separateObjects)

    # Update startIndex, length in x_info and determine x_start
    FloatType = getFloatType(m)
    x_start = []
    startIndex = 1
    for xe_info in m.equationInfo.x_info
        xe_info.startIndex = startIndex
        
        # Determine init/start value in m.p[1]
        xe_value = get_value(m.p[1], xe_info.x_name)
        if ismissing(xe_value)
            error("Missing start/init value ", xe_info.x_name)
        end
        if hasParticles(xe_value)
            len = 1
            push!(x_start, xe_value)
        else
            len = length(xe_value)
            push!(x_start, xe_value...)
        end 
        if len != xe_info.length
            printstyled("Model error: ", bold=true, color=:red)  
            printstyled("Length of ", xe_info.x_name, " shall be changed from ",
                        xe_info.length, " to $len\n",
                        "This is currently not support in TinyModia.", bold=true, color=:red)
            return false
        end        
        startIndex += xe_info.length
    end
    nx = startIndex - 1
    m.equationInfo.nx = nx
    
    # Temporarily remove units from x_start
    # (TODO: should first transform to the var_unit units and then remove)    
    converted_x_start = convert(Vector{FloatType}, [ustrip(v) for v in x_start])  # ustrip.(x_start) does not work for MonteCarloMeasurements
    m.x_start = deepcopy(converted_x_start)

    if logStates
        # List init/start values
        x_table = DataFrames.DataFrame(state=String[], init=Any[], unit=String[], nominal=Float64[])   
        for xe_info in m.equationInfo.x_info
            xe_init = get_xe(m.x_start, xe_info)
            if hasParticles(xe_init)
                xe_init = string(minimum(xe_init)) * " .. " * string(maximum(xe_init))
            end
            push!(x_table, (xe_info.x_name, xe_init, xe_info.unit, xe_info.nominal))
        end
        show(stdout, x_table; allrows=true, allcols=true, rowlabel = Symbol("#"), summary=false, eltypes=false)
        println("\n")
    end
    
    m.der_x = zeros(FloatType, nx)
        
    # Initialize model, linearEquations and compute and store all variables at the initial time
    if log
        println("      Initialization at time = ", startTime, " s")  
    end        
    m.nGetDerivatives = 0
    m.isInitial   = true
    m.storeResult = true
#    m.getDerivatives!(m.der_x, m.x_start, m, startTime)
    Base.invokelatest(m.getDerivatives!, m.der_x, m.x_start, m, startTime)
    m.isInitial   = false
    m.storeResult = false
    # Check vSolvedWithInitValuesAndUnit
    if length(m.vSolvedWithInitValuesAndUnit) > 0
        names = String[]
        valuesBeforeInit = Any[]
        valuesAfterInit  = Any[]

        for (name, valueBefore) in m.vSolvedWithInitValuesAndUnit
            valueAfterInit  = get_lastValue(m, name, unit=false)
            valueBeforeInit = ustrip.(valueBefore)
            if !isnothing(valueAfterInit) && abs(valueBeforeInit - valueAfterInit) >= max(abs(valueBeforeInit),abs(valueAfterInit),0.01*tolerance)*tolerance
                push!(names, name)
                push!(valuesBeforeInit, valueBeforeInit)
                push!(valuesAfterInit , valueAfterInit)
            end
        end

        if length(names) > 0
            v_table = DataFrames.DataFrame(name=names, beforeInit=valuesBeforeInit, afterInit=valuesAfterInit)
            #show(stderr, v_table; allrows=true, allcols=true, summary=false, eltypes=false)
            #print("\n\n")
            ioTemp = IOBuffer();
            show(ioTemp, v_table; allrows=true, allcols=true, rowlabel = Symbol("#"), summary=false, eltypes=false)
            str = String(take!(ioTemp))
            close(ioTemp)
            printstyled("Model error: ", bold=true, color=:red)  
            printstyled("The following variables are explicitly solved for, have init-values defined\n",
                        "and after initialization the init-values are not respected\n",
                        "(remove the init-values in the model or change them to start-values):\n",
                        str, bold=true, color=:red)
            return false
        end
    end
    return true
end


"""
    terminate!(m::SimulationModel, x, time)

Terminate model.
"""
function terminate!(m::SimulationModel, x, t)::Nothing
    m.isTerminal = true
    Base.invokelatest(m.getDerivatives!, m.der_x, x, m, t)
    m.isTerminal = false
    return nothing
end
   
   
"""
    outputs!(x, t, integrator)

DifferentialEquations FunctionCallingCallback function for `SimulationModel`
that is used to store results at communication points.
"""
function outputs!(x, t, integrator)::Nothing
    m = integrator.p
    m.storeResult = true
#    m.getDerivatives!(m.der_x, x, m, t)
    Base.invokelatest(m.getDerivatives!, m.der_x, x, m, t)

    m.storeResult = false
    return nothing
end


"""
    derivatives!(derx, x, m, t)

DifferentialEquations callback function to get the derivatives.
"""
function derivatives!(der_x, x, m, t)::Nothing
#    m.getDerivatives!(der_x, x, m, t)
    Base.invokelatest(m.getDerivatives!, der_x, x, m, t)
    return nothing
end


   
"""
    addToResult!(simulationModel, variableValues...)

Add `variableValues...` to `simulationModel::SimulationModel`.
It is assumed that the first variable in `variableValues` is `time`.
"""
function addToResult!(m::SimulationModel, variableValues...)::Nothing
    push!(m.result, variableValues)
    return nothing
end




"""
    code = generate_getDerivatives!(AST, equationInfo, parameters, variables, functionName;
                                    hasUnits=false)

Return the code of the `getDerivatives!` function as `Expr` using the
Symbol `functionName` as function name. By `eval(code)` or
`fc = @RuntimeGeneratedFunction(code)` the function is compiled and can afterwards be called.

# Arguments

- `AST::Vector{Expr}`: Abstract Syntax Tree of the equations as vector of `Expr`.

- `equationInfo::ModiaBase.EquationInfo`: Data structure returned by `ModiaBase.getSortedAndSolvedAST
            holding information about the states.

- `parameters`: Vector of parameter names (as vector of symbols)

- `variables`: Vector of variable names (as vector of symbols). The first entry is expected to be time, so `variables[1] = :time`.

- `functionName::Function`: The name of the function that shall be generated.


# Optional Arguments

- `hasUnits::Bool`: = true, if variables have units. Note, the units of the state vector are defined in equationinfo.
"""
function generate_getDerivatives!(AST::Vector{Expr}, equationInfo::ModiaBase.EquationInfo,
                                  parameters, variables, functionName::Symbol;
                                  hasUnits=false)

    # Generate code to copy x to struct and struct to der_x
    x_info     = equationInfo.x_info
    code_x     = Expr[]
    code_der_x = Expr[]
    code_p     = Expr[]

    if length(x_info) == 1 && x_info[1].x_name == "" && x_info[1].der_x_name == ""
        # Explicitly solved pure algebraic variables. Introduce dummy equation
        push!(code_der_x, :( _der_x[1] = -_x[1] ))
    else
        i1 = 0
        i2 = 0
        for (i, xe) in enumerate(x_info)
            i1 = i2 + 1
            i2 = i1 + xe.length - 1
            indexRange = i1 == i2 ? :($i1) :  :( $i1:$i2 )
            x_name     = xe.x_name_julia
            der_x_name = xe.der_x_name_julia
            # x_name     = Meta.parse("m."*xe.x_name)
            # der_x_name = Meta.parse("m."*replace(xe.der_x_name, r"der\(.*\)" => s"var\"\g<0>\""))
            if !hasUnits || xe.unit == ""
                push!(code_x, :( $x_name = _x[$indexRange] ) )
            else
                x_unit = xe.unit
                push!(code_x, :( $x_name = _x[$indexRange]*@u_str($x_unit)) )
            end
            if hasUnits
                push!(code_der_x, :( _der_x[$indexRange] = ustrip( $der_x_name )) )
            else
                push!(code_der_x, :( _der_x[$indexRange] = $der_x_name ))
            end
        end
    end
    for (i,pi) in enumerate(parameters)
        push!(code_p, :( $pi = _m.p[$i] ) )
    end

    timeName = variables[1]
    if hasUnits
        code_time = :( $timeName = _time*u"s" )
    else
        code_time = :( $timeName = _time )
    end

    # Generate code of the function
    code = quote
                function $functionName(_der_x, _x, _m, _time)::Nothing
                    _m.time = _time
                    _m.nGetDerivatives += 1
                    simulationModel = _m
                    $code_time
                    $(code_p...)
                    $(code_x...)
                    $(AST...)
                    $(code_der_x...)

                    if _m.storeResult
                        TinyModia.addToResult!(_m, $(variables...))
                    end

                    return nothing
                end
            end
    return code
end

# License for this file: MIT (expat)
# Copyright 2020-2021, DLR Institute of System Dynamics and Control


using  ModiaBase
using  Unitful
using  Measurements
import MonteCarloMeasurements
using  DataStructures: OrderedDict
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
            modelName, getDerivatives!, equationInfo, x_startValues,
            parameters, variableNames;
            v_init::OrderedDict{String,Any}(),
            vEliminated::Vector{Int}=Int[],
            vProperty::Vector{Int}=Int[],
            var_name::Function = v->nothing)


# Arguments

- parameters: A dictionary of (key, value) pairs. A key can be a Symbol or a String.
- variableNames: A vector of variable names. A name can be a Symbol or a String.
"""
mutable struct SimulationModel{FloatType}
    modelName::String
    getDerivatives!::Function
    equationInfo::ModiaBase.EquationInfo
    linearEquations::Vector{ModiaBase.LinearEquations{FloatType}}
    variables::OrderedDict{String,Int}                       # Dictionary of variables and their result indices (negated alias has negativ index)
    parametersAndConstantVariables::OrderedDict{String,Any}  # Dictionary of parameters and constant variables with their values
    vSolvedWithInitValuesAndUnit::OrderedDict{String,Any}     # Dictionary of (names, init values with units) for all explicitly solved variables with init-values defined
    p::AbstractVector               # Parameter values that are copied into the code
    storeResult::Bool
    isInitial::Bool
    nGetDerivatives::Int            # Number of getDerivatives! calls
    x_start::Vector{FloatType}
    der_x::Vector{FloatType}
    result::Vector{Tuple}
    algorithmType::DataType         # Type of integration algorithm (used in default-heading of plot)

    function SimulationModel{FloatType}(modelName, getDerivatives!, equationInfo, x_startValues,
                                        parameters, variableNames;
                                        vSolvedWithInitValuesAndUnit::AbstractDict = OrderedDict{String,Any}(),
                                        vEliminated::Vector{Int} = Int[],
                                        vProperty::Vector{Int}   = Int[],
                                        var_name::Function       = v -> nothing) where {FloatType}

        # Construct result dictionaries
        variables = OrderedDict{String,Int}()
        parametersAndConstantVariables = OrderedDict{String,Any}( [(string(key),parameters[key]) for key in keys(parameters)] )
        vSolvedWithInitValuesAndUnit2 = OrderedDict{String,Any}( [(string(key),vSolvedWithInitValuesAndUnit[key]) for key in keys(vSolvedWithInitValuesAndUnit)] )
            # Store variables
            for (i, name) in enumerate(variableNames)
                variables[string(name)] = i
            end

            # Store eliminated variables
            for v in vEliminated
                name = var_name(v)
                if ModiaBase.isZero(vProperty, v)
                    parametersAndConstantVariables[name] = 0.0
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

        # Set startIndex in x_info and compute nx
        startIndex = 1
        for xi_info in equationInfo.x_info
            xi_info.startIndex = startIndex
            startIndex += xi_info.length
        end
        nx = startIndex - 1
        equationInfo.nx = nx
        @assert(nx == length(x_startValues))
        x_start = deepcopy(x_startValues)

        isInitial       = true
        storeResult     = false
        nGetDerivatives = 0

        new(modelName, getDerivatives!, equationInfo, linearEquations, variables, parametersAndConstantVariables,
            vSolvedWithInitValuesAndUnit2, parameterValues,
            storeResult, isInitial, nGetDerivatives, x_start, zeros(FloatType,nx), Tuple[])
    end
end


"""
    floatType = getFloatType(simulationModel::SimulationModel)

Return the floating point type with which `simulationModel` is parameterized
(for example returns: `Float64, Float32, DoubleFloat, Measurements.Measurement{Float64}`).
"""
getFloatType(m::SimulationModel{FloatType}) where {FloatType} = FloatType


function get_states!(equationInfo::ModiaBase.EquationInfo,
                     states::Vector{FloatType},
                     x_names::Vector{Symbol},
                     x_startValues::AbstractVector)::Nothing where {FloatType}
    states .= 0
    x_info = equationInfo.x_info

 #=
    for xi_info in equationInfo.x_info
        if xi_info.x_name != ""
            (component, key) = ModiaBase.get_modelValuesAndName(modelValues, xi_info.x_name)
            if !isdefined(component,key)
                error("From ModiaBase:\nState ", key, " in component ", typeof(component), " has no value.")
            end
            value = ustrip( getfield(component,key) )
            istart = xi_info.startIndex
            for i = 1:xi_info.length
                states[istart+i-1] = value[i]
            end
        end
    end
=#
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
    elseif haskey(m.parametersAndConstantVariables, name)
        value = eval( m.parametersAndConstantVariables[name] )
    else
        @info "get_lastValue: $name is not known and is ignored."
        return nothing;
    end

    return unit ? value : ustrip(value)
end



"""
    init!(simulationModel, startTime, tolerance)


Initialize `simulationModel::SimulationModel` at `startTime`. In particular:

- Empty result data structure
- Call simulationModel.getDerivatives! once with isInitial = true to
  compute and store all variables in the result data structure at `startTime`
  and initialize simulationModel.linearEquations.
- Check whether explicitly solved variables that have init-values defined,
  have the required value after initialization (-> otherwise error).
"""
function init!(m::SimulationModel, startTime, tolerance)::Nothing
    empty!(m.result)

    # Call getDerivatives once to compute and store all variables and initialize linearEquations
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

        for (name, valueBeforeInit) in m.vSolvedWithInitValuesAndUnit
            valueAfterInit = get_lastValue(m, name)
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
            error("The following variables are explicitly solved for, have init-values defined\n",
                  "and after initialization the init-values are not respected\n",
                  "(remove the init-values in the model or change them to start-values):\n",
                  str)
        end
    end

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

                    _m.nGetDerivatives += 1
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

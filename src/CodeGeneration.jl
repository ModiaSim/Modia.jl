# License for this file: MIT (expat)
# Copyright 2020-2021, DLR Institute of System Dynamics and Control


using  ModiaBase
using  Unitful
using  Measurements
import MonteCarloMeasurements
using  DataStructures: OrderedDict, OrderedSet
using  DataFrames

export SimulationModel, measurementToString, get_lastValue
export positive, negative, previous, edge, after, reinit, pre
export initial, terminal, isInitial, isTerminal
export get_xNames
export registerExtraSimulateKeywordArguments
export get_extraSimulateKeywordArgumentsDict


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

Base.floatmax(::Type{MonteCarloMeasurements.Particles{T,N}})       where {T<:AbstractFloat,N} = Base.floatmax(T)
Base.floatmax(::Type{MonteCarloMeasurements.StaticParticles{T,N}}) where {T<:AbstractFloat,N} = Base.floatmax(T)


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


#= No longer needed
function get_x_start!(FloatType, equationInfo, parameters)
    x_start = []
    startIndex = 1
    for xe_info in equationInfo.x_info
        xe_info.startIndex = startIndex
        
        # Determine init/start value in m.parameters
        xe_value = get_value(parameters, xe_info.x_name)
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
    equationInfo.nx = startIndex - 1
    
    # Temporarily remove units from x_start
    # (TODO: should first transform to the var_unit units and then remove)    
    converted_x_start = convert(Vector{FloatType}, [stripUnit(v) for v in x_start])  # stripUnit.(x_start) does not work for MonteCarloMeasurements
    x_start2 = deepcopy(converted_x_start)
    return x_start2
end
=#

const BasicSimulationKeywordArguments = OrderedSet{Symbol}(
        [:merge, :tolerance, :startTime, :stopTime, :interval, :interp_points, :adaptive, :log, :logStates, :logEvents, :logParameters, :logEvaluatedParameters, :requiredFinalStates])
const RegisteredExtraSimulateKeywordArguments = OrderedSet{Symbol}()

function registerExtraSimulateKeywordArguments(keys)::Nothing
    for key in keys
        if key in BasicSimulationKeywordArguments
            @error "registerExtraSimulateKeywordArguments:\nKeyword $key is ignored, since standard keyword of simulate!(..) (use another keyword)\n\n"
        end
        push!(RegisteredExtraSimulateKeywordArguments, key)
    end
    return nothing
end


"""
    convertTimeVariable(TimeType, t)
    
If `t` has a unit, it is transformed to u"s", the unit is stripped off,
converted to `TimeType` and returned. Otherwise `t` is converted to `TimeType` and returned.

# Example
```julia
convertTimeVariable(Float32, 2.0u"hr")  # = 7200.0f0
```
"""
convertTimeVariable(TimeType, t) = typeof(t) <: Unitful.AbstractQuantity ? convert(TimeType, stripUnit(t)) : convert(TimeType, t)


struct SimulationOptions{FloatType,TimeType}
    merge::NamedTuple
    tolerance::Float64
    startTime::TimeType   # u"s"
    stopTime::TimeType    # u"s"
    interval::TimeType    # u"s"
    desiredResultTimeUnit    
    interp_points::Int
    adaptive::Bool
    log::Bool
    logStates::Bool
    logEvents::Bool
    logParameters::Bool
    logEvaluatedParameters::Bool
    requiredFinalStates::Union{Nothing, Vector{FloatType}}
    extra_kwargs::OrderedDict{Symbol,Any}
    
    function SimulationOptions{FloatType,TimeType}(merge, errorMessagePrefix=""; kwargs...) where {FloatType,TimeType}  
        success   = true       
        #merge     = get(kwargs, :merge, NamedTuple())
        tolerance = get(kwargs, :tolerance, 1e-6)
        if tolerance <= 0.0
            printstyled(errorMessagePrefix, "tolerance (= $(tolerance)) must be > 0\n\n", bold=true, color=:red)
            success = false 
        end
        startTime = convertTimeVariable(TimeType, get(kwargs, :startTime, 0.0) )
        rawStopTime       = get(kwargs, :stopTime, 1.0)
        stopTime  = convertTimeVariable(TimeType, rawStopTime)
        interval  = convertTimeVariable(TimeType, get(kwargs, :interval , (stopTime - startTime)/500.0) )
        desiredResultTimeUnit = unit(rawStopTime)          
        interp_points = get(kwargs, :interp_points, 0)
        if interp_points < 0
            printstyled(errorMessagePrefix, "interp_points (= $(interp_points)) must be > 0\n\n", bold=true, color=:red)
            success = false    
        elseif interp_points == 1
            # DifferentialEquations.jl crashes
            interp_points = 2
        end      
        adaptive      = get(kwargs, :adaptive     , true)
        log           = get(kwargs, :log          , false)
        logStates     = get(kwargs, :logStates    , false)
        logEvents     = get(kwargs, :logEvents    , false)
        logParameters = get(kwargs, :logParameters, false)
        logEvaluatedParameters = get(kwargs, :logEvaluatedParameters, false)
        requiredFinalStates    = get(kwargs, :requiredFinalStates, nothing)

        extra_kwargs = OrderedDict{Symbol,Any}()        
        for option in kwargs
            key = option.first
            if key in BasicSimulationKeywordArguments
                nothing;
            elseif key in RegisteredExtraSimulateKeywordArguments
                extra_kwargs[key] = option.second
            else
                printstyled(errorMessagePrefix, "$key is unknown (keyword argument: $key = $(option.second))\n\n", bold=true, color=:red)
                success = false
            end
        end    

        obj = new(isnothing(merge) ? NamedTuple() : merge, tolerance, startTime, stopTime, interval, desiredResultTimeUnit, interp_points,
                  adaptive, log, logStates, logEvents, logParameters, logEvaluatedParameters,
                  requiredFinalStates, extra_kwargs)
        return success ? obj : nothing
    end

    SimulationOptions{FloatType,TimeType}(; kwargs...) where {FloatType,TimeType} = 
        SimulationOptions{FloatType,TimeType}(NamedTuple(); kwargs...)
end



"""
    simulationModel = SimulationModel{FloatType,TimeType}(
            modelModule, modelName, getDerivatives!, equationInfo, x_startValues,
            parameters, variableNames;
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
mutable struct SimulationModel{FloatType,TimeType}
    modelModule::Module
    modelName::String
    options::SimulationOptions
    getDerivatives!::Function
    equationInfo::ModiaBase.EquationInfo
    linearEquations::Vector{ModiaBase.LinearEquations{FloatType}}
    eventHandler::EventHandler{FloatType,TimeType}      
    variables::OrderedDict{String,Int}      # Dictionary of variables and their result indices (negated alias has negativ index)
    zeroVariables::OrderedSet{String}       # Set of variables that are identically to zero
    vSolvedWithInitValuesAndUnit::OrderedDict{String,Any}   # Dictionary of (names, init values with units) for all explicitly solved variables with init-values defined
    #p::AbstractVector                       # Parameter and init/start values 
                                            # (parameter values are used in getDerivatives!,
                                            #  init/start values are extracted and stored in x_start) 
    parameters::NamedTuple
    evaluatedParameters::NamedTuple   
    previous::AbstractVector                # previous[i] is the value of previous(...., i)
    nextPrevious::AbstractVector            # nextPrevious[i] is the current value of the variable identified by previous(...., i)
    previous_names::Vector{String}          # previous_names[i] is the name of previous-variable i
    previous_dict::OrderedDict{String,Int}  # previous_dict[name] is the index of previous-variable name

    pre::AbstractVector
    nextPre::AbstractVector
    pre_names::Vector{String}
    pre_dict::OrderedDict{String,Int}

    hold::AbstractVector
    nextHold::AbstractVector
    hold_names::Vector{String}
    hold_dict::OrderedDict{String,Int}
    
    separateObjects::OrderedDict{Int,Any}   # Dictionary of separate objects   
    isInitial::Bool    
    storeResult::Bool
    time::TimeType
    nGetDerivatives::Int                    # Number of getDerivatives! calls
    x_start::Vector{FloatType}              # States x before first event iteration (before initialization)
    x_init::Vector{FloatType}               # States x after initialization (and before integrator is started)
    der_x::Vector{FloatType}                # Derivatives of states x or x_init 
    result::Vector{Tuple}                   # Simulation result
    algorithmType::Union{DataType,Missing}  # Type of integration algorithm (used in default-heading of plot)
    solution::Union{Any,Missing}            # Return value of DifferentialEquations.solve
    save_x_in_solution::Bool                # = true, if the states are stored in solution
     

    function SimulationModel{FloatType,TimeType}(modelModule, modelName, getDerivatives!, equationInfo, x_startValues,
                                        previousVars, preVars, holdVars,
                                        parameterDefinition, variableNames;
                                        nz::Int = 0,
                                        nAfter::Int = 0,
                                        vSolvedWithInitValuesAndUnit::AbstractDict = OrderedDict{String,Any}(),
                                        vEliminated::Vector{Int} = Int[],
                                        vProperty::Vector{Int}   = Int[],
                                        var_name::Function       = v -> nothing) where {FloatType,TimeType}       
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
            
            # Build dictionary of x_names and set start indices of x-vector
            ModiaBase.updateEquationInfo!(equationInfo)

            # Build previous-arrays
            previous       = Vector{Any}(missing, length(previousVars))
            previous_names = string.(previousVars)
            previous_dict  = OrderedDict{String,Int}(zip(previous_names, 1:length(previousVars)))
            
            # Build pre-arrays
            pre       = Vector{Any}(missing, length(preVars))
            pre_names = string.(preVars)
            pre_dict  = OrderedDict{String,Int}(zip(pre_names, 1:length(preVars)))
          
            # Build hold-arrays
            hold       = Vector{Any}(missing, length(holdVars))
            hold_names = string.(holdVars)
            hold_dict  = OrderedDict{String,Int}(zip(hold_names, 1:length(holdVars)))
            
        # Construct parameter values that are copied into the code
        #parameterValues = [eval(p) for p in values(parameters)]
        #@show typeof(parameterValues)
        #@show parameterValues
        parameters = deepcopy(parameterDefinition[:_p])
        
        # Determine x_start and previous values
        nx = equationInfo.nx
        x_start = Vector{FloatType}(undef,nx)
        evaluatedParameters = propagateEvaluateAndInstantiate!(modelModule, parameters, equationInfo, x_start, previous_dict, previous, pre_dict, pre, hold_dict, hold) 
        if isnothing(evaluatedParameters)
            return nothing
        end
        nextPrevious = deepcopy(previous)
        nextPre      = deepcopy(pre)
        nextHold     = deepcopy(hold)
        
        # Construct data structure for linear equations
        linearEquations = ModiaBase.LinearEquations{FloatType}[]
        for leq in equationInfo.linearEquations
            push!(linearEquations, ModiaBase.LinearEquations{FloatType}(leq...))
        end
        
        # Construct dictionary for separate objects
        separateObjects = OrderedDict{Int,Any}()
                
        # Initialize execution flags
        eventHandler = EventHandler{FloatType,TimeType}(nz=nz, nAfter=nAfter)
        eventHandler.initial = true
        isInitial   = true
        storeResult = false
        nGetDerivatives = 0
        
        new(modelModule, modelName, SimulationOptions{FloatType,TimeType}(), getDerivatives!, equationInfo, linearEquations, 
            eventHandler, variables, zeroVariables,
            vSolvedWithInitValuesAndUnit2, parameters, evaluatedParameters, #parameterValues,
            previous, nextPrevious, previous_names, previous_dict,
            pre, nextPre, pre_names, pre_dict,   
            hold, nextHold, hold_names, hold_dict,              
            separateObjects, isInitial, storeResult, convert(TimeType, 0), nGetDerivatives, 
            x_start, zeros(FloatType,nx), zeros(FloatType,nx), Tuple[], missing, missing, false)
    end
    
    
    function SimulationModel{FloatType,TimeType}(m::SimulationModel) where {FloatType,TimeType}       
        # Construct data structure for linear equations
        linearEquations = ModiaBase.LinearEquations{FloatType}[]
        for leq in m.equationInfo.linearEquations
            push!(linearEquations, ModiaBase.LinearEquations{FloatType}(leq...))
        end
        
        # Construct dictionary for separate objects
        separateObjects = OrderedDict{Int,Any}()
                
        # Initialize execution flags
        eventHandler = EventHandler{FloatType,TimeType}()
        eventHandler.initial = true
        isInitial   = true
        storeResult = false
        nGetDerivatives = 0
        nx = m.equationInfo.nx   
        
         new(m.modelModule, m.modelName, m.options, m.getDerivatives!, m.equationInfo, linearEquations, 
            eventHandler, m.variables, m.zeroVariables,
            m.vSolvedWithInitValuesAndUnit, deepcopy(m.parameters), deepcopy(m.evaluatedParameters), 
            deepcopy(m.previous), deepcopy(m.nextPrevious), m.previous_names, m.previous_dict,
            deepcopy(m.pre), deepcopy(m.nextPre), m.pre_names, m.pre_dict,
            deepcopy(m.hold), deepcopy(m.nextHold), m.hold_names, m.hold_dict,            
            separateObjects, isInitial, storeResult, convert(TimeType, 0), nGetDerivatives, 
            convert(Vector{FloatType}, m.x_start), zeros(FloatType,nx), zeros(FloatType,nx), 
            Tuple[], missing, missing, false)       
    end    
end

# Default constructors
SimulationModel(args...; kwargs...) = SimulationModel{Float64,Float64}(args...; kwargs...)
  
SimulationModel{Measurements.Measurement{T}}(args...; kwargs...) where {T} = SimulationModel{Measurements.Measurement{T},T}(args...; kwargs...)
SimulationModel{MonteCarloMeasurements.Particles{T,N}}(args...; kwargs...) where {T,N} = SimulationModel{MonteCarloMeasurements.Particles{T,N},T}(args...; kwargs...)
SimulationModel{MonteCarloMeasurements.StaticParticles{T,N}}(args...; kwargs...) where {T,N} = SimulationModel{MonteCarloMeasurements.StaticParticles{T,N},T}(args...; kwargs...)
SimulationModel{FloatType}(args...; kwargs...) where {FloatType} = SimulationModel{FloatType,FloatType}(args...; kwargs...)

positive(m::SimulationModel, args...; kwargs...) = TinyModia.positive!(m.eventHandler, args...; kwargs...)
negative(m::SimulationModel, args...; kwargs...) = TinyModia.negative!(m.eventHandler, args...; kwargs...)
change(  m::SimulationModel, args...; kwargs...) = TinyModia.change!(  m.eventHandler, args...; kwargs...)
edge(    m::SimulationModel, args...; kwargs...) = TinyModia.edge!(    m.eventHandler, args...; kwargs...)
after(   m::SimulationModel, args...; kwargs...) = TinyModia.after!(   m.eventHandler, args...; kwargs...)
pre(     m::SimulationModel, i)                  = m.pre[i]


"""
    v_zero = reinit(instantiatedModel, _x, j, v_new, _leqMode; v_threshold=0.01)
    
Re-initialize state j with v_new, that is set x[i] = v_new, with i = instantiatedModel.equationInfo.x_info[j].startIndex.
If v_new <= v_threshold, set x[i] to the floating point number that is closest to zero, so that x[i] > 0 and return
v_zero = true. Otherwise return v_zero = false.

An error is triggered if the function is not called during an event phase
or if leqMode >= 0 (which means that reinit is called inside the for-loop
to solve a linear equation system - this is not yet supported).

# Implementation notes

At the beginning of getDerivatives!(..), assignments of _x to appropriate variables are performed.
Therefore, setting _x[i] afterwards via reinit, has no immediate effect on this model evaluation.
With reinit, instantiatedModel.eventHandler.newEventIteration = true is set, to force a new
event iteration. At the next event iteration, the new value v_new is copied from _x and therefore
has then an effect.
"""
function reinit(m::SimulationModel, x, j, v_new, leqMode; v_threshold=0.01)
    if leqMode >= 0
        error("reinit(..) of model ", m.modelName, " is called when solving a linear equation system (this is not supported)")
    elseif !isEvent(m)
        error("reinit(..) of model ", m.modelName, " is called outside of an event phase (this is not allowed)")
    end
   
    eh = m.eventHandler   
    eh.restart = max(eh.restart, Restart)
    eh.newEventIteration = true
    
    i = m.equationInfo.x_info[j].startIndex
    if v_new <= v_threshold
        x[i] = nextfloat(convert(typeof(v_new), 0))
        if eh.logEvents
            println("        State ", m.equationInfo.x_info[j].x_name, " reinitialized to ", x[i], " (reinit returns true)")
        end  
        return true
    else
        x[i] = v_new
        if eh.logEvents
            println("        State ", m.equationInfo.x_info[j].x_name, " reinitialized to ", v_new, " (reinit returns false)")
        end        
        return false
    end
end




"""
    floatType = getFloatType(simulationModel::SimulationModel)

Return the floating point type with which `simulationModel` is parameterized
(for example returns: `Float64, Float32, DoubleFloat, Measurements.Measurement{Float64}`).
"""
getFloatType(m::SimulationModel{FloatType,TimeType}) where {FloatType,TimeType} = FloatType


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
        value = get_value(m.evaluatedParameters, name)
        if ismissing(value)
            @info "get_lastValue: $name is not known and is ignored."
            return nothing;
        end
    end

    return unit ? value : stripUnit(value)
end


"""
    get_extraSimulateKeywordArgumentsDict(instantiateModel)
    
Return the dictionary, in which extra keyword arguments of the last 
`simulate!(instantiatedModel,...)` call are stored.
"""
get_extraSimulateKeywordArgumentsDict(m::SimulationModel) = m.options.extra_kwargs



"""
    terminateEventIteration!(instantiatedModel)
    
Return true, if event iteration shall be terminated and simulation started.

Return false, if event iteration shall be continued. Reasons for continuing are:

1. pre != nextPre

2. positive(..), negative(..), change(..), edge(..) trigger a new iteration.

3. triggerEventAfterInitial = true

Note

- isInitial:
  When simulation shall be started, isInitial = true. When neither (1) nor (2) hold, isInitial = false
  and isInitial remains false for the rest of the simulation.

- firstInitialOfAllSegments: 
  When simulation shall be started, firstInitialOfAllSegments = true. After the first iteration,
  this flag is set to false. This flag can be used to initialize devices (e.g. renderer), that should be
  initialized once and not several times.
  
- firstEventIteration:
  This flag is true during the first iteration at an event.
  
- firstEventIterationDirectlyAfterInitial:
  This flag is true during the first iteration directly at the time event after initialization
  (at simulation startTime).
"""
function terminateEventIteration!(m::SimulationModel{FloatType,TimeType})::Bool where {FloatType,TimeType}
    h = m.eventHandler
    h.firstInitialOfAllSegments = false          
    h.firstEventIteration       = false    
    h.firstEventIterationDirectlyAfterInitial = false

    # pre-iteration
    pre_iteration = false
    for i = 1:length(m.pre)
        if m.pre[i] != m.nextPre[i]
            if m.options.logEvents
                println("        ", m.pre_names[i], " changed from ", m.pre[i], " to ", m.nextPre[i])
            end
            m.pre[i] = m.nextPre[i]
            pre_iteration = true
        end 
    end 
    
    if pre_iteration || h.newEventIteration 
        h.newEventIteration = false
        return false   # continue event iteration
    end
            
    if h.restart == Terminate || h.restart == FullRestart
        h.initial = false
        return true
        
    elseif h.triggerEventDirectlyAfterInitial
        h.triggerEventDirectlyAfterInitial = false
        if h.initial
            h.firstEventIterationDirectlyAfterInitial = true
        end
        h.initial = false
        return false # continue event iteration
    end
    
    h.initial = false
    return true
end


function eventIteration!(m::SimulationModel{FloatType,TimeType}, x::Vector{FloatType}, t_event::TimeType)::Nothing where {FloatType,TimeType}
    eh = m.eventHandler

    # Initialize event iteration
    initEventIteration!(eh, t_event)

    # Perform event iteration
	iter_max = 20
    iter     = 0
	success  = false
    eh.event = true
    while !success && iter <= iter_max
        iter += 1
        Base.invokelatest(m.getDerivatives!, m.der_x, x, m, t_event)      
        success = terminateEventIteration!(m)
        
        if eh.firstEventIterationDirectlyAfterInitial
            iter = 0  # reset iteration counter, since new event
            if m.options.logEvents
                println("\n      Time event at time = ", eh.time, " s")
            end
            eh.nTimeEvents += 1 
            eh.nextEventTime = m.options.stopTime
        end
    end
    eh.event = false        

    if !success
        error("Maximum number of event iterations (= $iter_max) reached")
    end
        
    return nothing
end

"""
    get_xNames(instantiatedModel)
    
Return the names of the elements of the x-vector in a Vector{String}.
"""
get_xNames(m::SimulationModel) = ModiaBase.get_xNames(m.equationInfo)



"""
    isInitial(instantiatedModel)

Return true, if **initialization phase** of simulation.
"""
isInitial(m::SimulationModel) = m.eventHandler.initial
initial(  m::SimulationModel) = m.eventHandler.initial


"""
    isTerminal(instantiatedModel)

Return true, if **terminal phase** of simulation.
"""
isTerminal(m::SimulationModel) = m.eventHandler.terminal
terminal(  m::SimulationModel) = m.eventHandler.terminal


"""
    isEvent(instantiatedModel)

Return true, if **event phase** of simulation (including initialization).
"""
isEvent(m::SimulationModel) = m.eventHandler.event


"""
    isFirstEventIteration(instantiatedModel)

Return true, if **event phase** of simulation (including initialization) and
during the first iteration of the event iteration. 
"""
isFirstEventIteration(m::SimulationModel) = m.eventHandler.firstEventIteration


"""
    isFirstEventIterationDirectlyAfterInitial(instantiatedModel)

Return true, if first iteration directly after initialization where initial=true
(so at the startTime of the simulation). 
"""
isFirstEventIterationDirectlyAfterInitial(m::SimulationModel) = m.eventHandler.firstEventIterationDirectlyAfterInitial


"""
    isAfterSimulationStart(instantiatedModel)

Return true, if **after start of simulation** (returns false during initialization).
"""
isAfterSimulationStart(m::SimulationModel) = m.eventHandler.afterSimulationStart


"""
    isZeroCrossing(instantiatedModel)

Return true, if **event indicators (zero crossings) shall be computed**.
"""
isZeroCrossing(m::SimulationModel) = m.eventHandler.crossing


"""
    storeResults(instantiatedModel)

Return true, if **results shall be stored**.
"""
storeResults(m::SimulationModel) = m.storeResult


isFirstInitialOfAllSegments(m::SimulationModel) = m.eventHandler.firstInitialOfAllSegments
isTerminalOfAllSegments(m::SimulationModel)     = m.eventHandler.isTerminalOfAllSegments


"""
    setNextEvent!(instantiatedModel, nextEventTime)

At an event instant, set the next time event to `nextEventTime`.
"""
setNextEvent!(m::SimulationModel{FloatType,TimeType}, nextEventTime) where {FloatType,TimeType} = 
        setNextEvent!(m.eventHandler, convert(TimeType,nextEventTime)) 


"""
    tCurrent = getTime(instantiatedModel)

Return current simulation time.
"""
getTime(m::SimulationModel) = m.time 


"""
    zStartIndex = addZeroCrossings(instantiatedModel, nz)

Add nz new zero crossing functions and return the start index with respect to 
instantiatedModel.eventHandler.z.
"""
function addZeroCrossings(m::SimulationModel, nz::Int)::Int
    eh = m.eventHandler
    zStartIndex = eh.nz + 1
    eh.nz += nz
    resize!(eh.z, eh.nz)
    resize!(eh.zPositive, eh.nz)
    return zStartIndex
end


get_xe(x, xe_info) = xe_info.length == 1 ? x[xe_info.startIndex] : x[xe_info.startIndex:xe_info.startIndex + xe_info.length-1]

#function set_xe!(x, xe_info, value)::Nothing
#    if xe_info.length == 1 
#        x[xe_info.startIndex] = value
#    else
#        x[xe_info.startIndex:xe_info.startIndex + xe_info.length-1] = value
#    end
#    return nothing
#end



"""
    success = init!(simulationModel)


Initialize `simulationModel::SimulationModel` at `startTime`. In particular:

- Empty result data structure.

- Merge parameter and init/start values into simulationModel.

- Construct x_start.

- Call simulationModel.getDerivatives! once with isInitial(simulationModel) = true to
  compute and store all variables in the result data structure at `startTime`
  and initialize simulationModel.linearEquations.

- Check whether explicitly solved variables that have init-values defined,
  have the required value after initialization (-> otherwise error).
  
If initialization is successful return true, otherwise false.
"""
function init!(m::SimulationModel)::Bool
    empty!(m.result)
    eh = m.eventHandler
    reinitEventHandler(eh, m.options.stopTime, m.options.logEvents)
    
    # Initialize auxiliary arrays for event iteration
    m.x_init .= 0
    m.der_x  .= 0
    
	# Apply updates from merge Map and propagate/instantiate/evaluate the resulting evaluatedParameters 
    if !isnothing(merge)
        m.parameters = recursiveMerge(m.parameters, m.options.merge)
        m.evaluatedParameters = propagateEvaluateAndInstantiate!(m.modelModule, m.parameters, m.equationInfo, m.x_start, m.previous_dict, m.previous, m.pre_dict, m.pre, m.hold_dict, m.hold)
        if isnothing(m.evaluatedParameters)
            return false
        end
    end       
        
    # Re-initialize dictionary of separate objects
    empty!(m.separateObjects)
    
    # Log parameters
    if m.options.logParameters
        parameters = m.parameters
        @showModel parameters
    end
    if m.options.logEvaluatedParameters
        evaluatedParameters = m.evaluatedParameters
        @showModel evaluatedParameters
    end
	
    if m.options.logStates
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

    # Initialize model, linearEquations and compute and store all variables at the initial time
    if m.options.log
        println("      Initialization at time = ", m.options.startTime, " s")  
    end  
 
    # Perform initial event iteration 
    m.nGetDerivatives = 0
    m.isInitial   = true
    eh.initial    = true
#    m.storeResult = true
#    m.getDerivatives!(m.der_x, m.x_start, m, startTime)
#    Base.invokelatest(m.getDerivatives!, m.der_x, m.x_start, m, startTime)
    for i in eachindex(m.x_init)
        m.x_init[i] = deepcopy(m.x_start[i])
    end
    eventIteration!(m, m.x_init, m.options.startTime)
    eh.initial    = false
    m.isInitial   = false
    m.storeResult = false   
    eh.afterSimulationStart = true
    
    # Check vSolvedWithInitValuesAndUnit
    if length(m.vSolvedWithInitValuesAndUnit) > 0
        # Store variables after initialization
        m.storeResult = true
        m.getDerivatives!(m.der_x, m.x_start, m, m.options.startTime) 
        m.storeResult = false 
        
        names = String[]
        valuesBeforeInit = Any[]
        valuesAfterInit  = Any[]
        tolerance = m.options.tolerance

        for (name, valueBefore) in m.vSolvedWithInitValuesAndUnit
            valueAfterInit  = get_lastValue(m, name, unit=false)
            valueBeforeInit = stripUnit.(valueBefore)
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
            println("\n")
            return false
        end
        
        # Remove result 
        empty!(m.result)
    end
    return true
end


"""
    terminate!(m::SimulationModel, x, time)

Terminate model.
"""
function terminate!(m::SimulationModel, x, t)::Nothing
    eh = m.eventHandler
    eh.terminal = true
    m.storeResult = true
    Base.invokelatest(m.getDerivatives!, m.der_x, x, m, t)
    m.storeResult = false
    eh.terminal = false
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
#    println("... Store result at time = $t: s = ", x[1], ", v = ", x[2])    
#    m.getDerivatives!(m.der_x, x, m, t)
    Base.invokelatest(m.getDerivatives!, m.der_x, x, m, t)

    m.storeResult = false
    return nothing
end


"""
    affect_outputs!(integrator)

DifferentialEquations PresetTimeCallback function for `SimulationModel`
that is used to store results at communication points.
"""
function affect_outputs!(integrator)::Nothing
    m = integrator.p
    m.storeResult = true
#    m.getDerivatives!(m.der_x, x, m, t)
    Base.invokelatest(m.getDerivatives!, m.der_x, integrator.u, m, integrator.t)
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
    affectEvent!(integrator)
    
Called when a time or state event is triggered   
"""
function affectEvent!(integrator)::Nothing
    m  = integrator.p
    eh = m.eventHandler
    time = integrator.t
    #println("... begin affect: time = ", time, ", nextEventTime = ", eh.nextEventTime)
                        
    # Compute and store outputs before processing the event
    outputs!(integrator.u, time, integrator)

    # Event iteration 
    eventIteration!(m, integrator.u, time)

    if eh.restart == Restart || eh.restart == FullRestart
        eh.nRestartEvents += 1
    end
    if eh.logEvents
        println("        restart = ", eh.restart)
    end
                
    # Adapt step size    
    if eh.restart != NoRestart && supertype(typeof(integrator.alg)) == DifferentialEquations.OrdinaryDiffEq.OrdinaryDiffEqAdaptiveAlgorithm
        DifferentialEquations.auto_dt_reset!(integrator)
        DifferentialEquations.set_proposed_dt!(integrator, integrator.dt)
    end
    
    # Compute outputs and store them after the event occurred
    outputs!(integrator.u, time, integrator)
    
    # Set next time event
    if abs(eh.nextEventTime - m.options.stopTime) < 1e-10
        eh.nextEventTime = m.options.stopTime
    end
    #println("... end affect: time ", time,", nextEventTime = ", eh.nextEventTime)
    if eh.nextEventTime <= m.options.stopTime
        DifferentialEquations.add_tstop!(integrator, eh.nextEventTime) 
    end
    return nothing
end



"""
    stateEventCondition!(z, x, t, integrator)
    
Called by integrator to compute zero crossings
"""
function stateEventCondition!(z, x, t, integrator)::Nothing
    m = integrator.p
    eh = m.eventHandler
    eh.nZeroCrossings += 1
    eh.crossing = true
    Base.invokelatest(m.getDerivatives!, m.der_x, x, m, t)
    eh.crossing = false
    for i = 1:eh.nz
        z[i] = eh.z[i]
    end 
    #println("... time = ", t, ", z = ", z)
    return nothing
end


"""
    affectStateEvent!(integrator, event_index)
    
Called by integrator when a state event is triggered   
"""
function affectStateEvent!(integrator, event_index)::Nothing
    m  = integrator.p
    eh = m.eventHandler
    time = integrator.t
    
    # Event iteration   
    if eh.logEvents
        println("\n      State event (zero-crossing) at time = ", time, " s")
    end
    affectEvent!(integrator)
    eh.nStateEvents += 1      
    return nothing
end


"""
    timeEventCondition!(u, t, integrator)
    
Called by integrator to check if a time event occurred
"""
function timeEventCondition!(u, t, integrator)::Bool
    m  = integrator.p
    eh = m.eventHandler
    if t >= eh.nextEventTime
        eh.nextEventTime = m.options.stopTime
        return true
    end
    return false   
end


"""
    affectTimeEvent!(integrator)
    
Called by integrator when a time event is triggered   
"""
function affectTimeEvent!(integrator)::Nothing
    m  = integrator.p
    eh = m.eventHandler
    time = integrator.t
    
    # Event iteration   
    if eh.logEvents
        println("\n      Time event at time = ", time, " s")
    end
    affectEvent!(integrator)
    eh.nTimeEvents += 1     
    
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

- pre:Vector{Symbol}: pre-variable names

- `hasUnits::Bool`: = true, if variables have units. Note, the units of the state vector are defined in equationinfo.
"""
function generate_getDerivatives!(AST::Vector{Expr}, equationInfo::ModiaBase.EquationInfo,
                                  parameters, variables, previousVars, preVars, holdVars, functionName::Symbol;
                                  pre::Vector{Symbol} = Symbol[], hasUnits=false)

    # Generate code to copy x to struct and struct to der_x
    x_info     = equationInfo.x_info
    code_x     = Expr[]
    code_der_x = Expr[]
    #code_p     = Expr[]

    if length(x_info) == 1 && x_info[1].x_name == "_dummy_x" && x_info[1].der_x_name == "der(_dummy_x)"
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
                push!(code_der_x, :( _der_x[$indexRange] = TinyModia.stripUnit( $der_x_name )) )
            else
                push!(code_der_x, :( _der_x[$indexRange] = $der_x_name ))
            end
        end
    end
    #for (i,pi) in enumerate(parameters)
    #    push!(code_p, :( $pi = _m.p[$i] ) )
    #end

    timeName = variables[1]
    if hasUnits
        code_time = :( $timeName = _time*upreferred(u"s") )
    else
        code_time = :( $timeName = _time )
    end

    # Code for previous
    code_previous = []
    if length(previousVars) > 0  
        code_previous2 = Expr[]
        for (i, value) in enumerate(previousVars)
            previousName = previousVars[i]
            push!(code_previous2, :( _m.nextPrevious[$i] = $previousName ))        
        end
        code_previous3 = quote
             if TinyModia.isFirstEventIteration(_m) && !TinyModia.isInitial(_m)
                 $(code_previous2...)
             end
        end
        push!(code_previous, code_previous3)
    end

    # Code for pre
    code_pre = Expr[]    
    for (i, value) in enumerate(preVars)
        preName = preVars[i]
        push!(code_pre, :( _m.nextPre[$i] = $preName ))        
    end
    
    # Generate code of the function
    code = quote
                function $functionName(_der_x, _x, _m, _time)::Nothing
                    _m.time = TinyModia.getValue(_time)
                    _m.nGetDerivatives += 1
                    instantiatedModel = _m
                    _p = _m.evaluatedParameters
                    _leq_mode = nothing
                    $code_time
                    $(code_x...)                  
                    $(AST...)
                    $(code_der_x...)
                    $(code_previous...)
                    $(code_pre...)

                    if _m.storeResult
                        TinyModia.addToResult!(_m, $(variables...))
                    end

                    return nothing
                end
            end
    return code
end

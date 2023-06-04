# License for this file: MIT (expat)
# Copyright 2020-2021, DLR Institute of System Dynamics and Control

using  OrderedCollections: OrderedDict, OrderedSet
using  DataFrames

#=
fieldnames(typeof(integrator)) = (:sol, :u, :du, :k, :t, :dt, :f, :p, :uprev, :uprev2, :duprev, :tprev, :alg,
                                  :dtcache, :dtchangeable, :dtpropose, :tdir, :eigen_est, :EEst, :qold, :q11,
                                  :erracc, :dtacc, :success_iter, :iter, :saveiter, :saveiter_dense, :cache,
                                  :callback_cache, :kshortsize, :force_stepfail, :last_stepfail, :just_hit_tstop,
                                  :do_error_check, :event_last_time, :vector_event_last_time, :last_event_error,
                                  :accept_step, :isout, :reeval_fsal, :u_modified, :reinitialize, :isdae, :opts,
                                  :destats, :initializealg, :fsalfirst, :fsallast)
=#

"""
    baseType(T)

Return the base type of a type T, according to the following definition.

```julia
baseType(::Type{T})                                           where {T}     = T
baseType(::Type{Unitful.Quantity{T,D,U}})                     where {T,D,U} = T
baseType(::Type{Measurements.Measurement{T}})                 where {T}     = T
baseType(::Type{MonteCarloMeasurements.Particles{T,N}})       where {T,N}   = T
baseType(::Type{MonteCarloMeasurements.StaticParticles{T,N}}) where {T,N}   = T
```

# Examples
```
baseType(Float32)                # Float32
baseType(Measurement{Float64})   # Float64
```
"""
baseType(::Type{T})                                           where {T}     = T
baseType(::Type{Unitful.Quantity{T,D,U}})                     where {T,D,U} = T
baseType(::Type{Measurements.Measurement{T}})                 where {T}     = T
baseType(::Type{MonteCarloMeasurements.Particles{T,N}})       where {T,N}   = T
baseType(::Type{MonteCarloMeasurements.StaticParticles{T,N}}) where {T,N}   = T

isQuantity(              ::Type{T}) where {T} = T <: Unitful.Quantity         || T <: MonteCarloMeasurements.AbstractParticles && baseType(T) <: Unitful.Quantity
isMeasurements(          ::Type{T}) where {T} = T <: Measurements.Measurement || T <: Unitful.Quantity && baseType(T) <: Measurements.Measurement
isMonteCarloMeasurements(::Type{T}) where {T} = T <: MonteCarloMeasurements.AbstractParticles

Base.floatmax(::Type{Unitful.Quantity{T,D,U}})                     where {T,D,U} = Base.floatmax(T)
Base.floatmax(::Type{Measurements.Measurement{T}})                 where {T}     = Base.floatmax(T)
Base.floatmax(::Type{MonteCarloMeasurements.Particles{T,N}})       where {T,N}   = Base.floatmax(T)
Base.floatmax(::Type{MonteCarloMeasurements.StaticParticles{T,N}}) where {T,N}   = Base.floatmax(T)


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
                        "This is currently not support in Modia.", bold=true, color=:red)
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
        [:merge, :tolerance, :startTime, :stopTime, :interval, :interp_points, :dtmax, :adaptive, :nlinearMinForDAE,
         :log, :logStates, :logEvents, :logProgress, :logTiming, :logParameters, :logEvaluatedParameters,
         :requiredFinalStates, :requiredFinalStates_rtol, :requiredFinalStates_atol, :useRecursiveFactorizationUptoSize])
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


mutable struct SimulationOptions{FloatType,TimeType}
    merge::OrderedDict{Symbol,Any}
    tolerance::Float64
    startTimeFirstSegment::TimeType   # u"s"; startTime of first segment
    startTime::TimeType               # u"s"; startTime of actual segment
    stopTime::TimeType                # u"s"
    interval::TimeType                # u"s"
    desiredResultTimeUnit
    interp_points::Int
    dtmax::Float64
    adaptive::Bool
    nlinearMinForDAE::Int
    log::Bool
    logStates::Bool
    logEvents::Bool
    logProgress::Bool
    logTiming::Bool
    logParameters::Bool
    logEvaluatedParameters::Bool
    requiredFinalStates::Union{Missing, Vector{FloatType}}
    requiredFinalStates_rtol::Float64
    requiredFinalStates_atol::Float64
    useRecursiveFactorizationUptoSize::Int
    extra_kwargs::OrderedDict{Symbol,Any}

    function SimulationOptions{FloatType,TimeType}(merge, errorMessagePrefix=""; _dummy=false, kwargs...) where {FloatType,TimeType}
        success   = true
        adaptive  = get(kwargs, :adaptive, true)
        tolerance = _dummy ? Inf : get(kwargs, :tolerance, 1e-6)
        if tolerance <= 0.0
            printstyled(errorMessagePrefix, "tolerance (= $(tolerance)) must be > 0\n\n", bold=true, color=:red)
            success = false
        elseif tolerance < 100*eps(FloatType) && adaptive
            newTolerance = max(tolerance, 100*eps(FloatType))
            printstyled("Warning from Modia.simulate!(..):\n"*
                        "tolerance (= $(tolerance)) is too small for FloatType = $FloatType (eps(FloatType) = $(eps(FloatType))).\n" *
                        "tolerance changed to $newTolerance.\n\n", bold=true, color=:red)
            tolerance = newTolerance
        end
        startTimeFirstSegment = convertTimeVariable(TimeType, get(kwargs, :startTime, 0.0) )
        startTime             = startTimeFirstSegment
        rawStopTime = get(kwargs, :stopTime, startTime)
        stopTime    = convertTimeVariable(TimeType, rawStopTime)
        interval    = convertTimeVariable(TimeType, get(kwargs, :interval , (stopTime - startTime)/500.0) )
        dtmax       = get(kwargs, :dtmax, 100*getValueOnly(interval))
        if ismissing(dtmax) || isnothing(dtmax)
            dtmax = 100*getValueOnly(interval)
        end
        dtmax = convert(Float64, dtmax)
        desiredResultTimeUnit = unit(rawStopTime)
        interp_points = get(kwargs, :interp_points, 0)
        if interp_points < 0
            printstyled(errorMessagePrefix, "interp_points (= $(interp_points)) must be > 0\n\n", bold=true, color=:red)
            success = false
        elseif interp_points == 1
            # DifferentialEquations.jl crashes
            interp_points = 2
        end
        nlinearMinForDAE = max(1, get(kwargs, :nlinearMinForDAE, 10))   # >= 10
        adaptive         = get(kwargs, :adaptive     , true)
        log              = get(kwargs, :log          , false)
        logStates        = get(kwargs, :logStates    , false)
        logEvents        = get(kwargs, :logEvents    , false)
        logProgress      = get(kwargs, :logProgress  , false)
        logTiming        = get(kwargs, :logTiming    , false)
        logParameters    = get(kwargs, :logParameters, false)
        logEvaluatedParameters   = get(kwargs, :logEvaluatedParameters  , false)
        requiredFinalStates      = get(kwargs, :requiredFinalStates     , missing)
        if isnothing(requiredFinalStates)
            requiredFinalStates = missing
        end
        requiredFinalStates_rtol = get(kwargs, :requiredFinalStates_rtol, 1e-3)
        requiredFinalStates_atol = get(kwargs, :requiredFinalStates_atol, 0.0)
        useRecursiveFactorizationUptoSize = get(kwargs, :useRecursiveFactorizationUptoSize, 0)
        extra_kwargs = OrderedDict{Symbol,Any}()
        for option in kwargs
            key = option.first
            if key in BasicSimulationKeywordArguments
                nothing;
            elseif key in RegisteredExtraSimulateKeywordArguments
                extra_kwargs[key] = option.second
            else
                printstyled(errorMessagePrefix, "simulate!(..., $key=$(option.second)): $key is unknown.\n\n", bold=true, color=:red)
                success = false
            end
        end

#        obj = new(isnothing(merge) ? NamedTuple() : merge, tolerance, startTime, stopTime, interval, desiredResultTimeUnit, interp_points,
        obj = new(ismissing(merge) || isnothing(merge) ? OrderedDict{Symbol,Any}() : merge, tolerance, startTimeFirstSegment, startTime, stopTime, interval, desiredResultTimeUnit, interp_points,
                  dtmax, adaptive, nlinearMinForDAE, log, logStates, logEvents, logProgress, logTiming, logParameters, logEvaluatedParameters,
                  requiredFinalStates, requiredFinalStates_rtol, requiredFinalStates_atol, useRecursiveFactorizationUptoSize, extra_kwargs)
        return success ? obj : nothing
    end

    SimulationOptions{FloatType,TimeType}() where {FloatType,TimeType} =
        SimulationOptions{FloatType,TimeType}(OrderedDict{Symbol,Any}(); _dummy=true)
end


struct LinearEquationsCopyInfoForDAEMode
    ileq::Int           # Index of LinearEquations()
    index::Vector{Int}  # Copy: leq.x[i]                = dae_der_x[index[i]]
                        #       dae_residuals[index[i]] = leq.residuals[i]

    LinearEquationsCopyInfoForDAEMode(ileq) = new(ileq, Int[])
end


"""
    simulationModel = InstantiatedModel{FloatType,TimeType}(
            modelModule, modelName, getDerivatives!, equationInfo, x_startValues,
            parameters, timeName, w_invariant_names;
            vSolvedWithInitValuesAndUnit::OrderedDict{String,Any}(),
            vEliminated::Vector{Int}=Int[],
            vProperty::Vector{Int}=Int[],
            var_name::Function = v->nothing)


# Arguments

- `modelModule`: Module in which `@instantiateModel` is invoked (it is used for `Core.eval(modelModule, ...)`),
  that is evaluation of expressions in the environment of the user.
- `modelName::String`: Name of the model
- `getDerivatives::Function`: Function that is used to evaluate the model equations,
  typically generated with [`Modia.generate_getDerivatives!`].
- `equationInfo::Modia.EquationInfo`: Information about the states and the equations.
- `x_startValues`:: Deprecated (is no longer used).
- `parameters`: A hierarchical NamedTuple of (key, value) pairs defining the parameter and init/start values.
- `timeName`: Name of time (as Symbol)
- `w_invariant_names`: A vector of variable names (as vector of symbols or Expr)
"""
mutable struct InstantiatedModel{FloatType,TimeType}
    # Available before propagateEvaluateAndInstantiate!(..) called (= partiallyInstantedModel)
    modelModule::Module
    modelName::String
    buildDict::OrderedDict{String,Any}
    timer::TimerOutputs.TimerOutput
    cpuFirst::UInt64                  # cpu time of start of simulation
    cpuLast::UInt64                   # Last time from time_ns()
    options::SimulationOptions{FloatType,TimeType}
    getDerivatives!::Function
    initialEquationInfo::Modia.EquationInfo
    linearEquations::Vector{Modia.LinearEquations{FloatType}}
    eventHandler::EventHandler{FloatType,TimeType}
    vSolvedWithInitValuesAndUnit::OrderedDict{String,Any}   # Dictionary of (names, init values with units) for all explicitly solved variables with init-values defined

    previous::AbstractVector                # previous[i] is the value of previous(...., i)
    previous_names::Vector{String}          # previous_names[i] is the name of previous-variable i
    previous_dict::OrderedDict{String,Int}  # previous_dict[name] is the index of previous-variable name

    pre::AbstractVector
    pre_names::Vector{String}
    pre_dict::OrderedDict{String,Int}

    hold::AbstractVector
    hold_names::Vector{String}
    hold_dict::OrderedDict{String,Int}

    isInitial::Bool
    solve_leq::Bool                         # = true, if linear equations 0 = A*x-b shall be solved
                                            # = false, if leq.x is provided by DAE solver and leq.residuals is used by the DAE solver.
    odeMode::Bool                           # = false: copy der(x) into linear equation systems that have leq.odeMode=false and do not solve these equation systems
    storeResult::Bool
    time::TimeType
    nf_total::Int                               # Number of getDerivatives! calls
    nf_integrator::Int                          # Number of getDerivatives! calls from integrator (without zero-crossing calls)
    odeIntegrator::Bool                         # = true , if ODE integrator used
                                                # = false, if DAE integrator used

    daeCopyInfo::Vector{LinearEquationsCopyInfoForDAEMode}  # Info to efficiently copy between DAE and linear equation systems
    algorithmName::Union{String,Missing}        # Name of integration algorithm as string (used in default-heading of plot)
    sundials::Bool                              # = true, if algorithm is a Sundials integrator
    addEventPointsDueToDEBug::Bool              # = true, if event points are explicitly stored for Sundials integrators, due to bug in DifferentialEquations
                                                #         (https://github.com/SciML/Sundials.jl/issues/309)
    success::Bool                               # = true, if after first outputs!(..) call and no error was triggered
                                                # = false, either before first outputs!(..) call or at first outputs!(..) after init!(..) and
                                                #          an error was triggered and simulate!(..) should be returned with nothing.
    unitless::Bool                              # = true, if simulation is performed without units.

    timeName::String
    w_invariant_names::Vector{String}
    hideResult_names::Vector{String}            # Names of hidden variables
    vEliminated::Vector{Int}
    vProperty::Vector{Int}
    var_name::Function
    result::Union{Result,Missing}               # Result data structure upto current time instant

    parameters::OrderedDict{Symbol,Any}         # Parameters as provided to InstantiatedModelconstructor
    equationInfo::Modia.EquationInfo            # Invariant part of equations are available
    x_terminate::Vector{FloatType}              # States x used at the last terminate!(..) call or [], if terminate!(..) not yet called.

    # Available before init!(..) called
    statistics::OrderedDict{Symbol,Any}         # Statistics of simulation run

    # Available after propagateEvaluateAndInstantiate!(..) called
    instantiateFunctions::Vector{Tuple{Union{Expr,Symbol},OrderedDict{Symbol,Any},String}}
                                                # All definitions `_initSegmentFunction = Par(functionName = XXX)` in the model to call
                                                # `XXX(instantiatedModel, submodel, submodelPath)` in the order occurring  during evaluation
                                                # of the parameters where, instantiatedFunctions[i] = (XXX, submodel, submodelPath)
    nsegments::Int                               # Current simulation segment
    evaluatedParameters::OrderedDict{Symbol,Any}     # Evaluated parameters
    nextPrevious::AbstractVector                     # nextPrevious[i] is the current value of the variable identified by previous(...., i)
    nextPre::AbstractVector
    nextHold::AbstractVector

    x_vec::Vector{Vector{FloatType}}            # x_vec[i] holds the actual values of (visible) state vector element
                                                # equationInfo.x_info[equationInfo.nx_info_fixedLength+i:equationInfo.nx_info_invariant]
    x_start::Vector{FloatType}                  # States x before first event iteration (before initialization)
    x_init::Vector{FloatType}                   # States x after initialization (and before integrator is started)
    x_segmented::Vector{FloatType}              # A copy of the current segmented states

    der_x_invariant::Vector{FloatType}          # Derivatives of states x or x_init that correspond to invariant states
                                                # This vector is filled with derivatives of invariants states with appendVariable!(m.der_x_invariant, ...) calls,
                                                # including derivatives of x_vec[i]
    der_x_segmented::Vector{FloatType}          # Derivatives of states x or x_init that correspond to segmented states (defined in functions and not visible in getDerivatives!(..))
    der_x::Vector{FloatType}                    # Derivatives of states x


    function InstantiatedModel{FloatType,TimeType}(modelModule, modelName, buildDict, getDerivatives!, equationInfo,
                                        previousVars, preVars, holdVars,
                                        parameterDefinition, timeName, w_invariant_names, hideResult_names;
                                        unitless::Bool=true,
                                        nz::Int = 0,
                                        nAfter::Int = 0,
                                        vSolvedWithInitValuesAndUnit::AbstractDict = OrderedDict{String,Any}(),
                                        vEliminated::Vector{Int} = Int[],
                                        vProperty::Vector{Int}   = Int[],
                                        var_name::Function       = v -> nothing) where {FloatType,TimeType}

        # Construct data structure for linear equations
        linearEquations = Modia.LinearEquations{FloatType}[]
        for leq in equationInfo.linearEquations
            push!(linearEquations, Modia.LinearEquations{FloatType}(leq...))
        end

        vSolvedWithInitValuesAndUnit2 = OrderedDict{String,Any}( [(string(key),vSolvedWithInitValuesAndUnit[key]) for key in keys(vSolvedWithInitValuesAndUnit)] )

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

        # Initialize execution flags
        eventHandler  = EventHandler{FloatType,TimeType}(nz=nz, nAfter=nAfter)
        isInitial     = true
        storeResult   = false
        solve_leq     = true
        nf_total      = 0
        nf_integrator = 0
        odeIntegrator = true
        daeCopyInfo = LinearEquationsCopyInfoForDAEMode[]
        algorithmName = missing
        sundials = false
        addEventPointsDueToDEBug = false
        success  = false

        w_invariant_names = String[string(name) for name in w_invariant_names]

        # Initialize other data
        result = missing
        instantiateResult = true
        newResultSegment = false
        parameters = deepcopy(parameterDefinition)
        x_terminate = FloatType[]

       new(modelModule, modelName, buildDict, TimerOutputs.TimerOutput(), UInt64(0), UInt64(0), SimulationOptions{FloatType,TimeType}(), getDerivatives!,
           equationInfo, linearEquations, eventHandler,
           vSolvedWithInitValuesAndUnit2,
           previous, previous_names, previous_dict,
           pre, pre_names, pre_dict,
           hold, hold_names, hold_dict,
           isInitial, solve_leq, true, storeResult, convert(TimeType, 0), nf_total, nf_integrator,
           odeIntegrator, daeCopyInfo, algorithmName, sundials, addEventPointsDueToDEBug, success, unitless,
           string(timeName), w_invariant_names, hideResult_names, vEliminated, vProperty, var_name, result,
           parameters, equationInfo, x_terminate)
    end

#=
    function InstantiatedModel{FloatType,TimeType}(m::InstantiatedModel) where {FloatType,TimeType}
        # Construct data structure for linear equations
        linearEquations = Modia.LinearEquations{FloatType}[]
        for leq in m.equationInfo.linearEquations
            push!(linearEquations, Modia.LinearEquations{FloatType}(leq...))
        end

        # Initialize execution flags
        eventHandler = EventHandler{FloatType,TimeType}()
        eventHandler.initial = true
        isInitial     = true
        storeResult   = false
        solve_leq     = true
        nf_total      = 0
        nf_integrator = 0
        odeIntegrator = true
        daeCopyInfo = LinearEquationsCopyInfoForDAEMode[]
        algorithmName = missing
        success = false
        nx          = m.equationInfo.nx
        nxInvariant = m.equationInfo.nxInvariant
        nxSegmented   = nx-nxInvariant

        emptyResult!(m.result)
        result   = deepcopy(m.result)
        nsegments = 1

        x_vec           = [zeros(FloatType, equationInfo.x_info[i].length) for i in equationInfo.nx_info_fixedLength+1:equationInfo.nx_info_invariant]
        x_start         = convert(Vector{FloatType}, m.x_start)
        x_init          = zeros(FloatType, nx)
        x_segmented       = zeros(FloatType, nxSegmented)
        der_x_invariant = zeros(FloatType, nxInvariant)
        der_x_segmented   = zeros(FloatType, nxSegmented)
        der_x           = zeros(FloatType, nx)

        new(m.modelModule, m.modelName, deepcopy(m.buildDict), TimerOutputs.TimerOutput(), UInt64(0), UInt64(0), deepcopy(m.options), m.getDerivatives!,
            deepcopy(m.equationInfo), deepcopy(linearEquations),
            deepcopy(eventHandler),
            m.vSolvedWithInitValuesAndUnit,
            deepcopy(m.previous), m.previous_names, m.previous_dict,
            deepcopy(m.pre), m.pre_names, m.pre_dict,
            deepcopy(m.hold), m.hold_names, m.hold_dict,
            isInitial, solve_leq, true, storeResult, convert(TimeType, 0), nf_total, nf_integrator,
            true, LinearEquationsCopyInfoForDAEMode[],
            odeIntegrator, daeCopyInfo, m.sundials, m.addEventPointsDueToDEBug, success, m.unitless,
            result, nsegments,
            deepcopy(m.parameters), deepcopy(m.instantiateFunctions), deepcopy(m.evaluatedParameters),
            deepcopy(m.nextPrevious), deepcopy(m.nextPre), deepcopy(m.nextHold),
            x_vec, x_start, x_init, x_segmented, der_x_invariant, der_x_segmented, der_x)
    end
=#
end


# Default constructors
InstantiatedModel{FloatType}(args...; kwargs...) where {FloatType} = InstantiatedModel{FloatType,FloatType}(args...; kwargs...)

InstantiatedModel{Measurements.Measurement{T},}(args...; kwargs...) where {T} = InstantiatedModel{Measurements.Measurement{T},T}(args...; kwargs...)
InstantiatedModel{MonteCarloMeasurements.Particles{T,N}}(args...; kwargs...) where {T,N,} = InstantiatedModel{MonteCarloMeasurements.Particles{T,N},T}(args...; kwargs...)
InstantiatedModel{MonteCarloMeasurements.StaticParticles{T,N}}(args...; kwargs...) where {T,N} = InstantiatedModel{MonteCarloMeasurements.StaticParticles{T,N},T}(args...; kwargs...)

timeType(m::InstantiatedModel{FloatType,TimeType}) where {FloatType,TimeType} = TimeType

# The following rule is important for DiffEqBase version 6.91.6 and later
# (https://github.com/SciML/DiffEqBase.jl/issues/791)
if Base.isdefined(DiffEqBase, :anyeltypedual)
    DiffEqBase.anyeltypedual(::InstantiatedModel) = Any
end

positive(m::InstantiatedModel, args...; kwargs...) = Modia.positive!(m.eventHandler, args...; kwargs...)
negative(m::InstantiatedModel, args...; kwargs...) = Modia.negative!(m.eventHandler, args...; kwargs...)
change(  m::InstantiatedModel, args...; kwargs...) = Modia.change!(  m.eventHandler, args...; kwargs...)
edge(    m::InstantiatedModel, args...; kwargs...) = Modia.edge!(    m.eventHandler, args...; kwargs...)
after(   m::InstantiatedModel, args...; kwargs...) = Modia.after!(   m.eventHandler, args...; kwargs...)
pre(     m::InstantiatedModel, i)                  = m.pre[i]


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
function reinit(m::InstantiatedModel, x, j, v_new, leqMode; v_threshold=0.01)
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
    floatType = getFloatType(simulationModel::InstantiatedModel)

Return the floating point type with which `simulationModel` is parameterized
(for example returns: `Float64, Float32, DoubleFloat, Measurements.Measurement{Float64}`).
"""
getFloatType(m::InstantiatedModel{FloatType,TimeType}) where {FloatType,TimeType} = FloatType


"""
    hasParticles(value)

Return true, if `value` is of type `MonteCarloMeasurements.StaticParticles` or
`MonteCarloMeasurements.Particles`.
"""
hasParticles(value) = typeof(value) <: MonteCarloMeasurements.StaticParticles ||
                      typeof(value) <: MonteCarloMeasurements.Particles


"""
    get_value(obj, name::String)

Return the value identified by `name` from the potentially hierarchically
dictionary `obj`. If `name` is not in `obj`, the function returns `missing`.

# Examples
```julia
s1 = Map(a = 1, b = 2, c = 3)
s2 = Map(v1 = s1, v2 = (d = 4, e = 5))
s3 = Map(v3 = s2, v4 = s1)

@show get_value(s3, "v3.v1.b")   # returns 2
@show get_value(s3, "v3.v2.e")   # returns 5
@show get_value(s3, "v3.v1.e")   # returns missing
```
"""
function get_value(obj::OrderedDict, name::String)
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
        if haskey(obj,key) && typeof(obj[key]) <: OrderedDict && length(name) > j
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
    get_names(obj)

Return `Vector{String}` containing all the names present in `obj`

# Examples
```julia
s1 = Map(a = 1, b = 2, c = 3)
s2 = Map(v1 = s1, v2 = (d = 4, e = 5))
s3 = Map(v3 = s2, v4 = s1)

@show get_names(s3)
```
"""
function get_names(obj) # ::NamedTuple)
    names = String[]
    get_names!(obj, names, "")
    return names
end
function get_names!(obj #= ::NamedTuple =# , names::Vector{String}, path::String)::Nothing
    for (key,value) in obj # zip(keys(obj), obj)
        if key != :_class
            name = appendName(path, key)
            if typeof(value) <: OrderedDict
                get_names!(value, names, name)
            else
                push!(names, name)
            end
        end
    end
    return nothing
end


"""
    getLastValue(model::InstantiatedModel, name::String; unit=true)

Return the last stored value of variable `name` from `model`.
If `unit=true` return the value with its unit, otherwise with stripped unit.

If `name` is not known or no result values yet available, an info message is printed
and the function returns `nothing`.

`name` can be a time-varying variable or a parameter.
"""
function getLastValue(m::InstantiatedModel{FloatType,TimeType}, name::String; unit::Bool=true) where {FloatType,TimeType}
    if isnothing(m) || ismissing(m.result)
        @info "getLastValue(model,\"$name\"): No results yet available."
        return nothing
    end

    result = m.result
    if haskey(result.info, name)
        # Time varying variable stored in m.result
        resInfo = result.info[name]
        if length(result.t) == 0
            @info "getLastValue(model,\"$name\"): No results yet available."
            return nothing
        end

        if resInfo.kind == RESULT_ELIMINATED
            value = getLastValue(m, resInfo.aliasName, unit=unit)
            if resInfo.aliasNegate
                value *= -1
            end

        elseif resInfo.kind == RESULT_CONSTANT
            value = resInfo.value
            if !unit
                value = stripUnit(value)
            end

        elseif resInfo.kind == RESULT_T
            value = result.t[end][end]
            if unit && resInfo.unit != ""
                value *= uparse(resInfo.unit)
            end

        elseif resInfo.kind == RESULT_X
            id = resInfo.id[end]
            segment = id.segment
            if segment < length(result.t)
                return missing
            end
            ibeg = id.index
            iend = ibeg + prod(id.size) - 1
            value = ibeg == iend ? result.x[segment][end][ibeg] : result.x[segment][end][ibeg:iend]
            if unit && resInfo.unit != ""
                value *= uparse(resInfo.unit)
            end

        elseif resInfo.kind == RESULT_DER_X
            id = resInfo.id[end]
            segment = id.segment
            if segment < length(result.t)
                return missing
            end
            ibeg = id.index
            iend = ibeg + prod(id.size) - 1
            value = ibeg == iend ? result.der_x[segment][end][ibeg] : result.der_x[segment][end][ibeg:iend]
            if unit && resInfo.unit != ""
                value *= uparse(resInfo.unit)
            end

        elseif resInfo.kind == RESULT_W_INVARIANT
            id = resInfo.id[end]
            value = result.w_invariant[end][end][id.index]
            if !unit
                value = stripUnit(value)
            end

        elseif resInfo.kind == RESULT_W_SEGMENTED
            id = resInfo.id[end]
            segment = id.segment
            if segment < length(result.t)
                return missing
            end
            value = result.w_segmented[segment][end][id.index]
            if unit && resInfo.unit != "" && eltype(value) <: AbstractFloat
                value *= uparse(resInfo.unit)
            end

        else
            error("Bug in getLastValue(...), name = $name, resInfo.kind = $(resInfo.kind).")
        end

    else
        # Parameter stored in m.evaluatedParameters
        value = get_value(m.evaluatedParameters, name)
        if ismissing(value)
            @info "getLastValue: $name is not known and is ignored."
            return nothing;
        end
    end

    return value
end
get_lastValue(m, name; unit=true) = getLastValue(m, name, unit=unit)


"""
    get_extraSimulateKeywordArgumentsDict(instantiateModel)

Return the dictionary, in which extra keyword arguments of the last
`simulate!(instantiatedModel,...)` call are stored.
"""
get_extraSimulateKeywordArgumentsDict(m::InstantiatedModel) = m.options.extra_kwargs



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
function terminateEventIteration!(m::InstantiatedModel)::Bool
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


function eventIteration!(m::InstantiatedModel, x, t_event)::Nothing
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
        #Base.invokelatest(m.getDerivatives!, m.der_x, x, m, t_event)
        invokelatest_getDerivatives_without_der_x!(x, m, t_event)
        eh.firstInitialOfAllSegments = false
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
get_xNames(m::InstantiatedModel) = Modia.get_xNames(m.equationInfo)


"""
    isInitial(instantiatedModel)

Return true, if **initialization phase** of simulation
(of the current segment of a segmented simulation).
"""
isInitial(m::InstantiatedModel) = m.eventHandler.initial
initial(  m::InstantiatedModel) = m.eventHandler.initial


"""
    isFirstInitialOfAllSegments(instantiatedModel)

Return true, if **initialization phase** of simulation of the **first segment**
of a segmented simulation.
"""
isFirstInitialOfAllSegments(m::InstantiatedModel) = m.eventHandler.firstInitialOfAllSegments


"""
    isTerminal(instantiatedModel)

Return true, if **terminal phase** of simulation
(of the current segment of a segmented simulation).
"""
isTerminal(m::InstantiatedModel) = m.eventHandler.terminal
terminal(  m::InstantiatedModel) = m.eventHandler.terminal


"""
    isTerminalOfAllSegments(instantiatedModel)

Return true, if **terminal phase** of simulation of the **last segment**
of a segmented simulation.
"""
isTerminalOfAllSegments(m::InstantiatedModel) = m.eventHandler.terminalOfAllSegments


"""
    isEvent(instantiatedModel)

Return true, if **event phase** of simulation (including initialization).
"""
isEvent(m::InstantiatedModel) = m.eventHandler.event


"""
    isFirstEventIteration(instantiatedModel)

Return true, if **event phase** of simulation (including initialization) and
during the first iteration of the event iteration.
"""
isFirstEventIteration(m::InstantiatedModel) = m.eventHandler.firstEventIteration


"""
    isFirstEventIterationDirectlyAfterInitial(instantiatedModel)

Return true, if first iteration directly after initialization where initial=true
(so at the startTime of the simulation).
"""
isFirstEventIterationDirectlyAfterInitial(m::InstantiatedModel) = m.eventHandler.firstEventIterationDirectlyAfterInitial


"""
    isFullRestart(instantiatedModel)

Return true, if **FullRestart event** of a segmented simulation.
"""
isFullRestart(m::InstantiatedModel) = m.eventHandler.fullRestart


"""
    isAfterSimulationStart(instantiatedModel)

Return true, if **after start of simulation** (returns false during initialization).
"""
isAfterSimulationStart(m::InstantiatedModel) = m.eventHandler.afterSimulationStart


"""
    isZeroCrossing(instantiatedModel)

Return true, if **event indicators (zero crossings) shall be computed**.
"""
isZeroCrossing(m::InstantiatedModel) = m.eventHandler.crossing


"""
    storeResults(instantiatedModel)

Return true, if **results shall be stored**.
"""
storeResults(m::InstantiatedModel) = m.storeResult


"""
    setNextEvent!(instantiatedModel, nextEventTime)

At an event instant, set the next time event to `nextEventTime`.
"""
setNextEvent!(m::InstantiatedModel{FloatType,TimeType}, nextEventTime) where {FloatType,TimeType} =
        setNextEvent!(m.eventHandler, convert(TimeType,nextEventTime))


function setFullRestartEvent!(m::InstantiatedModel)
    m.eventHandler.restart = FullRestart
    return nothing
end


"""
    tCurrent = getTime(instantiatedModel)

Return current simulation time.
"""
getTime(m::InstantiatedModel) = m.time


get_xe(x, xe_info) = xe_info.length == 1 ? x[xe_info.startIndex] : x[xe_info.startIndex:xe_info.startIndex + xe_info.length-1]

#function set_xe!(x, xe_info, value)::Nothing
#    if xe_info.length == 1
#        x[xe_info.startIndex] = value
#    else
#        x[xe_info.startIndex:xe_info.startIndex + xe_info.length-1] = value
#    end
#    return nothing
#end
import Printf

function invokelatest_getDerivatives_without_der_x!(x, m, t)::Nothing
    TimerOutputs.@timeit m.timer "Modia getDerivatives!" begin
        if m.options.logProgress && m.cpuLast != UInt64(0)
            cpuNew = time_ns()
            if (cpuNew - m.cpuLast) * 1e-9 > 5.0
                m.cpuLast = cpuNew
                Printf.@printf("      progress: integrated up to time = %.3g s (in cpu-time = %.3g s)\n", t, (cpuNew-m.cpuFirst)*1e-9)
            end
        end
        if length(m.x_vec) > 0
            # copy vector-valued x-elements from x to m.x_vec
            eqInfo = m.equationInfo
            x_vec  = m.x_vec
            j = 1
            for i in eqInfo.nx_info_fixedLength+1:eqInfo.nx_info_invariant
                xe = eqInfo.x_info[i]
                copyto!(x_vec[j], 1, x, xe.startIndex, xe.length)
                j += 1
            end
        end
        copyto!(m.x_segmented, 1, x, m.equationInfo.nxInvariant+1, m.equationInfo.nxSegmented)
        empty!(m.der_x_invariant)
        m.der_x_segmented .= 0   # Handles also the case for a dummy differential equation (if model has no states): der(_dummy_x) = 0.0
        Base.invokelatest(m.getDerivatives!, x, m, t)

        @assert(length(m.der_x_invariant) + length(m.der_x_segmented) == length(x))
    end
    return nothing
end


"""
    success = init!(simulationModel)

Initialize `simulationModel::InstantiatedModel` at `startTime`. In particular:

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
function init!(m::InstantiatedModel{FloatType,TimeType})::Bool where {FloatType,TimeType}
    # Initialize model, linearEquations and compute and store all variables at the initial time
    if m.options.log
        println("      Initialization at time = ", m.options.startTime, " s")
    end

    m.equationInfo         = deepcopy(m.initialEquationInfo)
    equationInfo           = m.equationInfo
    eh                     = m.eventHandler
    m.instantiateFunctions = Tuple{Union{Expr,Symbol}, OrderedDict{Symbol,Any}, String}[]
    m.nsegments            = 1
    m.statistics           = newStatistics(m.options.startTime, get_algorithmName_for_heading(m), m.odeIntegrator, FloatType)
    m.result = Result{FloatType,TimeType}(m.equationInfo, m.timeName, m.w_invariant_names, m.vEliminated, m.vProperty, m.var_name)
    result   = m.result

    if length(m.options.merge) > 0
        m.parameters = mergeModels(m.parameters, m.options.merge)
    end
    if m.options.logParameters
        parameters = m.parameters
        @showModel parameters
    end
    evaluatedParameters = propagateEvaluateAndInstantiate!(m, log=false)
    if isnothing(evaluatedParameters)
        return false
    end
    if m.options.logEvaluatedParameters
        @showModel evaluatedParameters
    end
    m.evaluatedParameters = evaluatedParameters
    m.nextPrevious = deepcopy(m.previous)
    m.nextPre      = deepcopy(m.pre)
    m.nextHold     = deepcopy(m.hold)
    m.x_start      = initialStateVector!(m)
    useRecursiveFactorizationUptoSize = m.options.useRecursiveFactorizationUptoSize
    for leq in m.linearEquations
        leq.useRecursiveFactorizationUptoSize = useRecursiveFactorizationUptoSize
        leq.useRecursiveFactorization         = length(leq.x) <= useRecursiveFactorizationUptoSize && length(leq.x) > 0
    end

    # update equationInfo
    for xi_info in equationInfo.x_info
        resInfo = result.info[xi_info.x_name]
        resInfo.signal[:start] = xi_info.startOrInit
        id = xi_info.x_segmented_startIndex == -1 ? ValuesID(         -1, xi_info.startIndex, size(xi_info.startOrInit)) :
                                                    ValuesID(m.nsegments, xi_info.startIndex, size(xi_info.startOrInit))
        push!(resInfo.id, id)
        resInfo = result.info[xi_info.der_x_name]
        push!(resInfo.id, id)
    end

    # Provide storage for x and der_x utility vectors
    nx                = length(m.x_start)
    nxSegmented       = nx-equationInfo.nxInvariant
    m.x_vec           = [zeros(FloatType, equationInfo.x_info[i].length) for i in equationInfo.nx_info_fixedLength+1:equationInfo.nx_info_invariant]
    m.x_init          = zeros(FloatType,nx)
    m.x_segmented     = zeros(FloatType, nxSegmented)
    m.der_x_invariant = zeros(FloatType,equationInfo.nxInvariant)
    m.der_x_segmented = zeros(FloatType, nxSegmented)
    m.der_x           = zeros(FloatType,nx)
    updateStatistics_nStates!(m.statistics, nx)

    if m.options.logStates
        # List init/start values
        x_table = DataFrames.DataFrame(state=String[], init=Any[], unit=String[])   #, nominal=String[])
        for xe_info in m.equationInfo.x_info
            xe_init = get_xe(m.x_start, xe_info)
            if hasParticles(xe_init)
                xe_init = string(minimum(xe_init)) * " .. " * string(maximum(xe_init))
            end
            # xe_nominal = isnan(xe_info.nominal) ? "" : xe_info.nominal
            push!(x_table, (xe_info.x_name, xe_init, xe_info.unit))   #, xe_nominal))
        end
        show(stdout, x_table; allrows=true, allcols=true, rowlabel = Symbol("#"), summary=false, eltypes=false, truncate=60)
        println("\n")
    end

    # Perform initial event iteration
    m.nf_total      = 0
    m.nf_integrator = 0
    m.isInitial = true
    eh.initial  = true
    for i in eachindex(m.x_init)
        m.x_init[i] = deepcopy(m.x_start[i])
    end
    TimerOutputs.@timeit m.timer "Modia eventIteration!" eventIteration!(m, m.x_init, m.options.startTime)
    m.success     = false   # is set to true at the first outputs! call.
    eh.initial    = false
    m.isInitial   = false
    m.storeResult = false
    eh.afterSimulationStart = true

    return true
end


"""
    initFullRestart!(instantiatedModel)

Re-initialize `instantiatedModel` after a `FullRestart` event.
"""
function initFullRestart!(m::InstantiatedModel{FloatType,TimeType})::Nothing where {FloatType,TimeType}
    if m.options.logEvents
        println("      Reinitialization due to FullRestart at time = ", m.time, " s")
    end
    eh = m.eventHandler
    removeSegmentedStates!(m.equationInfo)
    m.nsegments    += 1
    m.nf_total      = 0
    m.nf_integrator = 0
    m.options.startTime = m.time
    reinitEventHandlerForFullRestart!(eh, m.time, m.options.stopTime, m.options.logEvents)
    newResultSegment!(m.result, m.equationInfo, m.nsegments)

    # Evaluate instantiate functions
    for fc in m.instantiateFunctions
        logInstantiatedFunctionCalls = false
        initSegment = fc[1]
        path        = fc[3]
        ID          = path
        parameters  = fc[2]        
        Core.eval(m.modelModule, :($initSegment($m, $path, $ID, $parameters, log=$logInstantiatedFunctionCalls)))
    end
    resizeLinearEquations!(m, m.evaluatedParameters, m.options.log)

    # Get initial state vector
    m.x_start = initialStateVector!(m)

    # update equationInfo
    x_info = m.equationInfo.x_info
    for i = m.equationInfo.nx_info_invariant+1:length(x_info)
        xi_info = x_info[i]
        resInfo = m.result.info[xi_info.x_name]
        id = ValuesID(m.nsegments, xi_info.startIndex, size(xi_info.startOrInit))
        push!(resInfo.id, id)
        resInfo = m.result.info[xi_info.der_x_name]
        push!(resInfo.id, id)
    end

    # Resize states and results
    nx          = length(m.x_start)
    nxSegmented = nx - m.equationInfo.nxInvariant
    resize!(m.x_init         , nx)
    resize!(m.x_segmented    , nxSegmented)
    resize!(m.der_x_invariant, m.equationInfo.nxInvariant)
    resize!(m.der_x_segmented, nxSegmented)
    resize!(m.der_x          , nx)
    m.x_init          .= FloatType(0)
    m.x_segmented     .= FloatType(0)
    m.der_x_invariant .= FloatType(0)
    m.der_x_segmented .= FloatType(0)
    m.der_x           .= FloatType(0)
    updateStatistics_nStates_segmented!(m.statistics, nx)

    if m.options.logStates
        # List init/start values
        x_table = DataFrames.DataFrame(state=String[], init=Any[], unit=String[])   #, nominal=String[])
        for xe_info in m.equationInfo.x_info
            xe_init = get_xe(m.x_start, xe_info)
            if hasParticles(xe_init)
                xe_init = string(minimum(xe_init)) * " .. " * string(maximum(xe_init))
            end
            # xe_nominal = isnan(xe_info.nominal) ? "" : xe_info.nominal
            push!(x_table, (xe_info.x_name, xe_init, xe_info.unit))   #, xe_nominal))
        end
        show(stdout, x_table; allrows=true, allcols=true, rowlabel = Symbol("#"), summary=false, eltypes=false)
        println("\n")
    end


    # Perform event iteration
    m.isInitial   = true
    eh.initial    = true
    for i in eachindex(m.x_init)
        m.x_init[i] = deepcopy(m.x_start[i])
    end
    eventIteration!(m, m.x_init, m.options.startTime)
    eh.fullRestart = false
    eh.initial     = false
    m.isInitial    = false
    m.storeResult  = false
    eh.afterSimulationStart = true
    return nothing
end


"""
    success = check_vSolvedWithInitValuesAndUnit(m::InstantiatedModel)

Check that m.vSolvedWithInitValuesAndUnit is consistent to calculated values.
If this is not the case, print an error message and return false.
"""
function check_vSolvedWithInitValuesAndUnit(m::InstantiatedModel)::Bool
    # Check vSolvedWithInitValuesAndUnit
    if length(m.vSolvedWithInitValuesAndUnit) > 0
        names = String[]
        valuesBeforeInit = Any[]
        valuesAfterInit  = Any[]
        tolerance = m.options.tolerance

        for (name, valueBefore) in m.vSolvedWithInitValuesAndUnit
            valueAfterInit  = getLastValue(m, name, unit=false)
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
    end
    return true
end



"""
    terminate!(m::InstantiatedModel, x, time)

Terminate model.
"""
function terminate!(m::InstantiatedModel, x, t)::Nothing
    #println("... terminate! called at time = $t")
    eh = m.eventHandler
    eh.terminal = true
    eh.terminalOfAllSegments = m.eventHandler.restart != Modia.FullRestart
    invokelatest_getDerivatives_without_der_x!(x, m, t)
    eh.terminal = false
    eh.terminalOfAllSegments = false
    m.x_terminate = deepcopy(x)
    return nothing
end

approxTime(t1,t2) = eps(typeof(t1))

"""
    outputs!(x, t, integrator)

DifferentialEquations FunctionCallingCallback function for `InstantiatedModel`
that is used to store results at communication points.
"""
function outputs!(x, t, integrator)::Nothing
    m = integrator.p
    m.storeResult = true
    #println("... Store result at time = $t")
    if m.odeMode
        invokelatest_getDerivatives_without_der_x!(x, m, t)
    else
        if t==m.options.startTime
            m.der_x .= integrator.du   # Since IDA gives an error for integrator(t, Val{1]}) at the initial time instant
        else
            integrator(m.der_x, t, Val{1})  # Compute derx by interpolation
        end

        # Copy derx to linearEquations
        for copyInfo in m.daeCopyInfo
            leq = m.linearEquations[ copyInfo.ileq ]
            for i in 1:length(copyInfo.index)
                leq.x[i] = m.der_x[ copyInfo.index[i] ]
            end
        end
        m.solve_leq = false
        invokelatest_getDerivatives_without_der_x!(x, m, t)
        m.solve_leq = true
    end
    m.storeResult = false
    if !m.success
        m.success = check_vSolvedWithInitValuesAndUnit(m)
        if !m.success
            DifferentialEquations.terminate!(integrator)
        end
    elseif m.eventHandler.restart == Terminate
        DifferentialEquations.terminate!(integrator)
    end
    return nothing
end


#=
"""
    affect_outputs!(integrator)

DifferentialEquations PresetTimeCallback function for `InstantiatedModel`
that is used to store results at communication points.
"""
function affect_outputs!(integrator)::Nothing
    m = integrator.p
    m.storeResult = true
    m.solve_leq = false
    getDerivatives!2(m.der_x, integrator.u, m, integrator.t)
    m.solve_leq = true
    m.storeResult = false
    if m.eventHandler.restart == Terminate
        DifferentialEquations.terminate!(integrator)
    end
    return nothing
end
=#


"""
    copyDerivatives!(der_x, der_x_invariant, der_x_segmented)

Copy der_x_invariant and der_x_segmented to der_x (der_x .= [der_x_invariant, der_x_segmented])
"""
@inline function copyDerivatives!(der_x, der_x_invariant, der_x_segmented)::Nothing
    if length(der_x_segmented) == 0
        der_x .= der_x_invariant
    else
        @assert(length(der_x) == length(der_x_invariant) + length(der_x_segmented))
        unsafe_copyto!(der_x, 1, der_x_invariant, 1, length(der_x_invariant))
        unsafe_copyto!(der_x, length(der_x_invariant)+1, der_x_segmented, 1, length(der_x_segmented))
    end
    return nothing
end


"""
    derivatives!(derx, x, m, t)

DifferentialEquations callback function to get the derivatives.
"""
function derivatives!(der_x, x, m, t)
    m.nf_integrator += 1
    invokelatest_getDerivatives_without_der_x!(x, m, t)
    #println("t = $t, m.der_x_segmented = ", m.der_x_segmented)
    copyDerivatives!(der_x, m.der_x_invariant, m.der_x_segmented)
    return nothing
end



"""
    DAEresidualsForODE!(residuals, derx, x, m, t)

DifferentialEquations callback function for DAE integrator for ODE model
"""
function DAEresidualsForODE!(residuals, derx, x, m, t)::Nothing
    m.nf_integrator += 1

    # Copy derx to linearEquations
    for copyInfo in m.daeCopyInfo
        leq = m.linearEquations[ copyInfo.ileq ]
        for i in 1:length(copyInfo.index)
            leq.x[i] = derx[ copyInfo.index[i] ]
        end
    end

    m.solve_leq = false
    invokelatest_getDerivatives_without_der_x!(x, m, t)
    m.solve_leq = true

    copyDerivatives!(m.der_x, m.der_x_invariant, m.der_x_segmented)
    residuals .= m.der_x .- derx

    # Get residuals from linearEquations
    for copyInfo in m.daeCopyInfo
        leq = m.linearEquations[ copyInfo.ileq ]
        for i in 1:length(copyInfo.index)
            residuals[  copyInfo.index[i] ] = leq.residuals[i]
        end
    end
    return nothing
end



"""
    affectEvent!(integrator, stateEvent, eventIndex)

Called when a time event (stateEvent=false) or state event (stateEvent=true) is triggered.
In case of stateEvent, eventIndex is the index of the crossing function that triggered the event.
"""
function affectEvent!(integrator, stateEvent::Bool, eventIndex::Int)::Nothing
    m  = integrator.p
    eh = m.eventHandler
    time = integrator.t
    #println("... begin affect: time = ", time, ", nextEventTime = ", eh.nextEventTime)

    if m.addEventPointsDueToDEBug
        push!(integrator.sol.t, deepcopy(integrator.t))
        push!(integrator.sol.u, deepcopy(integrator.u))
    end

    # Compute and store outputs before processing the event
    sol_t = integrator.sol.t
    sol_x = integrator.sol.u
    ilast = length(sol_t)
    for i = length(m.result.t[end])+1:ilast    # -1
        outputs!(sol_x[i], sol_t[i], integrator)
    end
    # A better solution is needed.
    #if sol_t[ilast] == sol_t[ilast-1]
    #    # Continuous time instant is present twice because "saveat" and "DiscreteCallback" occur at the same time instant
    #    # Do not compute the model again
    #    push!(m.result_code , m.result_code[end])
    #    push!(m.result_der_x, m.result_der_x[end])
    #else
    #    outputs!(sol_x[ilast], sol_t[ilast], integrator)
    #end

    if stateEvent
        # State event
        if eh.logEvents
            println("\n      State event (zero-crossing) at time = ", time, " s (due to z[$eventIndex])")
        end
        eh.nStateEvents += 1
    else
        # Time event
        if eh.logEvents
            println("\n      Time event at time = ", time, " s")
        end
        eh.nTimeEvents += 1
    end

    # Event iteration
    eventIteration!(m, integrator.u, time)

    if eh.restart == Restart || eh.restart == FullRestart
        eh.nRestartEvents += 1
    end
    if eh.logEvents
        if eh.restart == Restart
            println("        restart = ", eh.restart)
        else
            printstyled("        restart = ", eh.restart, "\n", color=:red)
        end
    end

    # Compute outputs and store them after the event occurred
    if m.addEventPointsDueToDEBug
        push!(integrator.sol.t, deepcopy(integrator.t))
        push!(integrator.sol.u, deepcopy(integrator.u))
    end

    outputs!(integrator.u, time, integrator)
    if eh.restart == Terminate || eh.restart == FullRestart
        DifferentialEquations.terminate!(integrator)
        return nothing
    end

    # Adapt step size
    if eh.restart != NoRestart && supertype(typeof(integrator.alg)) == DifferentialEquations.OrdinaryDiffEq.OrdinaryDiffEqAdaptiveAlgorithm
        DifferentialEquations.auto_dt_reset!(integrator)
        DifferentialEquations.set_proposed_dt!(integrator, integrator.dt)
    end

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
    zeroCrossings!(z, x, t, integrator)

Called by integrator to compute zero crossings
"""
function zeroCrossings!(z, x, t, integrator)::Nothing
    m = integrator.p
    eh = m.eventHandler
    eh.nZeroCrossings += 1
    eh.crossing = true
    m.solve_leq = false # has only an effect for DAE integrators
    invokelatest_getDerivatives_without_der_x!(x, m, t)
    m.solve_leq = true
    eh.crossing = false
    copyto!(z, eh.z)
    #println("... t = $t, z = $z")
    #println("... time = ", t, ", z = ", z)
    return nothing
end


"""
    affectStateEvent!(integrator, event_index)

Called by integrator when a state event is triggered
"""
affectStateEvent!(integrator, event_index) = affectEvent!(integrator, true, event_index)


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
affectTimeEvent!(integrator) = affectEvent!(integrator, false, 0)


"""
    leq = initLinearEquationsIteration!(m::InstantiatedModel, leq_index::Int)

Initialize iteration over linear equations system of index `leq_index` and return a reference to it
that can be used in the while-loop to construct and solve the linear equation system.
"""
initLinearEquationsIteration!(m::InstantiatedModel, leq_index::Int) = begin
    leq      = m.linearEquations[leq_index]
    leq.mode = -3
    return leq
end


"""
    resizeLinearEquations!(partiallyInstantiatedModel, evaluatedParameters, log::Bool)

Inspect all linear equations inside the instantiatedModel and resize the
internal storage, of the length of vector-valued elements of the iteration variables
has changed due to changed parameter values.
"""
function resizeLinearEquations!(m::InstantiatedModel{FloatType,TimeType}, evaluatedParameters, log::Bool)::Nothing where {FloatType,TimeType}
    for (i,leq) in enumerate(m.linearEquations)
        nx_vec = length(leq.x_names) - leq.nx_fixedLength
        if nx_vec > 0
            j = 0
            for i = leq.nx_fixedLength+1:length(leq.x_names)
                # Length  of element is determined by start or init value
                xi_name  = leq.x_names[i]
                xi_param = get_value(evaluatedParameters, xi_name)
                if ismissing(xi_param)
                    @error "resizeLinearEquations!(instantiatedModel,$i): $xi_name is not part of the evaluated parameters."
                end
                leq.x_lengths[i] = length(xi_param)
                j += 1
                if length(leq.x_vec[j]) != leq.x_lengths[i]
                    if m.options.logEvents
                        println("      Modia: Resize linear equations vector $xi_name from ", length(leq.x_vec[j]), " to ", leq.x_lengths[i] )
                    end
                    resize!(leq.x_vec[j], leq.x_lengths[i])
                end
            end
            nx = sum(leq.x_lengths)
            if length(leq.x) != nx
                # Resize internal arrays
                resize!(leq.x        , nx)
                resize!(leq.b        , nx)
                resize!(leq.pivots   , nx)
                resize!(leq.residuals, nx)
                leq.A = zeros(FloatType,nx,nx)
                leq.useRecursiveFactorization = nx <= leq.useRecursiveFactorizationUptoSize
            end
        end
    end

    return nothing
end


"""
    addToResult!(instantiatedModel, x, time, w_invariant...)

Add result of current time instant (`time, x, der_x, w_invariant, w_segmented`) to `instantiatedModel`.
"""
function addToResult!(m::InstantiatedModel{FloatType,TimeType}, x, time, w_invariant...)::Nothing where {FloatType,TimeType}
    @assert(length(w_invariant) == m.result.n_w_invariant)
    copyDerivatives!(m.der_x, m.der_x_invariant, m.der_x_segmented)
    result::Result = m.result
    push!(result.t[          end], time)
    push!(result.x[          end], deepcopy(x))
    push!(result.der_x[      end], deepcopy(m.der_x))
    push!(result.w_invariant[end], deepcopy(w_invariant))
    push!(result.w_segmented[end], deepcopy(m.result.w_segmented_temp))
    return nothing
end


"""
    obj = get_instantiatedSubmodel(instantiatedModel, ID)
    
Return reference `obj` to an instantiated submodel struct, given `intantiatedModel` and the `ID` of the submodel.    
"""
get_instantiatedSubmodel(instantiatedModel, ID) = instantiatedModel.buildDict[ID]



"""
    startIndex = new_x_segmented_variable!(
                    partiallyInstantiatedModel::InstantiatedModel,
                    x_name::String, der_x_name::String, startOrInit, x_unit::String="";
                    nominal::Float64 = NaN, unbounded::Bool = false)::Int

Generate new states (`x_segmented` and `der_x_segmented` variables) and return the
`startIndex` of the variables in order that actual values can be inquired or copied from the
state `x` and state derivative `der(x)`vectors via [`get_x_startIndex_from_x_segmented_startIndex`](@ref).
`startOrInit` contain the `start` or `init` values of the newly generated `x_segmented` variable.

Actual values of these new variables are stored in:

- `instantiatedModel.x_segmented[startIndex:startIndex+prod(dims(startOrInit))-1]`
- `instantiatedModel.der_x_segmented[startIndex:startIndex+prod(dims(startOrInit))-1]`

Value `startOrInit` is the start/init value used during re-initialization of the new segment with `initFullRestart!(..)`.
"""
function new_x_segmented_variable!(m::InstantiatedModel{FloatType,TimeType}, x_name::String, der_x_name::String, startOrInit, x_unit::String="";
                                   nominal::Float64 = NaN, unbounded::Bool = false)::Int where {FloatType,TimeType}
    eqInfo = m.equationInfo
    result = m.result
    @assert(eqInfo.status == EquationInfo_Initialized_Before_All_States_Are_Known)
    new_result_info = true
    if haskey(result.info, x_name)
        # State was already defined in one of the previous segments
        new_result_info = false
        x_info = result.info[x_name]
        @assert(x_info.kind == RESULT_X)
        @assert(eltypeOrType(startOrInit) == FloatType)
        @assert(length( x_info.id[end].dims ) == ndims(startOrInit))   # Number of dimensions cannot change
        #if typeof(startOrInit) <: Number
        #    @assert(xi_info.scalar)
        #elseif typeof(startOrInit) <: AbstractVector
        #    @assert(!xi_info.scalar)
        #else
        #    error("new_x_segmented_variable(.. $x_name, $der_x_name, startOrInit,...): typeof(startOrInit) is neither a Number nor an AbstractVector)")
        #end
        @assert(get(x_info.signal,:unit,"") == x_unit)
        @assert(haskey(result.info, der_x_name))
        der_x_info = result.info[der_x_name]
        @assert(der_x_info.kind == RESULT_DER_X)
    end

    x_segmented_startIndex = eqInfo.nxSegmented+1
    if isnothing(startOrInit)
        @info "State $x_name has no start or init value defined. Using start value = 0.0."
        fixed = false
        startOrInit = FloatType(0)
    else
        fixed = true
    end
    xi_info = StateElementInfo(x_name, Symbol(x_name), der_x_name, Symbol(der_x_name),
                               XD, x_unit, startOrInit, fixed, nominal, unbounded,
                               x_segmented_startIndex = x_segmented_startIndex)
    push!(eqInfo.x_info, xi_info)
    x_infoIndex = length(eqInfo.x_info)
    eqInfo.x_dict[x_name]         = x_infoIndex
    eqInfo.der_x_dict[der_x_name] = x_infoIndex
    eqInfo.nxSegmented += length(startOrInit)

    if new_result_info
        # result.info can be only partially instantiated, because x_startIndex is only known
        # after function initialStateVector!(...) was called.
        t_unit = get(result.info[result.timeName].signal, :unit, "")
        der_x_unit = x_unit == "" ? unitAsParseableString(1/uparse(t_unit)) : unitAsParseableString(uparse(x_unit)/uparse(t_unit))
        if x_unit == ""
            x_var = Var(start=xi_info.startOrInit, fixed=xi_info.fixed, state=true, der=xi_info.der_x_name)
        else
            x_var = Var(unit=x_unit, start=xi_info.startOrInit, fixed=xi_info.fixed, state=true, der=xi_info.der_x_name)
        end
        result.info[x_name]     = ResultInfo(RESULT_X, x_var, FloatType)
        if der_x_unit == ""
            result.info[der_x_name] = ResultInfo(RESULT_DER_X, Var(), FloatType)
        else
            result.info[der_x_name] = ResultInfo(RESULT_DER_X, Var(unit=der_x_unit), FloatType)
        end
    end
    return x_segmented_startIndex
end


"""
    x_startIndex = get_x_startIndex_from_x_segmented_startIndex(
                      instantiatedModel::InstantiatedModel, x_segmented_startIndex)

Return the startindex of an `x_segmented` state with respect to the `x`-vector,
given the startIndex with respect to the `x_segmented` vector
(`x_segmented_startIndex` is the return value of `new_x_segmented_variable!(..)`).
"""
get_x_startIndex_from_x_segmented_startIndex(m::InstantiatedModel, x_segmented_startIndex::Int) = m.equationInfo.nxInvariant + x_segmented_startIndex


"""
    index = new_w_segmented_variable!(
               partiallyInstantiatedModel::InstantiatedModel, name::String,
               w_segmented_default, unit::String="")::Int

Generate new local variable (`w_segmented` variable) and return the `index` of the variable 
in order that actual values can be inquired or copied from the result data structure.
New values of `w_segmented` variables need only to be computed at communication points.
Value w_segmented_default is stored as default value and defines type and (fixed) size of the variable
in this simulation segment.
"""
function new_w_segmented_variable!(m::InstantiatedModel, name::String, w_segmented_default, unit::String="")::Int
    result = m.result
    w_size = size(w_segmented_default)
    push!(result.w_segmented_names, name)
    push!(result.w_segmented_temp, deepcopy(w_segmented_default))
    w_index = length(result.w_segmented_temp)

    if haskey(result.info, name)
        # Variable was already defined in one of the previous segments
        v_info = result.info[name]
        @assert(v_info.kind == RESULT_W_SEGMENTED)
        w_unit = get(v_info.signal, :unit, "")
        if w_unit != unit
            error("Variable $name changed unit from \"$w_unit\" to \"$unit\".")
        end
        #@assert(v_info.type == typeof(w_segmented_default))
        push!(v_info.id, ValuesID(m.nsegments, w_index, w_size))
    else
        # Variable is defined the first time in the segmented simulation
        if unit == ""
            signal = Var()
        else
            signal = Var(unit=unit)
        end
        result.info[name] = ResultInfo(RESULT_W_SEGMENTED, signal, ValuesID(m.nsegments, w_index, w_size), eltypeOrType(w_segmented_default))
    end
    #println("new_w_segmented_variable: w_segmented_temp = ", result.w_segmented_temp)
    return w_index
end


"""
    new_alias_segmented_variable!(partiallyInstantiatedModel::InstantiatedModel, 
       name, aliasName, aliasNegate=false)

Define new alias variable.
"""
function new_alias_segmented_variable!(m::InstantiatedModel, name::String, aliasName::String, aliasNegate::Bool=false)::Int
    result = m.result
    w_size = size(w_segmented_default)

    if haskey(result.info, name)
        error("new_alias_segmented_variable!(...): $name is already defined")
    elseif !haskey(result.info, aliasName)
        error("new_alias_segmented_variable!(...): $name should be made an alias to $aliasName, but this name is not defined")
    end
    result.info[name] = ResultInfo(RESULT_ELIMINATED, aliasName, aliasNegate)
    push!(result.alias_segmented_names, name)
    return nothing
end


"""
    startIndex = new_z_segmented_variable!(instantiatedModel, nz)

Generate `nz` new zero crossing variables and return the startIndex to of the variables
in order that actual values can be copied into the vector of zero crossings.
"""
function new_z_segmented_variable!(m::InstantiatedModel{F,TimeType}, nz::Int)::Int where {F,TimeType}
    eh = m.eventHandler
    zStartIndex = eh.nz + 1
    eh.nz += nz
    resize!(eh.z, eh.nz)
    resize!(eh.zPositive, eh.nz)
    for i = zStartIndex:eh.nz
        eh.z[i] = convert(F, 1.0)
        eh.zPositive[i] = false
    end
    return zStartIndex
end


"""
    value = Modia.copy_scalar_x_segmented_value_from_state(instantiatedModel, startIndex)

Return value of scalar `x_segmented` variable from state vector `x` by providing its `startIndex`
(returned from `new_x_segmented_variable!(..)`).
"""
copy_scalar_x_segmented_value_from_state(m::InstantiatedModel, startIndex::Int) = m.x_segmented[startIndex]


"""
    value = Modia.copy_SVector3_x_segmented_value_from_state(instantiatedModel, startIndex)

Return value of `SVector{3,FloatType} x_segmented` variable from state vector `x` by providing its `startIndex`
(returned from `new_x_segmented_variable!(..)`).
"""
@inline copy_SVector3_x_segmented_value_from_state(m::InstantiatedModel{FloatType,TimeType}, startIndex::Int) where {FloatType,TimeType} = begin
    x_segmented = m.x_segmented
    return SVector{3,FloatType}(x_segmented[startIndex], x_segmented[startIndex+1], x_segmented[startIndex+2])
end

"""
    Modia.copy_Vector_x_segmented_value_from_state(
        instantiatedModel::InstantiatedModel, startIndex, xi::Vector{FloatType})::Nothing

Return value of `Vector{FloatType} x_segmented` variable from state vector `x` by providing its `startIndex`
(returned from `new_x_segmented_variable!(..)`) and copying it into the pre-allocated vector `xi`.
"""
@inline function copy_Vector_x_segmented_value_from_state(m::InstantiatedModel{FloatType,TimeType}, startIndex::Int, xi::Vector{FloatType})::Nothing where {FloatType,TimeType}
    copyto!(xi, 1, m.x_segmented, startIndex, length(xi))
    return nothing
end


"""
    Modia.copy_der_x_segmented_value_to_state(
       instantiatedModel, startIndex, 
       der_x_segmented_value::[FloatType|Vector{FloatType}])

Copy `der_x_segmented_value` to state derivative vector `der(x)` by providing its `startIndex`
(returned from `new_x_segmented_variable!(..)`) and copying it into the pre-allocated vector `der_x_segmented_value`.
"""
@inline function copy_der_x_segmented_value_to_state(m::InstantiatedModel{FloatType,TimeType}, startIndex::Int, der_x_segmented_value::FloatType)::Nothing where {FloatType,TimeType}
    m.der_x_segmented[startIndex] = der_x_segmented_value
    return nothing
end
@inline function copy_der_x_segmented_value_to_state(m::InstantiatedModel{FloatType,TimeType}, startIndex::Int, der_x_segmented_value::Vector{FloatType})::Nothing where {FloatType,TimeType}
    copyto!(m.der_x_segmented, startIndex, der_x_segmented_value, 1, length(der_x_segmented_value))
    return nothing
end


"""
    Modia.copy_w_segmented_value_to_result(
        instantiatedModel::InstantiatedModel, index::Int, 
        w_segmented_value)::Nothing

Copy value of local variable (`w-segmented`) to result by providing its `index`
(returned from `new_w_segmented_variable!`),
"""
@inline function copy_w_segmented_value_to_result(m::InstantiatedModel, index::Int, w_segmented_value)::Nothing
    w_segmented_temp = m.result.w_segmented_temp
    @assert(typeof(w_segmented_value) == typeof(w_segmented_temp[index]))
    @assert(size(  w_segmented_value) == size(  w_segmented_temp[index]))
    w_segmented_temp[index] = deepcopy(w_segmented_value)
    return nothing
end


function initialStateVector!(m::InstantiatedModel{FloatType,TimeType})::Vector{FloatType} where {FloatType,TimeType}
    if length(m.equationInfo.x_info) == 0
        # Handle systems with only algebraic variables, by introducing a dummy
        # differential equation der_x[1] = -x[1], with state name _dummy_x
        new_x_segmented_variable!(m, "_dummy_x", "der(_dummy_x)", FloatType(0))
    end
    return initialStateVector!(m.equationInfo, FloatType, !isFullRestart(m), m.x_terminate)
end


"""
    code = generate_getDerivatives!(AST, equationInfo, parameters, timeName, w_invariant_names, functionName;
                                    hasUnits=false)

Return the code of the `getDerivatives!` function as `Expr` using the
Symbol `functionName` as function name. By `eval(code)` or
`fc = @RuntimeGeneratedFunction(code)` the function is compiled and can afterwards be called.

# Arguments

- `AST::Vector{Expr}`: Abstract Syntax Tree of the equations as vector of `Expr`.

- `equationInfo::Modia.EquationInfo`: Data structure returned by `Modia.getSortedAndSolvedAST
            holding information about the states.

- `parameters`: Vector of parameter names (as vector of symbols)

- `timeName`: Name of time (as symbol)

- `w_invariant_names`: Vector of variable names (as vector of symbols or Expr).

- `functionName::Function`: The name of the function that shall be generated.


# Optional Arguments

- pre:Vector{Symbol}: pre-variable names

- `hasUnits::Bool`: = true, if variables have units. Note, the units of the state vector are defined in equationinfo.
"""
function generate_getDerivatives!(FloatType, TimeType, AST::Vector{Expr}, equationInfo::Modia.EquationInfo,
                                  parameters, timeName, w_invariant_names, previousVars, preVars, holdVars, functionName::Symbol;
                                  pre::Vector{Symbol} = Symbol[], hasUnits=false)

    # Generate code to copy x to struct and struct to der_x
    x_info     = equationInfo.x_info
    code_x     = Expr[]
    code_der_x = Expr[]
    #code_p     = Expr[]

    if length(x_info) > 0
        i1 = 1
        for i in 1:equationInfo.nx_info_fixedLength
            xe = x_info[i]
            x_name = xe.x_name_julia
            # x_name     = Meta.parse("m."*xe.x_name)
            # der_x_name = Meta.parse("m."*replace(xe.der_x_name, r"der\(.*\)" => s"var\"\g<0>\""))
            if xe.scalar
                # x-element is a scalar
                if !hasUnits || xe.unit == ""
                    push!(code_x, :( $x_name = _x[$i1] ) )
                else
                    x_unit = xe.unit
                    push!(code_x, :( $x_name = _x[$i1]*@u_str($x_unit)) )
                end
                i1 += 1
            else
                # x-element is a static vector
                i2 = i1 + xe.length - 1
                v_length = xe.length
                x_elements = Expr[]
                for i in i1:i2
                    push!(x_elements, :( _x[$i] ))
                end
                if !hasUnits || xe.unit == ""
                    push!(code_x, :( $x_name = Modia.SVector{$v_length,_FloatType}($(x_elements...)) ))
                else
                    x_unit = xe.unit
                    push!(code_x, :( $x_name = Modia.SVector{$v_length,_FloatType}($(x_elements...))*@u_str($x_unit)) )
                end
                i1 = i2 + 1
            end
        end

        i1 = 0
        for i in equationInfo.nx_info_fixedLength+1:equationInfo.nx_info_invariant
            # x-element is a dynamic vector (length can change before initialization)
            xe     = x_info[i]
            x_name = xe.x_name_julia
            i1    += 1
            if !hasUnits || xe.unit == ""
                push!(code_x, :( $x_name = _m.x_vec[$i1]) )
            else
                x_unit = xe.unit
                push!(code_x, :( $x_name = _m.x_vec[$i1]*@u_str($x_unit)) )
            end
        end

        for i in 1:equationInfo.nx_info_invariant
            xe = equationInfo.x_info[i]
            der_x_name = xe.der_x_name_julia
            if hasUnits
                push!(code_der_x, :( Modia.appendVariable!(_m.der_x_invariant, Modia.stripUnit( $der_x_name )) ))
            else
                push!(code_der_x, :( Modia.appendVariable!(_m.der_x_invariant, $der_x_name) ))
            end
        end
    end

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
             if Modia.isFirstEventIteration(_m) && !Modia.isInitial(_m)
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
    # temporarily removed: _m.time = $TimeType(Modia.getValueOnly(_time))
    code = quote
                function $functionName(_x, _m::Modia.InstantiatedModel{$FloatType,$TimeType}, _time::$TimeType)::Nothing
                    _FloatType = $FloatType
                    _TimeType = $TimeType
                    _m.time = _time
                    _m.nf_total += 1
                    instantiatedModel = _m
                    _p = _m.evaluatedParameters
                    _leq_mode = nothing
                    $code_time
                    $(code_x...)
                    $(AST...)
                    $(code_der_x...)
                    $(code_previous...)
                    $(code_pre...)

                    if Modia.storeResults(_m)
                        Modia.TimerOutputs.@timeit _m.timer "Modia addToResult!" Modia.addToResult!(_m, _x, _time, $(w_invariant_names...))
                    end
                    return nothing
                end
            end
    return code
end


# License for this file: MIT (expat)
# Copyright 2022, DLR Institute of System Dynamics and Control

#--------------------------------------------------------------------------------------------------
#        Provide the overloaded Abstract Signal Tables Interface for the results of InstantiatedModel
#--------------------------------------------------------------------------------------------------

@inline function checkMissingResult(m::InstantiatedModel, name::String)::Bool
    if isnothing(m) || ismissing(m) || ismissing(m.result)
        error("$name: No simulation results available.")
    end
    return true
end

SignalTables.isSignalTable(r::Result) = true
SignalTables.isSignalTable(m::InstantiatedModel) = true


"""
    getIndependentSignalNames(instantiatedModel::Modia.InstantiatedModel|result::Modia.Result)::Vector{String}

Return the name of the independent variable of the result stored in instantiatedModel or in result.
"""
SignalTables.getIndependentSignalNames(result::Result)     = [result.timeName]
SignalTables.getIndependentSignalNames(m::InstantiatedModel) = begin
    if ismissing(m.result)
        error("independentSignalName(..): No simulation results available in instantiated model of $(m.modelName)")
    end
    SignalTables.getIndependentSignalNames(m.result)
end



function getParameterNames!(parameters, names::Vector{String}, path::String)::Nothing
    for (key,value) in parameters
        name = path == "" ? string(key) : path*"."*string(key)
        if typeof(value) <: AbstractDict
            getParameterNames!(value, names, name)
        else
            push!(names, name)
        end
    end
    return nothing
end

function getParameterNames(parameters)::Vector{String}
    names = String[]
    getParameterNames!(parameters, names, "")
    return names
end


function headPartOfName(name::String)::String
    # Determine head part of name (before last ".")
    i = first(something(findlast(".", name), 0:-1))
    return i > 1 && i < length(name) ? name[1:i-1] : name
end

function sortUpperHierarchies(names::Vector{String})
    sortedIndices = sortperm(String[headPartOfName(name) for name in names])
    return names[sortedIndices]
end


"""
    getSignalNames(instantiatedModel::Modia.InstantiatedModel;
                   getVar=true, getPar=true, getMap=true)::Vector{String}

Returns a string vector of the variables of an
[`@instantiateModel`](@ref) that are present in the result or are (evaluated) parameters.
- If getVar=true, Var(..) variables are included.
- If getPar=true, Par(..) variables are included.
- If getMap=true, Map(..) variables are included.
"""
function SignalTables.getSignalNames(m::InstantiatedModel; getVar=true, getPar=true, getMap=true,
                                     var=nothing, par=nothing)::Vector{String}   # backwards compatibility
    # For backwards compatibility
        if !isnothing(var)
            getVar = var
        end
        if !isnothing(par)
            getPar = par
        end

    if ismissing(m.result)
        error("getSignalNames(..): No simulation results available in instantiated model of $(m.modelName)")
    end

    var_names = getVar ? collect(keys(m.result.info)) : String[]
    par_names = getPar ? setdiff(getParameterNames(m.evaluatedParameters), m.hideResult_names) : String[]  # parameters without hideResult_names
    map_names = getMap ? ["_attributes"] : String[]
    if length(var_names) > 1
        sig_names = vcat(map_names, var_names[1], sortUpperHierarchies( union(var_names[2:end], par_names) ) )
    else
        sig_names = vcat(map_names, var_names, sortUpperHierarchies(par_names) )
    end
    return sig_names
end


"""
    getSignalNames(result::Modia.Result;
                   getVar=true, getPar=true, getMap=true)::Vector{String}

Returns a string vector of the variables that are present in the result.
"""
SignalTables.getSignalNames(result::Result; getVar=true, getPar=true, getMap=true) = getVar ? collect(keys(result.info)) : String[]


"""
    getStateNames(instantiatedModel::Modia.InstantiatedModel)::Vector{String}

Returns a string vector of the states that are present in [`@instantiateModel`](@ref)
"""
function getStateNames(m::InstantiatedModel)::Vector{String}
    if ismissing(m.result)
        return collect(keys(m.equationInfo.x_dict))
    else
        stateNames = String[]
        for (key,resultInfo) in m.result.info
            if resultInfo.kind == RESULT_X
                push!(stateNames, key)
            end
        end
        return stateNames
    end
end




"""
    getIndependentSignalsSize(instantiatedModel::Modia.InstantiatedModel|result::Modia.Result)::Dims

Returns the lengths of the independent signal of the result as (len,).
"""
SignalTables.getIndependentSignalsSize(m::InstantiatedModel) = SignalTables.getIndependentSignalsSize(m.result)
SignalTables.getIndependentSignalsSize(result::Result) = (sum(length(sk) for sk in result.t), )

#=
using Pkg

function getPackageInfo(name::String)
    for (key,value) in Pkg.dependencies()
        if value.name == name && value.is_direct_dep
            return (uuid=key, version=value.version, source=value.source)
        end
    end
    return nothing
end
=#


"""
    getSignal(instantiatedModel::Modia.InstantiatedModel|result::Modia.Result, name::String)

Returns signal `name` of the result present in instantiatedModel
(that is a [`SignalTables.Var`](@ref), [`SignalTables.Par`](@ref)) or [`SignalTables.Map`](@ref))
If `name` does not exist, an error is raised.
"""
function SignalTables.getSignal(m::InstantiatedModel, name::String)
    checkMissingResult(m, name)
    result = m.result

    if haskey(result.info, name)
        # name is a result variable (time-varying variable stored in m.result)
        signal = getSignal(result, name)

    elseif name == "_attributes"
        signal = SignalTables.Map(model = SignalTables.Map(name = m.modelName),
                                  experiment = SignalTables.Map(startTime=m.options.startTimeFirstSegment,
                                                                stopTime=m.options.stopTime,
                                                                interval=m.options.interval,
                                                                tolerance=m.options.tolerance,
                                                                algorithm=get_algorithmName_for_heading(m),
                                                                dtmax=m.options.dtmax,
                                                                interp_points=m.options.interp_points,
                                                                adaptive=m.options.adaptive
                                                                ),
                                  statistics = deepcopy(m.statistics),
                                  unitFormat = "Unitful")

    else
        # name might be a parameter
        sigValue = get_value(m.evaluatedParameters, name)
        if ismissing(sigValue)
            error("getSignal(.., $name): name is not known")
        end
        sigUnit = unitAsParseableString(sigValue)
        if sigUnit == ""
            signal = SignalTables.Par(value = sigValue)
        else
            signal = SignalTables.Par(value = ustrip.(sigValue), unit=sigUnit)
        end
    end
    return signal
end

function SignalTables.getSignal(result::Result, name::String)
    resInfo = result.info[name]

    if haskey(resInfo.signal, :values)
        return resInfo.signal
    end

    if resInfo.kind == RESULT_ELIMINATED
        signal      = resInfo.signal
        aliasName   = resInfo.aliasName
        aliasNegate = resInfo.aliasNegate
        aliasSignal = SignalTables.getSignal(result, aliasName)
        signal      = merge(aliasSignal, signal)
        if aliasNegate
            signal[:values] = -deepcopy(aliasSignal[:values])
        else
            signal[:alias]  = aliasName
            signal[:values] = aliasSignal[:values]
        end
        return signal
    end

    if resInfo.kind == RESULT_T
        sigValues = signalResultValues(result.t, result.t, resInfo, result; name=name)
    elseif resInfo.kind == RESULT_X
        sigValues = signalResultValues(result.t, result.x, resInfo, result; name=name)
    elseif resInfo.kind == RESULT_DER_X
        sigValues = signalResultValues(result.t, result.der_x, resInfo, result; name=name)
    elseif resInfo.kind == RESULT_W_INVARIANT
        index   = resInfo.id[1].index
        w_value = result.w_invariant[1][1][index]
        w_unit  = SignalTables.unitAsParseableString(w_value)

        # w_invariant has potentially a unit defined - remove it
            if w_unit != ""
                resInfo.signal[:unit] = w_unit
                try
                    w_value = ustrip.(w_value)
                catch
                    nothing
                end
            end

        # w_invariant is defined in all segments and has no size information defined - add size information
            try
                resInfo.id[1] = ValuesID(index, size(w_value))
            catch
                resInfo.id[1] = ValuesID(index, ())
            end

        resInfo._eltypeOrType = eltypeOrType(w_value)
        sigValues = signalResultValues(result.t, result.w_invariant, resInfo, result; name=name)
    elseif resInfo.kind == RESULT_W_SEGMENTED
        sigValues = signalResultValues(result.t, result.w_segmented, resInfo, result; name=name)
    elseif resInfo.kind == RESULT_CONSTANT
        value = resInfo.value
        t     = getSignal(result, result.timeName)[:values]
        len_t = sum(length(tk) for tk in t)
        if typeof(value) <: AbstractArray
            sigValues = Array{eltype(value), ndims(values)}(undef, (len_t, size(value)...))
            for i = 1:length(sigValues)
                for j = 1:length(value)
                    sigValues[i,j] = value[j]
                end
            end
        else
            sigValues = fill(value, len_t)
        end
    else
        error("Bug in getSignal: name=\"$name\" has ResultInfo=$resInfo, but ResultInfo.kind = $(resInfo.kind) is not known.")
    end

    resInfo.signal[:values] = sigValues
    return resInfo.signal
end


"""
    hasSignal(instantiatedModel::Modia.InstantiatedModel, name::String)

Returns `true` if signal `name` is present in the instantiatedModel result or in the evaluated parameters.
"""
SignalTables.hasSignal(m::InstantiatedModel, name::String) = begin
    if isnothing(m) || ismissing(m) || ismissing(m.result)
        return false
    end
    return haskey(m.result.info, name) || !ismissing(get_value(m.evaluatedParameters, name))
end


"""
    hasSignal(result::Modia.Result, name::String)

Returns `true` if signal `name` is present in result.
"""
SignalTables.hasSignal(result::Result, name::String) = begin
    if isnothing(result) || ismissing(result)
        return false
    end
    return haskey(result.info, name)
end


"""
        getSignalInfo(instantiatedModel::Modia.InstantiatedModel, name::String)

Returns signal info of variables stored in the result and of parameters.
"""
function SignalTables.getSignalInfo(m::InstantiatedModel, name::String)
    result = m.result
    if haskey(result.info, name)
        # name is a result variable (time-varying variable stored in m.result)
        signalInfo = getSignalInfo(m.result, name)
    elseif name == "_attributes"
        signalInfo = getSignal(m, name)
    else
        # name might be a parameter
        sigValue = get_value(m.evaluatedParameters, name)
        if ismissing(sigValue)
            error("getSignalInfo(.., \"$name\"): signal name is not known")
        end
        sigUnit = unitAsParseableString(sigValue)
        if sigUnit == ""
            signalInfo = SignalTables.Par(_eltypeOrType = eltypeOrType(sigValue))
        else
            signalInfo = SignalTables.Par(_eltypeOrType = eltypeOrType( ustrip.(sigValue) ), unit=sigUnit)
        end
        _size = nothing
        size_available = false
        try
            _size = size(sigValue)
            size_available = true
        catch
            size_available = false
        end
        if size_available
            signalInfo[:_size] = _size
        end
    end
    return signalInfo
end


#=

"""
        getSignalInfo(result::Modia.result, name::String)

Returns signal info of variables stored in result.
"""
function SignalTables.getSignalInfo(result::Result, name::String)
    resInfo = result.info[name]
    if haskey(resInfo.signal, :values)
        signalInfo = copy(resInfo.signal)
        delete!(signalInfo, :values)
        signalValues = resInfo.signal[:values]
        signalInfo[:_eltypeOrType] = eltypeOrType(signalValues)
        signalInfo[:_size]     = size(signalValues)
        return signalInfo
    end
end
=#

#=

"""
        getSignalInfo(instantiatedModel::Modia.InstantiatedModel|result::Modia.Result, name::String)
"""
function getSignalInfo(instantiatedModel::Modia.InstantiatedModel|result::Modia.Result, name::String)
end
function getSignalInfo(instantiatedModel::Modia.InstantiatedModel|result::Modia.Result, name::String)
end
=#


function get_algorithmName_for_heading(m::InstantiatedModel)::String
    if ismissing(m.algorithmName)
        algorithmName = "???"
    else
        algorithmName = m.algorithmName
        i1 = findfirst("CompositeAlgorithm", algorithmName)
        if !isnothing(i1)
            i2 = findfirst("Vern" , algorithmName)
            i3 = findfirst("Rodas", algorithmName)
            success = false
            if !isnothing(i2) && !isnothing(i3)
                i2b = findnext(',', algorithmName, i2[1])
                i3b = findnext('{', algorithmName, i3[1])
                if !isnothing(i2b) && !isnothing(i3b)
                    algorithmName = algorithmName[i2[1]:i2b[1]-1] * "(" * algorithmName[i3[1]:i3b[1]-1] * "())"
                    success = true
                end
            end
            if !success
                algorithmName = "CompositeAlgorithm"
            end
        else
            i1 = findfirst('{', algorithmName)
            if !isnothing(i1)
                algorithmName = algorithmName[1:i1-1]
            end
            i1 = findlast('.', algorithmName)
            if !isnothing(i1)
                algorithmName = algorithmName[i1+1:end]
            end
        end
    end
    return algorithmName
end


"""
    leaveName = get_leaveName(pathName::String)

Return the `leaveName` of `pathName`.
"""
get_leaveName(pathName::String) =
    begin
        j = findlast('.', pathName);
        typeof(j) == Nothing || j >= length(pathName) ? pathName : pathName[j+1:end]
    end


"""
    SignalTables.getDefaultHeading(instantiatedModel::Modia.InstantiatedModel)

Return default heading of instantiatedModel as a string.
"""
function SignalTables.getDefaultHeading(m::InstantiatedModel{FloatType,TimeType}) where {FloatType,TimeType}
    if isnothing(m) || ismissing(m) || ismissing(m.result)
        return ""
    end
    algorithmName = get_algorithmName_for_heading(m)
    if FloatType == "Float64"
        heading = m.modelName * " (" * algorithmName * ")"
    else
        heading = m.modelName * " (" * algorithmName * ", " * string(FloatType) * ")"
    end
end


# For backwards compatibility

"""
    signal    = get_result(instantiatedModel, name; unit=true)
    dataFrame = get_result(instantiatedModel; onlyStates=false, extraNames=missing)

- First form: After a successful simulation of `instantiatedModel`, return
  the result for the signal `name::String` as vector of points
  together with its unit. The time vector has path name `"time"`.
  If `unit=false`, the signal is returned, **without unit**.

- Second form: Return the **complete result** in form of a DataFrame object.
  Therefore, the whole functionality of package [DataFrames](https://dataframes.juliadata.org/stable/)
  can be used, including storing the result on file in different formats.
  Only scalar Var(..) variables are included in the DataFrame object.
  If `onlyStates=true`, then only the states and the signals
  identified with `extraNames::Vector{String}` are stored in `dataFrame`.
  If `onlyStates=false` and `extraNames` given, then only the signals
  identified with `extraNames` are stored in `dataFrame`.

# Example

```julia
using Modia
@usingModiaPlot
using Unitful

include("\$(Modia.path)/examples/Pendulum.jl")
using  .Model_Pendulum

pendulum = simulationModel(Pendulum)
simulate!(pendulum, stopTime=7.0)

# Get one signal from the result and plot with the desired plot package
time = get_result(pendulum, "time")  # vector with unit u"s"
phi  = get_result(pendulum, "phi")   # vector with unit u"rad"

import PyPlot
PyPlot.figure(4)   # Change to figure 4 (or create it, if it does not exist)
PyPlot.clf()       # Clear current figure
PyPlot.plot(stripUnit(time), stripUnit(phi), "b--", label="phi in " * string(unit(phi[1])))
PyPlot.xlabel("time in " * string(unit(time[1])))
PyPlot.legend()

# Get complete result and plot one signal
result = get_result(pendulum)
plot(result, "phi")

# Get only states to be used as reference and compare result with reference
reference = get_result(pendulum, onlyStates=true)
(success, diff, diff_names, max_error, within_tolerance) =
    SignalTables.compareResults(result, reference, tolerance=0.01)
println("Check results: success = $success")
```
"""
get_result(m::InstantiatedModel, name::AbstractString; unit=true) = unit ? getValuesWithUnit(m,name) : getValues(m,name)

function get_result(m::InstantiatedModel; onlyStates=false, extraNames=missing)
    dataFrame = DataFrames.DataFrame()
    dataFrame[!,"time"] = getValues(m, "time")

    if onlyStates || !ismissing(extraNames)
        if onlyStates
            for name in getStateNames(m)
                dataFrame[!,name] = getValues(m, name)
            end
        end
        if !ismissing(extraNames)
            for name in extraNames
                dataFrame[!,name] = getValues(m, name)
            end
        end

    else
        for name in getSignalNames(m, getPar=false, getMap=false)
            if name != "time"
                dataFrame[!,name] = getValues(m,name)
            end
        end
    end
    return dataFrame
end

timeSignalName(  m::InstantiatedModel) = "time"
hasOneTimeSignal(m::InstantiatedModel) = true
signalNames(m::InstantiatedModel) = sort!( getSignalNames(m) )
printResultInfo(m::InstantiatedModel) = showInfo(m)

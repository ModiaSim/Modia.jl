# License for this file: MIT (expat)
# Copyright 2022, DLR Institute of System Dynamics and Control

#--------------------------------------------------------------------------------------------------
#        Provide the overloaded Abstract Signal Tables Interface for the results of SimulationModel
#--------------------------------------------------------------------------------------------------

@inline function checkMissingResult(m::SimulationModel, name::String)::Bool
    if isnothing(m) || ismissing(m) || ismissing(m.result)
        error("$name: No simulation results available.")
    end
    return true
end

SignalTables.isSignalTable(r::Result) = true
SignalTables.isSignalTable(m::SimulationModel) = true


"""
    independentSignalNames(instantiatedModel::Modia.SimulationModel|result::Modia.Result)::Vector{String}
    
Return the name of the independent variable of the result stored in instantiatedModel or in result.
"""
SignalTables.independentSignalNames(result::Result)     = [result.timeName]
SignalTables.independentSignalNames(m::SimulationModel) = begin
    if ismissing(m.result)
        error("independentSignalName(..): No simulation results available in instantiated model of $(m.modelName)")
    end
    SignalTables.independentSignalNames(m.result)
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


"""
    signalNames(instantiatedModel::Modia.SimulationModel; var=true, par=true)::Vector{String}

Returns a string vector of the variables of an
[`@instantiateModel`](@ref) that are present in the result or are (evaluated) parameters.
- If var=true, Var(..) variables are included.
- If par=true, Par(..) variables are included.
"""
function SignalTables.signalNames(m::SimulationModel; var=true, par=true)::Vector{String}
    if ismissing(m.result)
        error("signalNames(..): No simulation results available in instantiated model of $(m.modelName)")
    end
    names1 = collect(keys(m.result.info))
    if var && !par
        return names1
    end
    names2 = setdiff(getParameterNames(m.evaluatedParameters), m.hideResult_names) # parameters without hideResult_names
    if var && par
        return union(names1,names2)
    elseif par
        return setdiff(names2,names1)
    else
        return String[]
    end
end


"""
    signalNames(result::Modia.Result)::Vector{String}

Returns a string vector of the variables that are present in the result.
"""
SignalTables.signalNames(result::Result) = collect(keys(result.info))


"""
    getIndependentSignalsSize(instantiatedModel::Modia.SimulationModel|result::Modia.Result)::Dims

Returns the lengths of the independent signal of the result as (len,).
"""
SignalTables.getIndependentSignalsSize(m::SimulationModel) = SignalTables.getIndependentSignalsSize(m.result)
SignalTables.getIndependentSignalsSize(result::Result) = (sum(length(sk) for sk in result.t), )


"""
    getSignal(instantiatedModel::Modia.SimulationModel|result::Modia.Result, name::String)

Returns signal `name` of the result present in instantiatedModel
(that is a [`SignalTables.Var`](@ref) or a [`SignalTables.Par`](@ref)).
If `name` does not exist, an error is raised.
"""
function SignalTables.getSignal(m::SimulationModel, name::String)
    checkMissingResult(m, name)
    result = m.result
    
    if haskey(result.info, name)
        # name is a result variable (time-varying variable stored in m.result)
        signal = getSignal(result, name)
    else
        # name might be a parameter
        sigValue = get_value(m.evaluatedParameters, name)
        if ismissing(sigValue)
            error("getSignal(.., $name): name is not known")
        end
        sigUnit = unitAsParseableString(sigValue)   
        if sigUnit == ""
            signal = Par(value = sigValue)
        else
            signal = Par(value = ustrip.(sigValue), unit=sigUnit)
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
                w_value = ustrip.(w_value)
            end

        # w_invariant is defined in all segments and has no size information defined - add size information
            resInfo.id[1] = ValuesID(index, size(w_value))

        resInfo._basetype = basetype(w_value)            
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
    hasSignal(instantiatedModel::Modia.SimulationModel, name::String)

Returns `true` if signal `name` is present in the instantiatedModel result or in the evaluated parameters.
"""
SignalTables.hasSignal(m::SimulationModel, name::String) = begin
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
        getSignalInfo(instantiatedModel::Modia.SimulationModel, name::String)
        
Returns signal info of variables stored in the result and of parameters.
"""
function SignalTables.getSignalInfo(m::SimulationModel, name::String)
    result = m.result    
    if haskey(result.info, name)
        # name is a result variable (time-varying variable stored in m.result)
        signalInfo = getSignalInfo(m.result, name)
    else
        # name might be a parameter
        sigValue = get_value(m.evaluatedParameters, name)
        if ismissing(sigValue)
            error("getSignalInfo(.., $name): name is not known")
        end
        sigUnit = unitAsParseableString(sigValue)   
        if sigUnit == ""
            signalInfo = Par(_basetype = basetype(sigValue))
        else
            signalInfo = Par(_basetype = basetype( ustrip.(sigValue) ), unit=sigUnit)
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
        signalInfo[:_basetype] = basetype(signalValues)
        signalInfo[:_size]     = size(signalValues)
        return signalInfo
    end
end
=#

#=

"""
        getSignalInfo(instantiatedModel::Modia.SimulationModel|result::Modia.Result, name::String)
"""
function getSignalInfo(instantiatedModel::Modia.SimulationModel|result::Modia.Result, name::String)
end
function getSignalInfo(instantiatedModel::Modia.SimulationModel|result::Modia.Result, name::String)
end
=#


function get_algorithmName_for_heading(m::SimulationModel)::String
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
    SignalTables.getDefaultHeading(instantiatedModel::Modia.SimulationModel)
    
Return default heading of instantiatedModel as a string.
""" 
function SignalTables.getDefaultHeading(m::SimulationModel{FloatType,TimeType}) where {FloatType,TimeType}
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
  Furthermore, also plot can be used on dataFrame.
  Parameters and zero-value variables are stored as SignalTables.OneValueVector inside dataFrame
  (are treated as vectors, but actually only the value and the number
  of time points is stored). If `onlyStates=true`, then only the states and the signals
  identified with `extraNames::Vector{String}` are stored in `dataFrame`.
  If `onlyStates=false` and `extraNames` given, then only the signals
  identified with `extraNames` are stored in `dataFrame`.
  These keyword arguments are useful, if `dataFrame` shall be
  utilized as reference result used in compareResults(..).

In both cases, a **view** on the internal result memory is provided
(so result data is not copied).

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
function get_result(m::SimulationModel, name::AbstractString; unit=true)
    #(xsig, xsigLegend, ysig, ysigLegend, yIsConstant) = SignalTables.getPlotSignal(m, "time", name)

    #resIndex = m.variables[name]
    #ysig = ResultView(m.result, abs(resIndex), resIndex < 0)

    #if SignalTables.timeSignalName(m) != 1
    if length(m.result.t) > 1
        error("Error in Modia.get_result(\"$name\"), because function cannot be used for a segmented simulation with more as one segmented.")
    end

    (tsig2, ysig2, ysigType) = SignalTables.rawSignal(m, name)
    ysig = ysig2[1]
    ysig = unit ? ysig : stripUnit.(ysig)

    #=
    if yIsConstant
        if ndims(ysig) == 1
            ysig = fill(ysig[1], length(xsig))
        else
            ysig = fill(ysig[1,:], length(xsig))
        end
    end
    =#

    return ysig
end


function setEvaluatedParametersInDataFrame!(obj::OrderedDict{Symbol,Any}, result_info, dataFrame::DataFrames.DataFrame, path::String, nResult::Int)::Nothing
    for (key,value) in zip(keys(obj), obj)
        name = appendName(path, key)
        if typeof(value) <: OrderedDict{Symbol,Any}
            setEvaluatedParametersInDataFrame!(value, result_info, dataFrame, name, nResult)
        elseif !haskey(result_info, name)
            dataFrame[!,name] = SignalTables.OneValueVector(value,nResult)
        end
    end
    return nothing
end


function get_result(m::SimulationModel; onlyStates=false, extraNames=missing)
    error("get_result(instantiatedModel) is no longer supported")
    
    #=
    if length(m.result.t) > 1
        error("Error in Modia.get_result(...), because function cannot be used for a segmented simulation with more as one segmented.")
    end

    dataFrame = DataFrames.DataFrame()

    (timeSignal, signal, signalType) = SignalTables.rawSignal(m, "time")
    dataFrame[!,"time"] = timeSignal[1]

    if onlyStates || !ismissing(extraNames)
        if onlyStates
            for name in keys(m.equationInfo.x_dict)
                (timeSignal, signal, signalType) = SignalTables.rawSignal(m, name)
                dataFrame[!,name] = signal[1]
            end
        end
        if !ismissing(extraNames)
            for name in extraNames
                (timeSignal, signal, signalType) = SignalTables.rawSignal(m, name)
                dataFrame[!,name] = signal[1]
            end
        end

    else
        for name in keys(m.result.info)
            if name != "time"
                (timeSignal, signal, signalType) = SignalTables.rawSignal(m, name)
                dataFrame[!,name] = signal[1]
            end
        end

        setEvaluatedParametersInDataFrame!(m.evaluatedParameters, m.result.info, dataFrame, "", length(timeSignal[1]))
    end
    return dataFrame
    =#
end
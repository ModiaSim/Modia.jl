# License for this file: MIT (expat)
# Copyright 2022, DLR Institute of System Dynamics and Control

#------------------------------------------------------------------------------------------------
#        Provide the overloaded ModiaResult Abstract Interface for the results of SimulationModel
#------------------------------------------------------------------------------------------------

"""
    hasSignal(instantiatedModel::Modia.SimulationModel, name::AbstractString)

Return true if time-varying variable `name` (for example `name = "a.b.c"`)
is defined in the instantiateModel that can be accessed and can be used for plotting.
"""
ModiaResult.hasSignal(m::SimulationModel, name::AbstractString) = begin
    if isnothing(m) || ismissing(m) || ismissing(m.result)
        return false
    end
    haskey(m.result.info, name) #|| !ismissing(get_value(m.evaluatedParameters, name))
end


"""
    timeSignalName(instantiatedModel::Modia.SimulationModel)
    
Return the name of the independent variable of the result stored in instantiatedModel.
"""
ModiaResult.timeSignalName(m::SimulationModel) = m.result.timeName


"""
    names = signalNames(instantiatedModel::Modia.SimulationModel)

Return the names of the time-varying variables of an
[`@instantiateModel`](@ref) that are present in the result (e.g. can be accessed for plotting).
"""
ModiaResult.signalNames(m::SimulationModel) = collect(keys(m.result.info))


function getVariableType(result::Result, resInfo::ResultInfo)::DataType
    if ismissing(resInfo.VariableType)
        @assert(resInfo.kind == RESULT_W_INVARIANT)
        return typeof( result.w_invariant[1][1][resInfo.id[1].index] )
    else
        return resInfo.VariableType
    end
end


"""
    info = SignalInfo(instantiatedModel::SimulationModel, name::AbstractString)

Return information about signal `name` of the result present in instantiatedModel.
"""
function ModiaResult.SignalInfo(m::SimulationModel, name::AbstractString)
    resInfo = m.result.info[name]
    if resInfo.kind == RESULT_ELIMINATED
        aliasResInfo = m.result.info[resInfo.aliasName]
        VariableType = getVariableType(m.result, aliasResInfo)
        return ModiaResult.SignalInfo(ModiaResult.Eliminated, VariableType, aliasResInfo.unit, aliasResInfo.value, resInfo.aliasName, resInfo.aliasNegate)
    end

    kind = resInfo.kind == RESULT_T        ? ModiaResult.Independent :
          (resInfo.kind == RESULT_CONSTANT ? ModiaResult.Constant    :
          ( (resInfo.kind == RESULT_X || resInfo.kind == RESULT_DER_X) && isInvariant(resInfo.id) || (resInfo.W_INVARIANT ? ModiaResult.Invariant : ModiaResult.Segmenteded)))

    VariableType = getVariableType(m.result, resInfo)
    return ModiaResult.SignalInfo(kind, VariableType, resInfo.unit, resInfo.value, "", false)
end


"""
    s = signalValues(instantiatedModel::Modia.SimulationModel, name; unitless=false)
    
Return signal `name` of the result present in instantiatedModel.
If `unitless=false`, the signal is returned with its unit.
"""
function ModiaResult.signalValues(m::SimulationModel, name; unitless=false)
    if haskey(m.result.info, name)
        return signalResultValues(m.result, name; unitless=unitless)
    else
        value = get_value(m.evaluatedParameters, name)
        if ismissing(value)
            error("signal(..): \"$name\" not in result and not a parameter of model $(m.modelName))")
        end
        return ModiaResult.OneValueVector(value, sum(length(tk) for tk in m.result.t))
    end
end


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
    ModiaResult.defaultHeading(instantiatedModel::Modia.SimulationModel)
    
Return default heading of instantiatedModel result as a string 
(can be used as default heading for a plot).
""" 
function ModiaResult.defaultHeading(m::SimulationModel)
    FloatType = get_leaveName( string( typeof( m.x_start[1] ) ) )

    algorithmName = get_algorithmName_for_heading(m)
    if FloatType == "Float64"
        heading = m.modelName * " (" * algorithmName * ")"
    else
        heading = m.modelName * " (" * algorithmName * ", " * FloatType * ")"
    end
    return heading
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
  Parameters and zero-value variables are stored as ModiaResult.OneValueVector inside dataFrame
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
    ModiaResult.compareResults(result, reference, tolerance=0.01)
println("Check results: success = $success")
```
"""
function get_result(m::SimulationModel, name::AbstractString; unit=true)
    #(xsig, xsigLegend, ysig, ysigLegend, yIsConstant) = ModiaResult.getPlotSignal(m, "time", name)

    #resIndex = m.variables[name]
    #ysig = ResultView(m.result, abs(resIndex), resIndex < 0)

    #if ModiaResult.timeSignalName(m) != 1
    if length(m.result.t) > 1
        error("Error in Modia.get_result(\"$name\"), because function cannot be used for a segmenteded simulation with more as one segmented.")
    end

    (tsig2, ysig2, ysigType) = ModiaResult.rawSignal(m, name)
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
            dataFrame[!,name] = ModiaResult.OneValueVector(value,nResult)
        end
    end
    return nothing
end


function get_result(m::SimulationModel; onlyStates=false, extraNames=missing)
    if length(m.result.t) > 1
        error("Error in Modia.get_result(...), because function cannot be used for a segmenteded simulation with more as one segmented.")
    end

    dataFrame = DataFrames.DataFrame()

    (timeSignal, signal, signalType) = ModiaResult.rawSignal(m, "time")
    dataFrame[!,"time"] = timeSignal[1]

    if onlyStates || !ismissing(extraNames)
        if onlyStates
            for name in keys(m.equationInfo.x_dict)
                (timeSignal, signal, signalType) = ModiaResult.rawSignal(m, name)
                dataFrame[!,name] = signal[1]
            end
        end
        if !ismissing(extraNames)
            for name in extraNames
                (timeSignal, signal, signalType) = ModiaResult.rawSignal(m, name)
                dataFrame[!,name] = signal[1]
            end
        end

    else
        for name in keys(m.result.info)
            if name != "time"
                (timeSignal, signal, signalType) = ModiaResult.rawSignal(m, name)
                dataFrame[!,name] = signal[1]
            end
        end

        setEvaluatedParametersInDataFrame!(m.evaluatedParameters, m.result.info, dataFrame, "", length(timeSignal[1]))
    end
    return dataFrame
end
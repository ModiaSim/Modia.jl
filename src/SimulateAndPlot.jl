
export simulate!, linearize!, get_result

using  ModiaPlot
import ModiaBase
using  Measurements
import MonteCarloMeasurements
using  Unitful
using  Test
import DifferentialEquations
import DataFrames
import ForwardDiff
import FiniteDiff

#---------------------------------------------------------------------
#                          Simulation
#---------------------------------------------------------------------

                   
"""
    simulate!(instantiatedModel [, algorithm]; merge = nothing,
              tolerance = 1e-6, startTime = 0.0, stopTime = 1.0, interval = NaN,
              interp_points = 0, adaptive = true, log = false, logStates = false,
              logEvents = false, logParameters = false, logEvaluatedParameters = false, 
              requiredFinalStates = nothing)

Simulate `instantiatedModel::SimulationModel` with `algorithm`
(= `alg` of [ODE Solvers of DifferentialEquations.jl](https://diffeq.sciml.ai/stable/solvers/ode_solve/)).
If the `algorithm` argument is missing, a default algorithm will be chosen from DifferentialEquations
(for details see [https://arxiv.org/pdf/1807.06430](https://arxiv.org/pdf/1807.06430), Figure 3).

The simulation results stored in `model` can be plotted with ModiaPlot.plot and the result values
can be retrieved with [`get_result`](@ref).


# Optional Arguments

- `merge`: Define parameters and init/start values that shall be merged with the previous values
           stored in `model`, before simulation is started.  
- `tolerance`: Relative tolerance.
- `startTime`: Start time. If value is without unit, it is assumed to have unit [s].
- `stopTime`: Stop time. If value is without unit, it is assumed to have unit [s].
- `interval`: Interval to store result. If `interval=NaN`, it is internally selected as
              (stopTime-startTime)/500.
              If value is without unit, it is assumed to have unit [s].
- `interp_points`: If crossing functions defined, number of additional interpolation points
              in one step.
- `adaptive`: = true, if the `algorithm` should use step-size control (if available).
              = false, if the `algorithm` should use a fixed step-size of `interval` (if available).
- `log`: = true, to log the simulation.
- `logStates`: = true, to log the states, its init/start values and its units.
- `logEvents`: = true, to log events.
- `logParameters`: = true, to log parameters and init/start values defined in model.
- `logEvaluatedParameters`: = true, to log the evaluated parameter and init/start values.
- `requiredFinalStates`: is not `nothing`: Check whether the ODE state vector at the final time instant with `@test`
              is in agreement to vector `requiredFinalStates` with respect to some tolerance. If this is not the case, print the
              final state vector (so that it can be included with copy-and-paste in the simulate!(..) call).

# Examples

```julia
using ModiaPlot
using DifferentialEquations

# Define model
inputSignal(t) = sin(t)

FirstOrder = Model(
    T = 0.2,
    x = Var(init=0.3),
    equations = :[u = inputSignal(time/u"s"),
                  T * der(x) + x = u,
                  y = 2*x]
)

# Modify parameters and initial values of model
FirstOrder2 = FirstOrder | Map(T = 0.4, x = Var(init=0.6))

# Instantiate model
firstOrder = @instantiateModel(FirstOrder2, logCode=true)


# Simulate with automatically selected algorithm and 
# modified parameter and initial values
simulate!(firstOrder, stopTime = 1.0, merge = Map(T = 0.6, x = 0.9), logEvaluatedParameters=true)

# Plot variables "x", "u" in diagram 1, "der(x)" in diagram 2, both diagrams in figure 3
plot(firstOrder, [("x","u"), "der(x)"], figure=3)

# Retrieve "time" and "u" values:
get_result(firstOrder, "time")
get_result(firstOrder, "u")
    
    
# Simulate with Runge-Kutta 5/4 with step-size control
simulate!(firstOrder, Tsit5(), stopTime = 1.0)

# Simulate with Runge-Kutta 4 with fixed step size
simulate!(firstOrder, RK4(), stopTime = 1.0, adaptive=false)

# Simulate with algorithm that switches between 
# Verners Runge-Kutta 6/5 algorithm if non-stiff region and
# Rosenbrock 4 (= A-stable method) if stiff region with step-size control
simulate!(firstOrder, AutoVern6(Rodas4()), stopTime = 1.0)

# Simulate with Sundials CVODE (BDF method with variable order 1-5) with step-size control
using Sundials
simulate!(firstOrder, CVODE_BDF(), stopTime = 1.0)
```
"""
function simulate!(m::Nothing, args...; kwargs...)
    @info "The call of simulate!(..) is ignored, since the first argument is nothing."
    return nothing
end
function simulate!(m::SimulationModel{FloatType,ParType,EvaluatedParType,TimeType}, algorithm=missing; merge=nothing, kwargs...) where {FloatType,TimeType,ParType,EvaluatedParType}
        empty!(m.result)
        options = SimulationOptions{FloatType,TimeType}(merge; kwargs...)
        if isnothing(options)
            return nothing
        end
        m.options = options

        # Initialize/re-initialize SimulationModel
        if m.options.log || m.options.logEvaluatedParameters || m.options.logStates
            println("... Simulate model ", m.modelName)
        end
        
        cpuStart::UInt64 = time_ns()
        cpuLast::UInt64  = cpuStart
        cpuStartIntegration::UInt64 = cpuStart
        m.algorithmType = typeof(algorithm)       
        success = init!(m)
        if !success
            return nothing
        end
            
        if m.eventHandler.restart == Terminate
            finalStates = m.x_init
            terminate!(m, finalStates, m.options.startTime)

        else
            # Define problem and callbacks based on algorithm and model type
            if abs(m.options.interval) < abs(m.options.stopTime-m.options.startTime)
                tspan2 = m.options.startTime:m.options.interval:m.options.stopTime
            else
                tspan2 = [m.options.startTime, m.options.stopTime]
            end
           	tspan   = (m.options.startTime, m.options.stopTime)            
            abstol  = 0.1*m.options.tolerance
            problem = DifferentialEquations.ODEProblem(derivatives!, m.x_init, tspan, m) 
            
            callback2 = DifferentialEquations.DiscreteCallback(timeEventCondition!, affectTimeEvent!)   
            eh = m.eventHandler            
            if eh.nz > 0        
                # Due to DifferentialEquations bug https://github.com/SciML/DifferentialEquations.jl/issues/686
                # FunctionalCallingCallback(outputs!, ...) is not correctly called when zero crossings are present.
                # A temporary fix is to use time events at the communication points, but this slows down simulation.
                callback1 = DifferentialEquations.PresetTimeCallback(tspan2, affect_outputs!)  
                callback3 = DifferentialEquations.VectorContinuousCallback(stateEventCondition!, 
                                 affectStateEvent!, eh.nz, interp_points=m.options.interp_points)                       
                callbacks = DifferentialEquations.CallbackSet(callback1, callback2, callback3)
            else
                callback1 = DifferentialEquations.FunctionCallingCallback(outputs!, funcat=tspan2, tdir=sign(m.options.stopTime-m.options.startTime))
                callbacks = DifferentialEquations.CallbackSet(callback1, callback2)
            end
    
            # Initial step size (the default of DifferentialEquations is too large) + step-size of fixed-step algorithm
            dt = m.options.adaptive ? m.options.interval/10 : m.options.interval   # initial step-size
    
            # Compute solution 
            tstops = (m.eventHandler.nextEventTime,)
            solution = ismissing(algorithm) ? DifferentialEquations.solve(problem, reltol=m.options.tolerance, abstol=abstol,
                                                save_everystep=false, save_start=false, save_end=true,
                                                callback=callbacks, adaptive=m.options.adaptive, dt=dt, tstops = tstops) :
                                            DifferentialEquations.solve(problem, algorithm, reltol=m.options.tolerance, abstol=abstol,
                                                save_everystep=false, save_start=false, save_end=true,
                                                callback=callbacks, adaptive=m.options.adaptive, dt=dt, tstops = tstops)
            m.solution = solution
            if ismissing(algorithm)
                m.algorithmType = typeof(solution.alg)
            end
    
            # Terminate simulation
            finalStates = solution[:,end]
            terminate!(m, finalStates, solution.t[end])
        end
        
        if m.options.log
            cpuTimeInitialization = convert(Float64, (cpuStartIntegration - cpuStart) * 1e-9)
            cpuTimeIntegration    = convert(Float64, (time_ns() - cpuStartIntegration) * 1e-9)
            cpuTime               = cpuTimeInitialization + cpuTimeIntegration

            println("      Termination at time = ", solution.t[end], " s")
            println("        cpuTime         = ", round(cpuTime, sigdigits=3), " s")
            #println("        cpuTime         = ", round(cpuTime              , sigdigits=3), " s (init: ",
            #                                      round(cpuTimeInitialization, sigdigits=3), " s, integration: ",
            #                                      round(cpuTimeIntegration   , sigdigits=3), " s)")
            println("        algorithm       = ", m.algorithmType == Nothing ? Nothing : solution.alg)
            println("        FloatType       = ", FloatType)
            println("        interval        = ", m.options.interval, " s")
            println("        tolerance       = ", m.options.tolerance, " (relative tolerance)")
            println("        nEquations      = ", length(m.x_start))
            println("        nResults        = ", length(m.result))
            println("        nAcceptedSteps  = ", solution.destats.naccept)
            println("        nRejectedSteps  = ", solution.destats.nreject)
            println("        nGetDerivatives = ", m.nGetDerivatives, " (total number of getDerivatives! calls)")
            println("        nf              = ", solution.destats.nf, " (number of getDerivatives! calls from integrator)")
            println("        nZeroCrossings  = ", eh.nZeroCrossings, " (number of getDerivatives! calls for zero crossing detection)")
            println("        nJac            = ", solution.destats.njacs, " (number of Jacobian computations)")
            println("        nErrTestFails   = ", solution.destats.nreject)          
            println("        nTimeEvents     = ", eh.nTimeEvents)
            println("        nStateEvents    = ", eh.nStateEvents)
            println("        nRestartEvents  = ", eh.nRestartEvents)      
        end

        requiredFinalStates = m.options.requiredFinalStates
        if !isnothing(requiredFinalStates)
            if length(finalStates) != length(requiredFinalStates)
                success = false
            else
                success = finalStates == requiredFinalStates || isapprox(finalStates, requiredFinalStates, rtol=1e-3)
            end

            if success
                @test success
            else
                if length(requiredFinalStates) > 0 && typeof(requiredFinalStates[1]) <: Measurements.Measurement
                    println(  "\nrequiredFinalStates = ", measurementToString(requiredFinalStates))
                    printstyled("finalStates         = ", measurementToString(finalStates), "\n\n", bold=true, color=:red)
                else
                    println(  "\nrequiredFinalStates = ", requiredFinalStates)
                    printstyled("finalStates         = ", finalStates, "\n\n", bold=true, color=:red)
                end
                @test finalStates == requiredFinalStates  #|| isapprox(finalStates, requiredFinalStates, rtol=1e-3)
            end
        end
        return solution
end

#get_x_startIndexAndLength(m::SimulationModel, name) = ModiaBase.get_x_startIndexAndLength(m.equationInfo, name)


#---------------------------------------------------------------------
#                          Linearization
#---------------------------------------------------------------------

"""
    (A, finalStates) = linearize!(instantiatedModel [, algorithm];
                                  stopTime = 0.0,
                                  analytic = false,
                                  <all other keyword arguments of simulate!>)
    
Simulate until `stopTime` and linearize `instantiatedModel` at `finalStates`.
The names of the state vector can be inquired by `get_xNames(instantiatedModel)`.

By default, linearization is performed numerically with a central finite difference
approximation using package [FiniteDiff](https://github.com/JuliaDiff/FiniteDiff.jl).
When setting `analytic = true`, linearization is preformed analytically with package 
[ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl),
so is computed by symbolically differentiating the model.
`ForwardDiff` might not be compatible with some floating point types, such as
`Measurements` and Julia triggers an error that some overloaded
operations are ambiguous. So `analytic=true` will not work in such cases.

Analytic linearization returns matrix `A` in full precision whereas numeric linearization
returns `A` in reduced precision (if FloatType = Float64, analytic linearization results in about
15 correct digits and numeric linearization in about 10 correct digits in the result).
You can improve this situation, by using a larger
`FloatType` for `instantiatedModel`, in case this is critical (see example below).

# Output arguments

- `A::Matrix`: Matrix A of the linear ODE: ``\\dot{\\Delta x} = A*\\Delta x``.

- `finalStates::Vector`: Linearization point.


# Example

```julia
using TinyModia
using DoubleFloats
using Measurements

FirstOrder = Model(
    T = 0.4 ± 0.04,
    x = Var(init = 0.9 ± 0.09),
    equations = :[u = inputSignal(time/u"s"),
                  T * der(x) + x = u]
)

firstOrder1 = @instantiateModel(FirstOrder, FloatType = Measurement{Float64})

# Standard precision
(A1, finalStates1) = linearize!(firstOrder1)

# Higher precision
firstOrder2 = SimulationModel{Measurement{Double64}}(firstOrder1)
(A2, finalStates2) = linearize!(firstOrder2)

# Show results with 15 digits (default print with Measurements shows 3 digits)
println(IOContext(stdout, :error_digits=>15), "A1 = ", A1)
println(IOContext(stdout, :error_digits=>15), "A2 = ", A2)
```
"""
function linearize!(m::Nothing, args...; kwargs...)
    @info "The call of linearize!(..) is ignored, since the first argument is nothing."
    return   nothing
end
function linearize!(m::SimulationModel{FloatType,ParType,EvaluatedParType,TimeType}, 
                    algorithm=missing;
                    merge = nothing, stopTime = 0.0, analytic = false, kwargs...) where {FloatType,ParType,EvaluatedParType,TimeType}
    solution = simulate!(m, algorithm; merge=merge, stopTime=stopTime, kwargs...)
    finalStates = solution[:,end]
    
    # Function that shall be linearized
    function modelToLinearize!(der_x, x)
        Base.invokelatest(m.getDerivatives!, der_x, x, m, m.options.startTime)
        return nothing
    end
    
    # Linearize
    if analytic
        A = ForwardDiff.jacobian(modelToLinearize!, m.der_x, finalStates)
    else
        A = zeros(FloatType, length(finalStates), length(finalStates))
        FiniteDiff.finite_difference_jacobian!(A, modelToLinearize!, finalStates)
    end
    
    return (A, finalStates)
end




#---------------------------------------------------------------------
#        Provide ModiaPlot result access functions for SimulationModel
#---------------------------------------------------------------------

"""
    ModiaPlot.hasSignal(model::SimulationModel, name::String)
    
Return true if parameter or time-varying variable `name` (for example `a.b.c`)
is defined in the TinyModia SimulationModel (generated with [`TinyModia.@instantiateModel`](@ref)
that can be accessed and can be used for plotting.
"""
ModiaPlot.hasSignal(m::SimulationModel, name) = 
    # m.save_x_in_solution ? name == "time" || haskey(m.equationInfo.x_dict, name) :
     haskey(m.variables, name) || name in m.zeroVariables || !ismissing(get_value(m.evaluatedParameters, name))

        
"""
    ModiaPlot.getNames(model::SimulationModel)
    
Return the variable names (parameters, time-varying variables) of a TinyModia SimulationModel
(generated with [`TinyModia.@instantiateModel`](@ref) that can be accessed
and can be used for plotting.
"""
function ModiaPlot.getNames(m::SimulationModel)
    #if m.save_x_in_solution 
    #    names = ["time"]
    #    append!(names, collect( keys(m.equationInfo.x_dict) ))
    #else
        names = get_names(m.evaluatedParameters)
        append!(names, collect(m.zeroVariables))
        append!(names, collect( keys(m.variables) ) )
    #end
    return sort(names)
end



struct ResultView{T} <: AbstractArray{T, 1}
    v
    index::Int
    neg::Bool
    ResultView{T}(v,i) where {T} = new(v,abs(i),i < 0)
end
ResultView(v,i) = ResultView{typeof(v[1][abs(i)])}(v,i)

Base.getindex(sig::ResultView, i::Int) = sig.neg ? -sig.v[i][sig.index] : sig.v[i][sig.index]
Base.size(sig::ResultView)             = (length(sig.v),)
Base.IndexStyle(::Type{<:ResultView})  = IndexLinear()

#=
import ChainRules

function ChainRules.rrule(::typeof(ResultView), v, i)
    y = ResultView(v,i)
    
    function ResultView_pullback(ȳ)
        return ChainRules.NO_FIELDS, collect(y)...
    end
    
    return y, ResultView_pullback
end
=#


function ModiaPlot.getRawSignal(m::SimulationModel, name)
    if m.save_x_in_solution 
        if name == "time"
            return (false, m.solution.t)
        elseif haskey(m.equationInfo.x_dict, name)
            xe_info = m.equationInfo.x_info[ m.equationInfo.x_dict[name] ]
            @assert(xe_info.length == 1)   # temporarily only scalars are supported
            xe_index = xe_info.startIndex
            return (false, m.solution[xe_index,:])
        end
        
        #else
        #    error("getRawSignal(m, $name): only states can be inquired currently")
        #end
    end
    
    if haskey(m.variables, name)
        resIndex = m.variables[name]
        signal = ResultView(m.result, resIndex)       
        if name == "time" && !(m.options.desiredResultTimeUnit == NoUnits ||
                               m.options.desiredResultTimeUnit == u"s")
            signal  = uconvert.(m.options.desiredResultTimeUnit, signal)
        end
#=        
        negAlias = false
        if resIndex < 0
            resIndex = -resIndex
            negAlias = true
        end

        value  = m.result[1][resIndex]
        signal = Vector{typeof(value)}(undef, length(m.result))
        for (i, value_i) in enumerate(m.result)
            signal[i] = value_i[resIndex]
        end
        if negAlias
            signal .= -signal
        end
=#
        return (false, signal)
        
    elseif name in m.zeroVariables
        return (true, 0.0)

    else
        value = get_value(m.evaluatedParameters, name)
        if ismissing(value)
            error("ModiaPlot.getRawSignal: ", name, " not in result of model ", m.modelName)
        end
        return (true, value)
    end
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


function ModiaPlot.getDefaultHeading(m::SimulationModel)
    FloatType = get_leaveName( string( typeof( m.x_start[1] ) ) )

    if ismissing(m.algorithmType)
        algorithmName = "???"
    else
        algorithmName = string(m.algorithmType)
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

    if FloatType == "Float64"
        heading = m.modelName * " (" * algorithmName * ")"
    else
        heading = m.modelName * " (" * algorithmName * ", " * FloatType * ")"
    end
    return heading
end


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
  Furthermore, also ModiaPlot.plot can be used on dataFrame.
  Parameters and zero-value variables are stored as ModiaPlot.OneValueVector inside dataFrame
  (are treated as vectors, but actually only the value and the number
  of time points is stored). If `onlyStates=true`, then only the states and the signals
  identified with `extraNames::Vector{String}` are stored in `dataFrame`. 
  If `onlyStates=false` and `extraNames` given, then only the signals 
  identified with `extraNames` are stored in `dataFrame`.
  These keyword arguments are useful, if `dataFrame` shall be 
  utilized as reference result used in ModiaPlot.compareResults(..).

In both cases, a **view** on the internal result memory is provided
(so result data is not copied).

# Example

```julia
using TinyModia
using ModiaPlot
using Unitful

include("\$(ModiaBase.path)/demos/models/Model_Pendulum.jl")
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
ModiaPlot.plot(result, "phi")

# Get only states to be used as reference and compare result with reference
reference = get_result(pendulum, onlyStates=true)
(success, diff, diff_names, max_error, within_tolerance) = 
    ModiaPlot.compareResults(result, reference, tolerance=0.01)
println("Check results: success = $success")
```
"""
function get_result(m::SimulationModel, name::String; unit=true)
    #(xsig, xsigLegend, ysig, ysigLegend, yIsConstant) = ModiaPlot.getPlotSignal(m, "time", name)

    #resIndex = m.variables[name]
    #ysig = ResultView(m.result, abs(resIndex), resIndex < 0)
        
    (isConstant, ysig) = ModiaPlot.getRawSignal(m, name)

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


function setEvaluatedParametersInDataFrame!(obj::NamedTuple, variables, dataFrame::DataFrames.DataFrame, path::String, nResult::Int)::Nothing
    for (key,value) in zip(keys(obj), obj)
        name = appendName(path, key)
        if typeof(value) <: NamedTuple
            setEvaluatedParametersInDataFrame!(value, variables, dataFrame, name, nResult)
        elseif !haskey(variables, name)
            dataFrame[!,name] = ModiaPlot.OneValueVector(value,nResult)
        end
    end
    return nothing
end

function get_result(m::SimulationModel; onlyStates=false, extraNames=missing)
    #if m.save_x_in_solution 
    #    @error "get_result(instantiatedModel) not yet supported, if result is stored in DifferentialEquations.jl solution"
    #end
        
    dataFrame = DataFrames.DataFrame()
    
    if onlyStates || !ismissing(extraNames)
        dataFrame[!,"time"] = get_result(m, "time")
        if onlyStates
            for (name, dummy) in m.equationInfo.x_dict
                dataFrame[!,name] = ResultView(m.result, m.variables[name])
            end
        end
        if !ismissing(extraNames)
            for name in extraNames
                dataFrame[!,name] = get_result(m, name)
            end
        end
        
    else
        for (name, resIndex) in m.variables
            if name == "time"
                dataFrame[!,name] = get_result(m, "time")  # Takes care of conversion to unit m.options.desiredResultTimeUnit
            else
                dataFrame[!,name] = ResultView(m.result, resIndex) 
            end
        end
    
        zeroVariable = ModiaPlot.OneValueVector(0.0, length(m.result))
        for name in m.zeroVariables
            dataFrame[!,name] = zeroVariable
        end
    
        setEvaluatedParametersInDataFrame!(m.evaluatedParameters, m.variables, dataFrame, "", length(m.result))
    end
    return dataFrame
end

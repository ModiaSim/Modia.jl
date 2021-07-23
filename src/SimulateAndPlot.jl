
export simulate!, linearize!, get_result
export @usingModiaPlot, usePlotPackage, usePreviousPlotPackage, currentPlotPackage
export resultInfo, printResultInfo, rawSignal, getPlotSignal, defaultHeading
export signalNames, timeSignalName, hasOneTimeSignal, hasSignal
  
# For backwards compatibility
export getNames, hasName

import ModiaResult
import ModiaResult: @usingModiaPlot, usePlotPackage, usePreviousPlotPackage, currentPlotPackage
import ModiaResult: resultInfo, printResultInfo, rawSignal, getPlotSignal, defaultHeading
import ModiaResult: signalNames, timeSignalName, hasOneTimeSignal, hasSignal

import TimerOutputs

import ModiaBase
using  Measurements
import MonteCarloMeasurements
using  Unitful
using  Test
import Sundials
import DifferentialEquations
import DataFrames
import ForwardDiff
import FiniteDiff

const  CVODE_BDF = Sundials.CVODE_BDF
export CVODE_BDF

#---------------------------------------------------------------------
#                          Simulation
#---------------------------------------------------------------------

                   
"""
    simulate!(instantiatedModel [, algorithm]; 
              merge         = missing,  # change parameter/init/start values
              tolerance     = 1e-6,     # relative tolerance
              startTime     = 0.0, 
              stopTime      = 0.0,      # stopTime >= startTime required
              interval      = missing,  # = (stopTime-startTime)/500
              interp_points = 0, 
              dtmax         = missing,  # = 100*interval
              adaptive      = true, 
              log           = false, 
              logStates     = false,
              logEvents     = false,
              logTiming     = false,
              logParameters = false, 
              logEvaluatedParameters   = false, 
              requiredFinalStates      = missing
              requiredFinalStates_rtol = 1e-3)

Simulate `instantiatedModel::SimulationModel` with `algorithm`
(= `alg` of [ODE Solvers of DifferentialEquations.jl](https://diffeq.sciml.ai/stable/solvers/ode_solve/)).

If the `algorithm` argument is missing, `algorithm=Sundials.CVODE_BDF()` is used, provided
instantiatedModel has `FloatType = Float64`. Otherwise, a default algorithm will be chosen from DifferentialEquations
(for details see [https://arxiv.org/pdf/1807.06430](https://arxiv.org/pdf/1807.06430), Figure 3).
The symbol `CVODE_BDF` is exported from TinyModia, so that `simulate!(instantiatedModel, CVODE_BDF(), ...)`
can be used (instead of `import Sundials; simulate!(instantiatedModel, Sundials.CVODE_BDF(), ...)`).

The simulation results stored in `model` can be plotted with plot and the result values
can be retrieved with `rawSignal(..)` or `getPlotSignal(..)`
(for details see chapter [Results and Plotting](@ref)).


# Optional Arguments

- `merge`: Define parameters and init/start values that shall be merged with the previous values
           stored in `model`, before simulation is started.  
- `tolerance`: Relative tolerance.
- `startTime`: Start time. If value is without unit, it is assumed to have unit [s].
- `stopTime`: Stop time. If value is without unit, it is assumed to have unit [s].
- `interval`: Interval to store result. If `interval=missing`, it is internally selected as
              (stopTime-startTime)/500.
              If value is without unit, it is assumed to have unit [s].
- `interp_points`: If crossing functions defined, number of additional interpolation points
              in one step.
- `dtmax`: Maximum step size. If `dtmax==missing`, it is internally set to `100*interval`.
- `adaptive`: = true, if the `algorithm` should use step-size control (if available).
              = false, if the `algorithm` should use a fixed step-size of `interval` (if available).
- `log`: = true, to log the simulation.
- `logStates`: = true, to log the states, its init/start values and its units.
- `logEvents`: = true, to log events.
- `logTiming`: = true, to log the timing with `instantiatedModel.timer` which is an instance
               of [TimerOutputs](https://github.com/KristofferC/TimerOutputs.jl).TimerOutput.
               A user function can include its timing via\\
               `TimerOutputs.@timeit instantiatedModel.timer "My Timing" <statement>`.
- `logParameters`: = true, to log parameters and init/start values defined in model.
- `logEvaluatedParameters`: = true, to log the evaluated parameter and init/start values that
                            are used for initialization and during simulation.
- `requiredFinalStates`: is not `missing`: Test with `@test` whether the ODE state vector at the 
              final time instant is in agreement to vector `requiredFinalStates` with respect 
              to some relative tolerance `requiredFinalStates_rtol`. If this is not the case, print the
              final state vector (so that it can be included with copy-and-paste in the simulate!(..) call).
- `requiredFinalStates_rtol`: Relative tolerance used for `requiredFinalStates`.

# Examples

```julia
using TinyModia
using DifferentialEquations
using Unitful
@usingModiaPlot

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


# Simulate with automatically selected algorithm (Sundials.CVODE_BDF())
# and modified parameter and initial values
simulate!(firstOrder, stopTime = 1.0, merge = Map(T = 0.6, x = 0.9), logEvaluatedParameters=true)

# Plot variables "x", "u" in diagram 1, "der(x)" in diagram 2, both diagrams in figure 3
plot(firstOrder, [("x","u"), "der(x)"], figure=3)

# Retrieve "time" and "u" values:
usig = getPlotSignal(firstOrder, "x")
       # usig.xsig      : time vector
       # usig.xsigLegend: legend for time vector
       # usig.ysig      : "x" vector
       # usig.ysigLegend: legend for "x" vector
       # usig.ysigType  : ModiaResult.Continuous or ModiaResult.Clocked
       
# Simulate with Runge-Kutta 5/4 with step-size control
simulate!(firstOrder, Tsit5(), stopTime = 1.0)

# Simulate with Runge-Kutta 4 with fixed step size
simulate!(firstOrder, RK4(), stopTime = 1.0, adaptive=false)

# Simulate with algorithm that switches between 
# Verners Runge-Kutta 6/5 algorithm if non-stiff region and
# Rosenbrock 4 (= A-stable method) if stiff region with step-size control
simulate!(firstOrder, AutoVern6(Rodas4()), stopTime = 1.0)
```
"""
function simulate!(m::Nothing, args...; kwargs...)
    @info "The call of simulate!(..) is ignored, since the first argument is nothing."
    @test false
    return nothing
end
function simulate!(m::SimulationModel{FloatType,ParType,EvaluatedParType,TimeType}, algorithm=missing; merge=nothing, kwargs...) where {FloatType,TimeType,ParType,EvaluatedParType}
    enable_timer!(m.timer)
    reset_timer!(m.timer)
    TimerOutputs.@timeit m.timer "simulate!" begin 
        options = SimulationOptions{FloatType,TimeType}(merge; kwargs...)
        if isnothing(options)
            @test false
            return nothing
        end
        m.options = options
        if ismissing(algorithm) && FloatType == Float64
            algorithm = Sundials.CVODE_BDF()
        end
    
        # Initialize/re-initialize SimulationModel
        if m.options.log || m.options.logEvaluatedParameters || m.options.logStates
            println("... Simulate model ", m.modelName)
        end
        
        m.algorithmType = typeof(algorithm)       
        TimerOutputs.@timeit m.timer "init!" success = init!(m)
        if !success
            @test false
            return nothing
        end
    
        # Define problem and callbacks based on algorithm and model type
        interval = m.options.interval
        if  abs(m.options.stopTime - m.options.startTime) <= 0
            interval = 1.0
            tspan2   = [m.options.startTime]
        elseif abs(m.options.interval) < abs(m.options.stopTime-m.options.startTime)
            tspan2 = m.options.startTime:m.options.interval:m.options.stopTime 
        else
            tspan2 = [m.options.startTime, m.options.stopTime]
        end
        tspan   = (m.options.startTime, m.options.stopTime)         
        TimerOutputs.@timeit m.timer "ODEProblem"  problem = DifferentialEquations.ODEProblem{true}(derivatives!, m.x_init, tspan, m) 
        
        callback2 = DifferentialEquations.DiscreteCallback(timeEventCondition!, affectTimeEvent!)   
        eh = m.eventHandler            
        if eh.nz > 0        
            #println("\n!!! Callback set with crossing functions")
            # Due to DifferentialEquations bug https://github.com/SciML/DifferentialEquations.jl/issues/686
            # FunctionalCallingCallback(outputs!, ...) is not correctly called when zero crossings are present.
            # The fix is to call outputs!(..) from the previous to the current event, when an event occurs.
            # (alternativey: callback4 = DifferentialEquations.PresetTimeCallback(tspan2, affect_outputs!) )
            callback1 = DifferentialEquations.FunctionCallingCallback(outputs!, funcat=[m.options.startTime]) # call outputs!(..) at startTime
            callback3 = DifferentialEquations.VectorContinuousCallback(zeroCrossings!, 
                            affectStateEvent!, eh.nz, interp_points=m.options.interp_points)                                   
            #callback4 = DifferentialEquations.PresetTimeCallback(tspan2, affect_outputs!) 
            callbacks = DifferentialEquations.CallbackSet(callback1, callback2, callback3)   #, callback4)
        else
            #println("\n!!! Callback set without crossing functions")        
            callback1 = DifferentialEquations.FunctionCallingCallback(outputs!, funcat=tspan2)
            callbacks = DifferentialEquations.CallbackSet(callback1, callback2)
        end
        
        # Initial step size (the default of DifferentialEquations is too large) + step-size of fixed-step algorithm
        if !ismissing(algorithm) && typeof(algorithm) <: Sundials.CVODE_BDF
            cvode_bdf = true
        else
            cvode_bdf = false
            dt = m.options.adaptive ? m.options.interval/10 : m.options.interval   # initial step-size
        end
        m.addEventPointsDueToDEBug = cvode_bdf
        
        # Compute solution 
        abstol = 0.1*m.options.tolerance
        tstops = (m.eventHandler.nextEventTime,)
    
        if ismissing(algorithm)
            TimerOutputs.@timeit m.timer "solve" solution = DifferentialEquations.solve(problem, reltol=m.options.tolerance, abstol=abstol, save_everystep=false,
                                                callback=callbacks, adaptive=m.options.adaptive, saveat=tspan2, dt=dt, dtmax=m.options.dtmax, tstops = tstops)
        elseif cvode_bdf
            TimerOutputs.@timeit m.timer "solve" solution = DifferentialEquations.solve(problem, algorithm, reltol=m.options.tolerance, abstol=abstol, save_everystep=false,
                                                callback=callbacks, adaptive=m.options.adaptive, saveat=tspan2, dtmax=m.options.dtmax, tstops = tstops)
        else
            TimerOutputs.@timeit m.timer "solve" solution = DifferentialEquations.solve(problem, algorithm, reltol=m.options.tolerance, abstol=abstol, save_everystep=false,
                                                callback=callbacks, adaptive=m.options.adaptive, saveat=tspan2, dt=dt, dtmax=m.options.dtmax, tstops = tstops)
        end
        
        # Compute and store outputs from last event until final time
        sol_t = solution.t
        sol_x = solution.u
        m.storeResult = true        
        for i = length(m.result_vars)+1:length(sol_t)
            Base.invokelatest(m.getDerivatives!, m.der_x, sol_x[i], m, sol_t[i])
        end
        m.storeResult = false
        
        # Final update of instantiatedModel
        m.result_x = solution
        if ismissing(algorithm)
            m.algorithmType = typeof(solution.alg)
        end
        
        # Terminate simulation
        finalStates = solution.u[end]
        finalTime   = solution.t[end]
        terminate!(m, finalStates, finalTime)
    end
    disable_timer!(m.timer)  
    
    if !m.success
        return nothing
    end
    
    if m.options.log
        println("      Termination of ", m.modelName, " at time = ", finalTime, " s")
        println("        cpuTime         = ", round(TimerOutputs.time(m.timer["simulate!"])*1e-9, sigdigits=3), " s")
        println("        allocated       = ", round(TimerOutputs.allocated(m.timer["simulate!"])/1048576.0, sigdigits=3), " MiB")    
        println("        algorithm       = ", get_algorithmName(m))
        println("        FloatType       = ", FloatType)
        println("        interval        = ", m.options.interval, " s")
        println("        tolerance       = ", m.options.tolerance, " (relative tolerance)")
        println("        nEquations      = ", length(m.x_start))
        println("        nResults        = ", length(m.result_x.t))          
        println("        nGetDerivatives = ", m.nGetDerivatives, " (total number of getDerivatives! calls)")
        println("        nf              = ", m.nf, " (number of getDerivatives! calls from integrator)")  # solution.destats.nf
        println("        nZeroCrossings  = ", eh.nZeroCrossings, " (number of getDerivatives! calls for zero crossing detection)")

        if cvode_bdf && (eh.nTimeEvents > 0 || eh.nStateEvents > 0)
            # statistics is wrong, due to a bug in the Sundials.jl interface
            println("        nJac            = ??? (number of Jacobian computations)")
            println("        nAcceptedSteps  = ???")
            println("        nRejectedSteps  = ???")
            println("        nErrTestFails   = ???")
        else
            println("        nJac            = ", solution.destats.njacs, " (number of Jacobian computations)")
            println("        nAcceptedSteps  = ", solution.destats.naccept)
            println("        nRejectedSteps  = ", solution.destats.nreject)
            println("        nErrTestFails   = ", solution.destats.nreject)
        end
        println("        nTimeEvents     = ", eh.nTimeEvents)
        println("        nStateEvents    = ", eh.nStateEvents)
        println("        nRestartEvents  = ", eh.nRestartEvents)  
    end
    if m.options.logTiming
        println("\n... Timings for simulation of ", m.modelName,":")
        TimerOutputs.print_timer(m.timer; title="xxx")   #Timings for simulation of "* m.modelName * ":")
    end

    requiredFinalStates = m.options.requiredFinalStates
    if !ismissing(requiredFinalStates)
        if length(finalStates) != length(requiredFinalStates)
            success = false
        else
            success = finalStates == requiredFinalStates || isapprox(finalStates, requiredFinalStates, rtol=m.options.requiredFinalStates_rtol)
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
            @test finalStates == requiredFinalStates  || isapprox(finalStates, requiredFinalStates, rtol=m.options.requiredFinalStates_rtol)
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


#------------------------------------------------------------------------------------------------
#        Provide the overloaded ModiaResult Abstract Interface for the results of SimulationModel
#------------------------------------------------------------------------------------------------

ModiaResult.timeSignalName(  m::SimulationModel) = "time"
ModiaResult.hasOneTimeSignal(m::SimulationModel) = true


"""
    hasSignal(instantiatedModel, name::AbstractString)
    
Return true if parameter or time-varying variable `name` (for example `name = "a.b.c"`)
is defined in the instantiateModel that can be accessed and can be used for plotting.
"""
ModiaResult.hasSignal(m::SimulationModel, name::AbstractString) = 
    # m.save_x_in_solution ? name == "time" || haskey(m.equationInfo.x_dict, name) :
     haskey(m.result_info, name) || !ismissing(get_value(m.evaluatedParameters, name))

# For backwards compatibility
hasName(m::SimulationModel, name::AbstractString) = ModiaResult.hasSignal(m,name)



"""
    signalNames(instantiatedModel)
    
Return the variable names (parameters, time-varying variables) of an
[`@instantiateModel`](@ref) that can be accessed and can be used for plotting.
"""
function ModiaResult.signalNames(m::SimulationModel)
    #if m.save_x_in_solution 
    #    names = ["time"]
    #    append!(names, collect( keys(m.equationInfo.x_dict) ))
    #else
        all_names = get_names(m.evaluatedParameters)
        append!(all_names, collect( keys(m.result_info) ) )
    #end
    sort!(all_names)
    return all_names 
end

# For backwards compatibility
getNames(m::SimulationModel) = ModiaResult.signalNames(m)


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

function ModiaResult.rawSignal(m::SimulationModel, name::AbstractString)
    tsig = m.result_x.t
    if !m.unitless
        tsig = tsig*u"s"
        if !(m.options.desiredResultTimeUnit == NoUnits ||
            m.options.desiredResultTimeUnit == u"s")
            tsig = uconvert.(m.options.desiredResultTimeUnit, tsig)
        end
    end
        
    if name == "time"
        return ([tsig], [tsig], ModiaResult.Independent)
    end
    
    if haskey(m.result_info, name)
        resInfo = m.result_info[name]
        
        if resInfo.store == RESULT_X
            (ibeg,iend,xunit) = get_xinfo(m, resInfo.index) 
            if ibeg == iend
                xSig = ModiaResult.FlattenedSignalView(m.result_x.u, ibeg, (), resInfo.negate)
            else
                xSig = ModiaResult.FlattenedSignalView(m.result_x.u, ibeg, (iend-ibeg+1,), resInfo.negate)
            end
            if !m.unitless && xunit != ""
                xSig = xSig*uparse(xunit)
            end
            return ([tsig], [xSig], ModiaResult.Continuous)

        elseif resInfo.store == RESULT_DER_X
            (ibeg,iend,xunit) = get_xinfo(m, resInfo.index) 
            if ibeg == iend
                derxSig = ModiaResult.FlattenedSignalView(m.result_der_x, ibeg, (), resInfo.negate)
            else
                derxSig = ModiaResult.FlattenedSignalView(m.result_der_x, ibeg, (iend-ibeg+1,), resInfo.negate)
            end
            if !m.unitless
                if xunit == ""
                    derxSig = derxSig/u"s"
                else
                    derxSig = derxSig*(uparse(xunit)/u"s")
                end
            end            
            return ([tsig], [derxSig], ModiaResult.Continuous)   

        elseif resInfo.store == RESULT_VARS
            signal = ModiaResult.SignalView(m.result_vars, resInfo.index, resInfo.negate)
            if length(signal) != length(tsig)
                lens = length(signal)
                lent = length(tsig)
                error("Bug in SimulateAndPlot.jl (rawSignal(..)): name=\"$name\",\nlength(signal) = $lens, length(tsig) = $lent")
            end
            return ([tsig], [signal], ModiaResult.Continuous)
            
        elseif resInfo.store == RESULT_ZERO
            signal = ModiaResult.OneValueVector(0.0, length(tsig))
            return ([tsig], [signal], ModiaResult.Continuous)
            
        else
            error("Bug in SimulateAndPlot.jl (rawSignal(..)): name=\"$name\", resInfo=$resInfo")
        end  

    else
        value = get_value(m.evaluatedParameters, name)
        if ismissing(value)
            error("rawSignal: \"$name\" not in result of model $(m.modelName))")
        end
        signal = ModiaResult.OneValueVector(value, length(tsig))      
        return ([tsig], [signal], ModiaResult.Continuous)    
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


function get_algorithmName(m::SimulationModel)::String
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
    return algorithmName
end


function ModiaResult.defaultHeading(m::SimulationModel)
    FloatType = get_leaveName( string( typeof( m.x_start[1] ) ) )

    algorithmName = get_algorithmName(m)
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
using TinyModia
@usingModiaPlot
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


function setEvaluatedParametersInDataFrame!(obj::EvaluatedParType, result_info, dataFrame::DataFrames.DataFrame, path::String, nResult::Int)::Nothing where {EvaluatedParType}
    for (key,value) in zip(keys(obj), obj)
        name = appendName(path, key)
        if typeof(value) <: EvaluatedParType
            setEvaluatedParametersInDataFrame!(value, result_info, dataFrame, name, nResult)
        elseif !haskey(result_info, name)
            dataFrame[!,name] = ModiaResult.OneValueVector(value,nResult)
        end
    end
    return nothing
end

function get_result(m::SimulationModel; onlyStates=false, extraNames=missing) 
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
        for name in keys(m.result_info)
            if name != "time"
                (timeSignal, signal, signalType) = ModiaResult.rawSignal(m, name)
                dataFrame[!,name] = signal[1]
            end
        end
    
        setEvaluatedParametersInDataFrame!(m.evaluatedParameters, m.result_info, dataFrame, "", length(timeSignal[1]))
    end
    return dataFrame
end

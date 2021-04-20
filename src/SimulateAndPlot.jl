
export simulate!, get_result

using  ModiaPlot
using  Measurements
import MonteCarloMeasurements
using  Unitful
using  Test
import DifferentialEquations


#---------------------------------------------------------------------
#                          Simulation
#---------------------------------------------------------------------

"""
    convertTimeVariable(FloatType, t)

The function returns variable `t` in a normalized form:

1. If `t` has a unit, it is transformed to u"s" and the unit is stripped off.
2. `t` is converted to `FloatType`.
"""
convertTimeVariable(FloatType, t) = typeof(t) <: Unitful.AbstractQuantity ? convert(FloatType, ustrip(uconvert(u"s", t))) : convert(FloatType, t)


"""
    simulate!(model [, algorithm]; merge = nothing,
              tolerance = 1e-6, startTime = 0.0, stopTime = 1.0, interval = NaN,
              adaptive = true, log = true, logParameters = true, logStates = true,
              requiredFinalStates = nothing)

Simulate `model::SimulationModel` with `algorithm`
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
- `adaptive`: = true, if the `algorithm` should use step-size control (if available).
              = false, if the `algorithm` should use a fixed step-size of `interval` (if available).
- `log`: = true, to log the simulation.
- `logParameters`: = true, to log the parameter and init/start values
- `logStates` : = true, to log the states, its init/start values and its units.
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
simulate!(firstOrder, stopTime = 1.0, merge = Map(T = 0.6, x = 0.9), logParameters=true)

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
function simulate!(m::SimulationModel{FloatType,TimeType}, algorithm=missing;
                   tolerance = 1e-6,
                   startTime = 0.0,
                   stopTime  = 1.0,
                   interval  = NaN,
                   interp_points = 0,
                   nz = 0,
                   merge     = nothing,
                   adaptive::Bool      = true,
                   log::Bool           = false,
                   logParameters::Bool = false,
                   logStates::Bool     = false,
                   logEvents::Bool     = false,
                   logParameterExpressions::Bool = false,
                   requiredFinalStates           = nothing) where {FloatType,TimeType}
    #initialized = false
    #try
        eh = m.eventHandler
        if nz > 0
            eh.nz = nz
            eh.z  = ones(FloatType,nz)
            eh.zPositive = fill(false, nz)
        end
        if interp_points == 1
            # DifferentialEquations.jl crashes
            interp_points = 2
        end
        m.algorithmType = typeof(algorithm)
        startTime2 = convertTimeVariable(TimeType, startTime)
        stopTime2  = convertTimeVariable(TimeType, stopTime)
        interval2  = convertTimeVariable(TimeType, interval)
        interval2  = isnan(interval2) ? (stopTime2 - startTime2)/500.0 : interval2
        
        # Initialize/re-initialize SimulationModel
        if log || logParameters || logStates
            println("... Simulate model ", m.modelName)
        end
        cpuStart::UInt64 = time_ns()
        cpuLast::UInt64  = cpuStart
        cpuStartIntegration::UInt64 = cpuStart
        success = init!(m, startTime2, tolerance, merge, log, logParameterExpressions, logParameters, logStates, logEvents)
        if !success
            return nothing
        elseif m.eventHandler.restart == Terminate
            
        end
        #initialized = true
        
        if m.eventHandler.restart == Terminate
            finalStates = m.x_init
        else
            # Define problem and callbacks based on algorithm and model type
            if abs(interval2) < abs(stopTime2-startTime2)
                tspan2 = (startTime2+interval2):interval2:stopTime2
            else
                tspan2 = [stopTime2]
            end
           	tspan     = (startTime2, stopTime2)            
            abstol    = 0.1*tolerance
            problem   = DifferentialEquations.ODEProblem(derivatives!, m.x_init, tspan, m) 
            if logEvents
                println("      Number of zero crossing functions = ", eh.nz)
            end
            if eh.nz > 0
                # Due to DifferentialEquations bug https://github.com/SciML/DifferentialEquations.jl/issues/686
                # FunctionalCallingCallback(outputs!, ...) is not correctly called when zero crossings are present.
                # A temporary fix is to use time events at the communication points, but this slows down simulation.
                callback1 = DifferentialEquations.PresetTimeCallback(tspan2, affect_outputs!)  
                callback2 = DifferentialEquations.VectorContinuousCallback(conditions!, 
                                 affect!, eh.nz, interp_points=interp_points)                  
                callbacks = DifferentialEquations.CallbackSet(callback1, callback2)
            else
                callbacks = DifferentialEquations.FunctionCallingCallback(outputs!, funcat=tspan2, tdir=1)
            end
    
            # Initial step size (the default of DifferentialEquations is too large) + step-size of fixed-step algorithm
            dt = adaptive ? interval2/10 : interval2    # initial step-size
    
            # Compute solution
            solution = ismissing(algorithm) ? DifferentialEquations.solve(problem, reltol=tolerance, abstol=abstol,
                                                save_everystep=false, save_start=false, save_end=true,
                                                callback=callbacks, adaptive=adaptive, dt=dt) :
                                            DifferentialEquations.solve(problem, algorithm, reltol=tolerance, abstol=abstol,
                                                save_everystep=false, save_start=false, save_end=true,
                                                callback=callbacks, adaptive=adaptive, dt=dt)
            if ismissing(algorithm)
                m.algorithmType = typeof(solution.alg)
            end
    
            # Terminate simulation
            finalStates = solution[:,end]
            terminate!(m, finalStates, solution.t[end])
        end
        
        if log
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
            println("        interval        = ", interval2, " s")
            println("        tolerance       = ", tolerance, " (relative tolerance)")
            println("        nEquations      = ", length(m.x_start))
            println("        nResults        = ", length(m.result))
            println("        nAcceptedSteps  = ", solution.destats.naccept)
            println("        nRejectedSteps  = ", solution.destats.nreject)
            println("        nGetDerivatives = ", m.nGetDerivatives, " (total number of getDerivatives! calls)")
            println("        nf              = ", solution.destats.nf, " (number of getDerivatives! calls from integrator)")
            println("        nZeroCrossings  = ", eh.nZeroCrossings, " (number of getDerivatives! calls for zero crossing detection)")
            println("        nJac            = ", solution.destats.njacs, " (number of Jacobian computations)")
            println("        nErrTestFails   = ", solution.destats.nreject)          
           #println("        nTimeEvents     = ", eh.nTimeEvents)
            println("        nStateEvents    = ", eh.nStateEvents)
            println("        nRestartEvents  = ", eh.nRestartEvents)      
        end

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
#=
    catch e
        if initialized
            terminate!(m, m.x_start, m.time)
        end
        
        if isa(e, ErrorException)
            printstyled("\nError from simulate!:\n", color=:red)
            printstyled(e.msg, "\n\n", color=:red)
            printstyled("Aborting simulate!\n\n", color=:red)
            empty!(m.result)
        else
            Base.rethrow()
        end
    end

    return nothing
=#
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
    haskey(m.variables, name) || 
    name in m.zeroVariables ||
    !ismissing(get_value(m.parameters, name))


"""
    ModiaPlot.getNames(model::SimulationModel)
    
Return the variable names (parameters, time-varying variables) of a TinyModia SimulationModel
(generated with [`TinyModia.@instantiateModel`](@ref) that can be accessed
and can be used for plotting.
"""
function ModiaPlot.getNames(m::SimulationModel)
    names = get_names(m.parameters)
    append!(names, collect(m.zeroVariables))
    append!(names, collect( keys(m.variables) ) )
    return sort(names)
end

function ModiaPlot.getRawSignal(m::SimulationModel, name)
    if haskey(m.variables, name)
        resIndex = m.variables[name]
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
            signal = -signal
        end
        return (false, signal)
        
    elseif name in m.zeroVariables
        return (true, 0.0)

    else
        value = get_value(m.parameters, name)
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
    signal = get_result(model, name; unit=true)

After a successful simulation of `model::SimulationModel`, return
the result for the signal `name::String` as vector of points
together with its unit. The time vector has path name `"time"`.
If `unit=false`, the signal is returned, **without unit**.

# Example

```julia
using ModiaBase
using Unitful

include("\$(ModiaBase.path)/demos/models/Model_Pendulum.jl")
using  .Model_Pendulum

pendulum = simulationModel(Pendulum)
simulate!(pendulum, stopTime=7.0)

time = get_result(pendulum, "time")  # vector with unit u"s"
phi  = get_result(pendulum, "phi")   # vector with unit u"rad"

import PyPlot
PyPlot.figure(4)   # Change to figure 4 (or create it, if it does not exist)
PyPlot.clf()       # Clear current figure
PyPlot.plot(ustrip(time), ustrip(phi), "b--", label="phi in " * string(unit(phi[1])))
PyPlot.xlabel("time in " * string(unit(time[1])))
PyPlot.legend()
```
"""
function get_result(m::SimulationModel, name::String; unit=true)
    #(xsig, xsigLegend, ysig, ysigLegend, yIsConstant) = ModiaPlot.getPlotSignal(m, "time", name)

    (isConstant, ysig) = ModiaPlot.getRawSignal(m, name)

    ysig = unit ? ysig : ustrip.(ysig)


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

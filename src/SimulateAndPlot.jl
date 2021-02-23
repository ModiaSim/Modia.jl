
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
    convertTimeVariable(t)
    
The function returns variable `t` in a normalized form:

1. If `t` has a unit, it is transformed to u"s" and the unit is stripped off.
2. `t` is converted to `Float64`.
"""
convertTimeVariable(t) = typeof(t) <: Unitful.AbstractQuantity ? convert(Float64, ustrip(uconvert(u"s", t))) : convert(Float64, t)



"""
    simulate!(model [, algorithm];
              tolerance=1e-6, startTime=0.0, stopTime=1.0, interval=NaN, 
              adaptive=true, log=true, requiredFinalStates=nothing)

Simulate `model::SimulationModel` with `algorithm` 
(= `alg` of [ODE Solvers of DifferentialEquations.jl](https://diffeq.sciml.ai/stable/solvers/ode_solve/)).
If the `algorithm` argument is missing, a default algorithm will be chosen from DifferentialEquations
(for details see [https://arxiv.org/pdf/1807.06430](https://arxiv.org/pdf/1807.06430), Figure 3).

The simulation results stored in `model` can be plotted with ModiaPlot.plot and the result values
can be retrieved with [`get_result`](@ref).


# Optional Arguments

- `tolerance`: Relative tolerance.
- `startTime`: Start time. If value is without unit, it is assumed to have unit [s].
- `stopTime`: Stop time. If value is without unit, it is assumed to have unit [s].
- `interval`: Interval to store result. If `interval=NaN`, it is internally selected as
              (stopTime-startTime)/500.
              If value is without unit, it is assumed to have unit [s].
- `adaptive`: = true, if the `algorithm` should use step-size control (if available).
              = false, if the `algorithm` should use a fixed step-size of `interval` (if available).
- `log`: = true, to log the simulation.
- `requiredFinalStates`: is not `nothing`: Check whether the ODE state vector at the final time instant with `@test`
              is in agreement to vector `requiredFinalStates` with respect to some tolerance. If this is not the case, print the
              final state vector (so that it can be included with copy-and-paste in the simulate!(..) call).

# Examples

```julia
import ModiaPlot
using  DifferentialEquations

# Runge-Kutta 5/4 with step-size control
simulate!(model, DifferentialEquations.Tsit5(), stopTime = 1.0)

    # Plot variables "v1", "v2" in diagram 1, "v3" in diagram 2, both diagrams in figure 3   
    ModiaPlot.plot(model, [("v1","v2"), "v3"], figure=3) 

    # Retrieve "time" and "v1" values:
    get_result(model, "time")
    get_result(model, "v1")

# Runge-Kutta 4 with fixed step size
simulate!(model, DifferentialEquations.RK4(), stopTime = 1.0, adaptive=false)

# Switching between Verners Runge-Kutta 6/5 algorithm if non-stiff region and 
# Rosenbrock 4 (= A-stable method) if stiff region with step-size control
simulate!(model, AutoVern6(Rodas4()), stopTime = 1.0)

# Automatically selected algorithm
simulate!(model, stopTime = 1.0)

# Sundials CVODE (BDF method with variable order 1-5) with step-size control
using Sundials
simulate!(model, CVODE_BDF(), stopTime = 1.0)
```
"""
function simulate!(m::Nothing, args...; kwargs...)
    @info "The call of simulate!(..) is ignored, since the first argument is nothing."
    return
end  
function simulate!(m::SimulationModel, algorithm=missing;
                   tolerance=1e-6,
                   startTime=0.0,
                   stopTime=1.0,
                   interval=NaN,
				   adaptive::Bool=true,
                   log::Bool=false,
                   requiredFinalStates=nothing)
    try                
        m.algorithmType = typeof(algorithm)
        tolerance = convert(Float64, tolerance)    
        startTime = convertTimeVariable(startTime)
        stopTime  = convertTimeVariable(stopTime)
        interval  = convertTimeVariable(interval)
        interval  = isnan(interval) ? (stopTime - startTime)/500.0 : interval
    
        # Determine x_start (currently: x_start = 0)
        # println("x_start = ", m.x_start)
    
        # Target type
        FloatType = getFloatType(m)
        
        # Initialize/re-initialize SimulationModel
        if log
            println("... Simulate model ", m.modelName)        
            cpuStart::UInt64 = time_ns()
            cpuLast::UInt64  = cpuStart
            cpuStartIntegration::UInt64 = cpuStart      
            println("      Initialization at time = ", startTime, " s")
        end    
        init!(m, startTime, tolerance)
    
    
        # Define problem and callbacks based on algorithm and model type
        tspan    = (startTime, stopTime)
        tspan2   = startTime:interval:stopTime	
        tdir     = startTime <= stopTime ? 1 : -1   
        abstol   = 0.1*tolerance
        problem  = DifferentialEquations.ODEProblem(derivatives!, m.x_start, tspan, m)
        callback = DifferentialEquations.FunctionCallingCallback(outputs!, funcat=tspan2, tdir=tdir)
    
    
        # Initial step size (the default of DifferentialEquations is too large) + step-size of fixed-step algorithm
        dt = adaptive ? interval/10 : interval    # initial step-size
            
        # Compute solution
        solution = ismissing(algorithm) ? DifferentialEquations.solve(problem, reltol=tolerance, abstol=abstol,
                                            save_everystep=false, save_start=false, save_end=true, 
                                            callback=callback, adaptive=adaptive, dt=dt) :
                                        DifferentialEquations.solve(problem, algorithm, reltol=tolerance, abstol=abstol,
                                            save_everystep=false, save_start=false, save_end=true, 
                                            callback=callback, adaptive=adaptive, dt=dt)
        if ismissing(algorithm)
            m.algorithmType = typeof(solution.alg)
        end
        
        # Terminate simulation      
        if log
            cpuTimeInitialization = convert(Float64, (cpuStartIntegration - cpuStart) * 1e-9)
            cpuTimeIntegration    = convert(Float64, (time_ns() - cpuStartIntegration) * 1e-9)
            cpuTime               = cpuTimeInitialization + cpuTimeIntegration
            
            println("      Termination at time    = ", solution.t[end], " s")
            println("        cpuTime         = ", round(cpuTime, sigdigits=3), " s")        
            #println("        cpuTime         = ", round(cpuTime              , sigdigits=3), " s (init: ", 
            #                                      round(cpuTimeInitialization, sigdigits=3), " s, integration: ", 
            #                                      round(cpuTimeIntegration   , sigdigits=3), " s)")
            println("        algorithm       = ", typeof(solution.alg))
            println("        FloatType       = ", FloatType)        
            println("        interval        = ", interval, " s")
            println("        tolerance       = ", tolerance, " (relative tolerance)")
            println("        nEquations      = ", length(m.x_start))
            println("        nResults        = ", length(m.result))
            println("        nAcceptedSteps  = ", solution.destats.naccept)
            println("        nRejectedSteps  = ", solution.destats.nreject)        
            println("        nGetDerivatives = ", m.nGetDerivatives, " (number of getDerivatives! calls)")
            println("        nf              = ", solution.destats.nf, " (number of getDerivatives! calls from integrator)")        
            println("        nJac            = ", solution.destats.njacs, " (number of Jacobian computations)")
            println("        nErrTestFails   = ", solution.destats.nreject)       
        end                                           
    
        finalStates = solution[:,end] 
        
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

    catch e
        if isa(e, ErrorException)
            printstyled("\nError from simulate!:\n", color=:red)
            printstyled(e.msg, "\n\n", color=:red)
            printstyled("Aborting simulate!\n\n", color=:red)
        else
            Base.rethrow()
        end
    end
    return nothing
end


#---------------------------------------------------------------------
#        Provide ModiaPlot result access functions for SimulationModel
#---------------------------------------------------------------------
    
ModiaPlot.hasSignal(m::SimulationModel, name) = 
    haskey(m.variables, name) || haskey(m.parametersAndConstantVariables, name)

ModiaPlot.getNames(m::SimulationModel) = 
    append!(collect( keys(m.parametersAndConstantVariables) ),
            collect( keys(m.variables) ) ) 
                                     
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
    else
        return (true, eval( m.parametersAndConstantVariables[name] ) )
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

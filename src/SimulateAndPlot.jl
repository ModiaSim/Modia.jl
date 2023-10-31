using  Test
import DataFrames
import ForwardDiff
import FiniteDiff

#---------------------------------------------------------------------
#                          Simulation
#---------------------------------------------------------------------

function getAlgorithmName(algorithm)::String
    algorithmType = typeof(algorithm)
    if algorithmType == Missing
        return "???"
    end
    name = string(algorithmType)

    if algorithmType <: DifferentialEquations.OrdinaryDiffEq.QNDF
        if algorithm.kappa == tuple(0//1,0//1,0//1,0//1,0//1)
            name = replace(name, "QNDF" => "QBDF")
        end

    elseif algorithmType <: DifferentialEquations.OrdinaryDiffEq.QNDF1 ||
           algorithmType <: DifferentialEquations.OrdinaryDiffEq.QNDF2
        if algorithm.kappa == 0
            name = replace(name, "QNDF" => "QBDF")
        end
    end
    return name
end

"""
    solution = simulate!(instantiatedModel [, algorithm];
              merge            = missing,  # change parameter/init/start values
              tolerance        = 1e-6,     # relative tolerance
              startTime        = 0.0,
              stopTime         = 0.0,      # stopTime >= startTime required
              interval         = missing,  # = (stopTime-startTime)/500
              interp_points    = 0,
              dtmax            = missing,  # = 100*interval
              adaptive         = true,
              nlinearMinForDAE = 10,
              log              = false,
              logStates        = false,
              logEvents        = false,
              logProgress      = false,
              logTiming        = false,
              logParameters    = false,
              logEvaluatedParameters   = false,
              requiredFinalStates      = missing
              requiredFinalStates_rtol = 1e-3,
              requiredFinalStates_atol = 0.0,
              useRecursiveFactorizationUptoSize = 0)

Simulate `instantiatedModel::InstantiatedModel` with `algorithm`
(= `alg` of [ODE Solvers of DifferentialEquations.jl](https://diffeq.sciml.ai/stable/solvers/ode_solve/)
or [DAE Solvers of DifferentialEquations.jl](https://diffeq.sciml.ai/stable/solvers/dae_solve/)).

If the `algorithm` argument is missing, `algorithm=Sundials.CVODE_BDF()` is used, provided
instantiatedModel has `FloatType = Float64`. Otherwise, a default algorithm will be chosen from DifferentialEquations
(for details see [https://arxiv.org/pdf/1807.06430](https://arxiv.org/pdf/1807.06430), Figure 3).
The symbols `CVODE_BDF` and `IDA` are exported from Modia, so that `simulate!(instantiatedModel, CVODE_BDF(), ...)`
and `simulate!(instantiatedModel, IDA(), ...)`
can be used (instead of `import Sundials; simulate!(instantiatedModel, Sundials.xxx(), ...)`).

The simulation results are stored in `instantiatedModel` and can be plotted with
`plot(instantiatedModel, ...)`. The result values
can be retrieved with `getValues(..)` for Var(..) and `getValue(..)` for Par(..).
`showInfo(instantiatedModel)` prints information about the signals in the result.
For more details, see sections [Parameters/Init/Start](@ref), [Results and Plotting](@ref).

The return argument `solution` is the return argument from `DifferentialEquations.solve(..)` and
therefore all post-processing functionality from `DifferentialEqautions.jl` can be used. Especially,
- solution.t[i] # time-instant at storage point i (solution.t[end] = stopTime)
- solution.u[i] # states at storage point i

A simulation run can be aborted with `<CTRL> C` (SIGINT), provided `using PyPlot` or `import PyPlot` was
not called before (the signals in Python module matplotlib.pyplot intervene with Julias signals, see
[PyPlot.jl issue 305](https://github.com/JuliaPy/PyPlot.jl/issues/305)).

# Optional ArgumentsS

- `merge`: Define parameters and init/start values that shall be merged with the previous values
           stored in `model`, before simulation is started. If, say, an init value `phi = Var(init=1.0)`
           is defined in the model, a different init value can be provided with
           `merge = Map(phi=2.0)`.
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
- `nlinearMinForDAE`: If `algorithm` is a DAE integrator (e.g. `IDA()`) and the size of a linear equation system
              is `>= nlinearMinForDAE` and the iteration variables of this equation system are a subset of the
              DAE state derivatives, then during continuous integration (but not at events, including
              initialization) this equation system is not locally solved but is solved via the DAE integrator.
              Typically, for large linear equation systems, simulation efficiency is considerably improved
              in such a case.f
- `log`: = true, to log the simulation.
- `logStates`: = true, to log the states, its init/start values and its units.
- `logEvents`: = true, to log events.
- `logProgress` = true, to printout current simulation time every 5s.
- `logTiming`: = true, to log the timing with `instantiatedModel.timer` which is an instance
               of [TimerOutputs](https://github.com/KristofferC/TimerOutputs.jl).TimerOutput.
               A user function can include its timing via\\
               `TimerOutputs.@timeit instantiatedModel.timer "My Timing" <statement>`.
- `logParameters`: = true, to log parameters and init/start values defined in model.
- `logEvaluatedParameters`: = true, to log the evaluated parameter and init/start values that
                            are used for initialization and during simulation.
- `requiredFinalStates`: is not `missing`: Test with `@test` whether the ODE state vector at the
              final time instant is in agreement to vector `requiredFinalStates` with respect
              to tolerances `requiredFinalStates_rtol, requiredFinalStates_atol`. If this is not the case, print the
              final state vector (so that it can be included with copy-and-paste in the simulate!(..) call).
              If you checked that the result of the simulation is correct, use `requiredFinalStates = []` to get
              a printout of the required final states and then copy it in your test.
- `requiredFinalStates_rtol`: Relative tolerance used for `requiredFinalStates`.
- `requiredFinalStates_atol`: Absolute tolerance used for `requiredFinalStates` (see atol in `?isapprox`)
- `useRecursiveFactorizationUptoSize`: = 0: Linear equation systems A*v=b are solved with
               `RecursiveFactorization.jl` instead of the default `lu!(..)` and `ldiv!(..)`, if
               `length(v) <= useRecursiveFactorizationUptoSize`.
               According to `RecursiveFactorization.jl` docu, it is faster as `lu!(..)` with OpenBLAS,
               for `length(v) <= 500` (typically, more as a factor of two).
               Since there had been some cases where `lu!(..)!` was successful,
               but `RecursiveFactorization.jl` failed due to a singular system, the default is to use `lu!(..)!`.

# Examples

```julia
using Modia
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

function simulate!(m::InstantiatedModel{FloatType,TimeType}, algorithm=missing; merge=nothing, kwargs...) where {FloatType,TimeType}
    options = SimulationOptions{FloatType,TimeType}(merge; kwargs...)
    if isnothing(options)
        @test false
        return nothing
    end
    m.options   = options
    m.time      = options.startTime
    m.isInitial = true
    m.nsegments = 1
    reinitEventHandler!(m.eventHandler, m.options.stopTime, m.options.logEvents)

    if ismissing(algorithm) && FloatType == Float64
        algorithm = Sundials.CVODE_BDF()
    end
    m.algorithmName = getAlgorithmName(algorithm)
    m.sundials      = !ismissing(algorithm) && (typeof(algorithm) <: Sundials.SundialsODEAlgorithm || typeof(algorithm) <: Sundials.SundialsDAEAlgorithm)
    m.addEventPointsDueToDEBug = m.sundials
    m.odeIntegrator = !(typeof(algorithm) <: DifferentialEquations.DiffEqBase.AbstractDAEAlgorithm)

    # Initialize/re-initialize InstantiatedModel
    if m.options.log || m.options.logEvaluatedParameters || m.options.logStates
        println("\n... Simulate model ", m.modelName)
    end

    enable_timer!(m.timer)
    reset_timer!(m.timer)
    TimerOutputs.@timeit m.timer "Modia.init!"  begin
        if m.options.log || m.options.logTiming
            timedStats = @timed @time (success = init!(m); if m.options.log || m.options.logTiming; print("      Initialization finished within") end)
        else
            timedStats = @timed (success = init!(m))
        end
    end
    updateStatistics_initStat!(m.statistics, timedStats.time, timedStats.bytes*1e-6)
    if !success
        @test false
        return nothing
    end

    m.cpuLast  = time_ns()
    m.cpuFirst = m.cpuLast
    eh = m.eventHandler
    solution = nothing

    # --------------------------------------------------------- Core simulation loop ---------------------------------
    TimerOutputs.@timeit m.timer "Modia.simulate!" while true
        solution = simulateSegment!(m, algorithm; kwargs...)

        # Terminate simulation of current segment
        finalStates = solution.u[end]
        finalTime   = solution.t[end]
        #m.result_x = solution
        terminate!(m, finalStates, finalTime)

        # Update statistics
        if m.sundials && (eh.nTimeEvents > 0 || eh.nStateEvents > 0)
            nJac           = missing
            nAcceptedSteps = missing
            nRejectedSteps = missing
        else
            nJac           = solution.destats.njacs
            nAcceptedSteps = solution.destats.naccept
            nRejectedSteps = solution.destats.nreject
        end
        updateStatistics_performanceIndictors!(m.statistics, finalTime,
            nResultsCurrentSegment(m.result),
            m.nf_total, m.nf_integrator, eh.nZeroCrossings,
            nJac, nAcceptedSteps, nRejectedSteps,
            eh.nTimeEvents, eh.nStateEvents, eh.nRestartEvents,
            eh.nFullRestartEvents, m.linearEquations)

        # Raise an error, if simulation was not successful
        if !DifferentialEquations.SciMLBase.successful_retcode(solution.retcode)
            error("\nsolution = simulate!(", m.modelName, ", ...) failed with solution.retcode = :$(solution.retcode) at time = $finalTime.\n")
        end


        if m.eventHandler.restart != Modia.FullRestart
            break
        end

        # Re-initialize model for FullRestart
        initFullRestart!(m)
    end # --------------------------------------------------------- End of core simulation loop ------------------------------

    disable_timer!(m.timer)
    updateStatistics_final!(m.statistics,
                               TimerOutputs.time(m.timer["Modia.simulate!"])*1e-9,
                               TimerOutputs.allocated(m.timer["Modia.simulate!"])*1e-6)
    if !m.success
        return nothing
    end

    if m.options.log
        eh = m.eventHandler
        println("      Termination of ", m.modelName, " at time = ", solution.t[end], " s")
        #=
        useRecursiveFactorization = Bool[leq.useRecursiveFactorization for leq in m.linearEquations]
        println("        cpuTime (without init.)   = ", round(TimerOutputs.time(m.timer["Modia.simulate!"])*1e-9, sigdigits=3), " s")
        println("        allocated (without init.) = ", round(TimerOutputs.allocated(m.timer["Modia.simulate!"])/1048576.0, sigdigits=3), " MiB")
        println("        algorithm                 = ", get_algorithmName_for_heading(m))
        println("        FloatType                 = ", FloatType)
        println("        interval                  = ", m.options.interval, " s")
        println("        tolerance                 = ", m.options.tolerance, " (relative tolerance)")
        println("        nStates                   = ", length(m.x_start))
        println("        linearSystemSizes         = ", Int[length(leq.b) for leq in m.linearEquations])
        println("        useRecursiveFactorization = ", useRecursiveFactorization)
        println("        odeModeLinearSystems      = ", Bool[leq.odeMode for leq in m.linearEquations])
        println("        nResults                  = ", nResults(m.result))
        println("        nf_total                  = ", m.nf_total, " (total number of getDerivatives! calls)")
        println("        nf_integrator             = ", m.nf_integrator, " (number of getDerivatives! calls from integrator)")  # solution.destats.nf
        println("        nf_zeroCrossings          = ", eh.nZeroCrossings, " (number of getDerivatives! calls for zero crossing detection)")

        if m.sundials && (eh.nTimeEvents > 0 || eh.nStateEvents > 0)
            # statistics is wrong, due to a bug in the Sundials.jl interface
            println("        nJac                      = ??? (number of Jacobian computations)")
            println("        nAcceptedSteps            = ???")
            println("        nRejectedSteps            = ???")
        else
            println("        nJac                      = ", solution.destats.njacs, " (number of Jacobian computations)")
            println("        nAcceptedSteps            = ", solution.destats.naccept)
            println("        nRejectedSteps            = ", solution.destats.nreject)
        end
        println("        nTimeEvents               = ", eh.nTimeEvents)
        println("        nStateEvents              = ", eh.nStateEvents)
        if eh.nFullRestartEvents > 0
            println("        nFullRestartEvents        = ", eh.nFullRestartEvents)
        end
        println("        nRestartEvents            = ", eh.nRestartEvents)
        =#
        printStatistics(m.statistics, m.options.interval, m.options.tolerance)
    end
    if m.options.logTiming
        println("\n... Timings for simulation of ", m.modelName," (without initialization):")
        TimerOutputs.print_timer(TimerOutputs.flatten(m.timer), compact=true)
    end

    requiredFinalStates = m.options.requiredFinalStates
    if !ismissing(requiredFinalStates)
        finalStates = solution.u[end]
        rtol = m.options.requiredFinalStates_rtol
        atol = m.options.requiredFinalStates_atol
        if length(finalStates) != length(requiredFinalStates)
            success = false
        else
            success = isapprox(finalStates, requiredFinalStates, rtol=rtol, atol=atol)
        end

        if success
            @test success
        else
            println("\nrequiredFinalStates_rtol = $rtol")
            println("requiredFinalStates_atol = $atol")
            if length(requiredFinalStates) > 0 && typeof(requiredFinalStates[1]) <: Measurements.Measurement
                println(  "\nrequiredFinalStates   = ", measurementToString(requiredFinalStates))
                printstyled("finalStates           = ", measurementToString(finalStates), "\n\n", bold=true, color=:red)
                if length(finalStates) == length(requiredFinalStates)
                    printstyled("difference            = ", measurementToString(requiredFinalStates-finalStates), "\n\n", bold=true, color=:red)
                end
            else
                println(  "\nrequiredFinalStates   = ", requiredFinalStates)
                printstyled("finalStates           = ", finalStates, "\n\n", bold=true, color=:red)
                if length(finalStates) == length(requiredFinalStates)
                    printstyled("difference            = ", requiredFinalStates-finalStates, "\n\n", bold=true, color=:red)
                end
            end
            if length(finalStates) == length(requiredFinalStates)
                printstyled("maximum(|difference|) = ", maximum(abs.(requiredFinalStates-finalStates)), "\n\n", bold=true, color=:red)
            end
            @test isapprox(finalStates, requiredFinalStates, rtol=rtol, atol=atol)
        end
    end

    return solution

    #=
    catch e
        if isa(e, ErrorException)
            println()
            printstyled("Error during simulation at time = $(m.time) s:\n\n", bold=true, color=:red)
            printstyled(e.msg, "\n", bold=true, color=:red)
            printstyled("\nAborting simulate!(..) for model $(m.modelName) instantiated in file\n$(m.modelFile).\n", bold=true, color=:red)
            println()
            m.lastMessage = deepcopy(e.msg)
            #@test false
        elseif isa(e, InterruptException)
            println()
            m.lastMessage = "<ctrl> C interrupt during simulation at time = $(m.time) s.\n"
            printstyled(m.lastMessage, bold=true, color=:red)
            printstyled("\nAborting simulate!(..) for model $(m.modelName) instantiated in file\n$(m.modelFile).", bold=true, color=:red)
            println()
        else
            println("... in else branch")
            Base.rethrow()
        end
    end
    =#
end


function simulateSegment!(m::InstantiatedModel{FloatType,TimeType}, algorithm=missing; kwargs...) where {FloatType,TimeType}
    solution = nothing
    options  = m.options

    sizesOfLinearEquationSystems = Int[length(leq.b) for leq in m.linearEquations]

    # Define problem and callbacks based on algorithm and model type
    interval = options.interval
    if  abs(options.stopTime - options.startTime) <= 0
        interval = 1.0
        tspan2   = [options.startTime]
    elseif abs(options.interval) < abs(options.stopTime-options.startTime)
        if m.nsegments == 1
            tspan2 = options.startTime:interval:options.stopTime
        else
            i      = ceil( (options.startTime - options.startTimeFirstSegment)/interval )
            tnext  = options.startTimeFirstSegment + i*interval
            tspan2 = tnext:interval:options.stopTime
            if tspan2[1] > options.startTime
                tspan2 = [options.startTime, tspan2...]
            end
        end
        if tspan2[end] < options.stopTime
            tspan2 = [tspan2..., options.stopTime]
        end
    else
        tspan2 = [options.startTime, options.stopTime]
    end
    tspan = (options.startTime, options.stopTime)
    eh = m.eventHandler
    for leq in m.linearEquations
        leq.odeMode = true
    end
    if m.odeIntegrator
        # ODE integrator
        TimerOutputs.@timeit m.timer "DifferentialEquations.ODEProblem" problem = DifferentialEquations.ODEProblem{true,DifferentialEquations.SciMLBase.FullSpecialize}(derivatives!, m.x_init, tspan, m)
    else
        # DAE integrator
        nx = length(m.x_init)
        differential_vars = eh.nz > 0 ? fill(true, nx) : nothing    # due to DifferentialEquations issue #549
        copyDerivatives!(m.der_x, m.der_x_invariant, m.der_x_segmented)
        TimerOutputs.@timeit m.timer "DifferentialEquations.DAEProblem" problem = DifferentialEquations.DAEProblem{true}(DAEresidualsForODE!, m.der_x, m.x_init, tspan, m, differential_vars = differential_vars)
        empty!(m.daeCopyInfo)
        if length(sizesOfLinearEquationSystems) > 0 && maximum(sizesOfLinearEquationSystems) >= options.nlinearMinForDAE
            # Prepare data structure to efficiently perform copy operations for DAE integrator
            x_info      = m.equationInfo.x_info
            der_x_dict  = m.equationInfo.der_x_dict
            der_x_names = keys(der_x_dict)
            daeMode     = false
            for (ileq,leq) in enumerate(m.linearEquations)
                if sizesOfLinearEquationSystems[ileq] >= options.nlinearMinForDAE
                    answer  = leq.x_names .∈ Ref(der_x_names)
                    daeMode = true
                    for (index, val) in enumerate(answer)
                        if !val && leq.x_lengths[index] > 0
                            daeMode = false
                            break
                        end
                    end

                    if daeMode
                        # Linear equation shall be solved by DAE and all unknowns of the linear equation system are DAE derivatives
                        if eh.nz > 0
                            leq.odeMode = true
                            daeMode = false
                            if options.log
                                println("      No DAE mode for equation system $ileq because $(eh.nz) crossing function(s) defined (see issue #686 of DifferentialEquations.jl)")
                            end

                        else
                            leq.odeMode = false
                            leq_copy = LinearEquationsCopyInfoForDAEMode(ileq)
                            for ix in 1:length(leq.x_names)
                                x_length = leq.x_lengths[ix]
                                if x_length > 0
                                    x_name   = leq.x_names[ix]
                                    x_info_i = x_info[ der_x_dict[x_name] ]
                                    @assert(x_length == x_info_i.length)
                                    startIndex = x_info_i.startIndex
                                    endIndex   = startIndex + x_length - 1
                                    append!(leq_copy.index, startIndex:endIndex)
                                end
                            end
                            push!(m.daeCopyInfo, leq_copy)
                        end
                    else
                        if options.log
                            unknownsThatAreNoStateDerivatives = ""
                            first = true
                            for (index, val) in enumerate(answer)
                                if !val && leq.x_lengths[index] > 0
                                    if first
                                        first = false
                                        unknownsThatAreNoStateDerivatives = "\"" * leq.x_names[index] * "\""
                                    else
                                        unknownsThatAreNoStateDerivatives *= ",\"" * leq.x_names[index] * "\""
                                    end
                                end
                            end
                            println("      No DAE mode for equation system $ileq because the unknowns $unknownsThatAreNoStateDerivatives are no state derivatives!")
                        end
                        leq.odeMode = true
                    end
                else
                    leq.odeMode = true
                end
            end
        elseif length(sizesOfLinearEquationSystems) > 0
            (leqSizeMax, ileq) = findmax(sizesOfLinearEquationSystems)
            println("      No DAE mode for equation system $ileq because size of equation system (= $leqSizeMax) < nlinearMinForDAE (= $(options.nlinearMinForDAE)).")
        end
    end

    if length(m.linearEquations) == 0
        m.odeMode   = true
        m.solve_leq = true
    else
        m.odeMode   = false
        m.solve_leq = false
        for leq in m.linearEquations
            if leq.odeMode
                m.odeMode   = true
                m.solve_leq = true
                break
            end
        end
    end

    callback2 = DifferentialEquations.DiscreteCallback(timeEventCondition!, affectTimeEvent!)
    if eh.nz > 0
        #println("\n!!! Callback set with crossing functions")
        # Due to DifferentialEquations bug https://github.com/SciML/DifferentialEquations.jl/issues/686
        # FunctionalCallingCallback(outputs!, ...) is not correctly called when zero crossings are present.
        # The fix is to call outputs!(..) from the previous to the current event, when an event occurs.
        # (alternativey: callback4 = DifferentialEquations.PresetTimeCallback(tspan2, affect_outputs!) )
        callback1 = DifferentialEquations.FunctionCallingCallback(outputs!, funcat=[options.startTime]) # call outputs!(..) at startTime
        callback3 = DifferentialEquations.VectorContinuousCallback(zeroCrossings!,
                        affectStateEvent!, eh.nz, interp_points=options.interp_points, rootfind=DifferentialEquations.SciMLBase.RightRootFind)
        #callback4 = DifferentialEquations.PresetTimeCallback(tspan2, affect_outputs!)
        callbacks = DifferentialEquations.CallbackSet(callback1, callback2, callback3)   #, callback4)
    else
        #println("\n!!! Callback set without crossing functions")
        callback1 = DifferentialEquations.FunctionCallingCallback(outputs!, funcat=tspan2)
        callbacks = DifferentialEquations.CallbackSet(callback1, callback2)
    end

    # Initial step size (the default of DifferentialEquations integrators is too large) + step-size of fixed-step algorithm
    if !m.sundials
        dt = options.adaptive ? options.interval/10 : options.interval   # initial step-size
    end

    # Compute solution
    abstol = 0.1*options.tolerance
    tstops = [m.eventHandler.nextEventTime,]
    maxiters = Int(typemax(Int32))  # switch off maximum number of iterations (typemax(Int) gives an inexact error for Sundials)
    if ismissing(algorithm)
        TimerOutputs.@timeit m.timer "DifferentialEquations.solve" solution = DifferentialEquations.solve(problem, reltol=options.tolerance, abstol=abstol, save_everystep=false,
                                                                        callback=callbacks, adaptive=options.adaptive, saveat=tspan2, dt=dt, dtmax=options.dtmax, maxiters=maxiters, tstops = tstops,
                                                                        initializealg = DifferentialEquations.NoInit())
    elseif m.sundials
        TimerOutputs.@timeit m.timer "DifferentialEquations.solve" solution = DifferentialEquations.solve(problem, algorithm, reltol=options.tolerance, abstol=abstol, save_everystep=false,
                                                                        callback=callbacks, adaptive=options.adaptive, saveat=tspan2, dtmax=options.dtmax, maxiters=maxiters, tstops = tstops,
                                                                        initializealg = DifferentialEquations.NoInit())
    else
        TimerOutputs.@timeit m.timer "DifferentialEquations.solve" solution = DifferentialEquations.solve(problem, algorithm, reltol=options.tolerance, abstol=abstol, save_everystep=false,
                                                                        callback=callbacks, adaptive=options.adaptive, saveat=tspan2, dt=dt, dtmax=options.dtmax, maxiters=maxiters, tstops = tstops,
                                                                        initializealg = DifferentialEquations.NoInit())
    end

    # Compute and store outputs from last event until final time
    sol_t = solution.t
    sol_x = solution.u
    m.storeResult = true
#   for i = length(m.result_code)+1:length(sol_t)
    for i = length(m.result.t[end])+1:length(sol_t)
        invokelatest_getDerivatives_without_der_x!(sol_x[i], m, sol_t[i])
    end
    m.storeResult = false

    if ismissing(algorithm)
        m.algorithmName = getAlgorithmName(solution.alg)
    end

    return solution
end


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

- `A::Matrix`: Matrix A of the linear ODE: ``\\Delta \\dot{x} = A*\\Delta x``.

- `finalStates::Vector`: Linearization point.


# Example

```julia
using Modia
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
firstOrder2 = InstantiatedModel{Measurement{Double64}}(firstOrder1)
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

function linearize!(m::InstantiatedModel{FloatType,TimeType}, algorithm=missing;
                    merge = nothing, stopTime = 0.0, analytic = false, kwargs...) where {FloatType,TimeType}
    if analytic
        @info "linearize!(.., analytic=true) of model $(m.modelName) \nis modified to analytic=false, because analytic=true is currently not supported!"
        analytic = false
    end

    solution = simulate!(m, algorithm; merge=merge, stopTime=stopTime, kwargs...)
    finalStates = solution[:,end]

    # Function that shall be linearized
    function modelToLinearize!(der_x, x)
        derivatives!(der_x, x, m, m.options.startTime)
        return nothing
    end

    # Linearize
    if analytic
        der_x = zeros(FloatType, length(finalStates))
        A = ForwardDiff.jacobian(modelToLinearize!, der_x, finalStates)
    else
        A = zeros(FloatType, length(finalStates), length(finalStates))
        FiniteDiff.finite_difference_jacobian!(A, modelToLinearize!, finalStates)
    end

    return (A, finalStates)
end



"""
    hasParameter(instantiatedModel, name::AbstractString)

Return true if parameter `name` (for example `name = "a.b.c"`)
is defined in the instantiateModel.
"""
hasParameter(m::InstantiatedModel, name::AbstractString) = begin
    if isnothing(m) || ismissing(m) || ismissing(m.result)
        return false
    end
    !ismissing(get_value(m.evaluatedParameters, name))
end


"""
    getParameter(instantiatedModel, name::AbstractString)

Return the value of parameter or init/start value `name` (for example `name = "a.b.c"`).
If `name` is not known, `missing` is returned.
"""
getParameter(m::InstantiatedModel, name::AbstractString) = get_value(m.parameters, name)


"""
    getEvaluatedParameter(instantiatedModel, name::AbstractString)

Return the value of evaluated parameter or init/start value `name` (for example `name = "a.b.c"`).
If `name` is not known, `missing` is returned.
"""
getEvaluatedParameter(m::InstantiatedModel, name::AbstractString) = get_value(m.evaluatedParameters, name)


"""
    showParameters(instantiatedModel)

Print the parameters and the init/start values.
"""
function showParameters(m::InstantiatedModel)::Nothing
    parameters = m.parameters
    @showModel parameters
    return nothing
end


"""
    showEvaluatedParameters(instantiatedModel)

Print the evaluated parameters and the evaluated init/start values.
"""
function showEvaluatedParameters(m::InstantiatedModel)::Nothing
    evaluatedParameters = m.evaluatedParameters
    @showModel evaluatedParameters
    return nothing
end


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

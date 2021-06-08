export ModiaProblem, ModiaSolve

function ModiaProblem(m1::TinyModia.SimulationModel{FloatType1,FloatType1}, 
                      m2::TinyModia.SimulationModel{FloatType2,FloatType2};
                      p, merge=nothing, kwargs...) where {FloatType1,FloatType2}

    TimeType1 = FloatType1
    TimeType2 = FloatType2
    merge2 = Map(p = parameter | convert(Vector{FloatType2}, p))
    if !isnothing(merge)
        merge2 = merge | merge2
    end
       
    # Store states x of m1 at communication points in the DifferentialEquations result data structure solution
    m1.save_x_in_solution = true 
    m2.save_x_in_solution = true

    # Save kwargs of ModiaProblem in simulation models
    empty!(m1.result)
    empty!(m2.result)    
    options1 = SimulationOptions{FloatType1,TimeType1}(merge ; kwargs...)
    options2 = SimulationOptions{FloatType2,TimeType2}(merge2; kwargs...)    
    if isnothing(options1) || isnothing(options2)
        return nothing
    end
    m1.options = options1
    m2.options = options2        

   	tspan = (options1.startTime, options1.stopTime)

    # Initialize/re-initialize SimulationModel
    if options1.log || options1.logParameters || options1.logEvaluatedParameters || options1.logStates
        println("... Construct simulation problem of DifferentialEquations.jl for ", m.modelName)
    end
     
    success1 = TinyModia.init!(m1)
    if !success1|| m1.eventHandler.restart == TinyModia.Terminate
        return nothing
    end
    
    success2 = TinyModia.init!(m2)
    if !success2|| m2.eventHandler.restart == TinyModia.Terminate
        return nothing
    end
    
#=    
    # Change parameters p in m2 to FloatType2
    m2.parameters = TinyModia.recursiveMerge(m2.parameters, Map(p = convert(Vector{FloatType2}, p)))
    m2.evaluatedParameters = TinyModia.propagateEvaluateAndInstantiate!(m2.modelModule, m2.parameters, m2.equationInfo, m2.x_start)
    if isnothing(m2.evaluatedParameters)
        return false
    end 
    
    function x_init(pp, t0::FloatType1)
        if t0 <= 0.1
            println("... init m1 called")
        end    
        success = TinyModia.init!(m1, tolerance, merge, log, logParameters, logEvaluatedParameters, logStates, logEvents)
        if success 
            return m1.x_init
        else
            return nothing
        end
    end
    
    function x_init(pp, t0::FloatType2)
        if t0 <= 0.1
            println("... init m2 called")
        end       
        success = TinyModia.init!(m2, tolerance, merge, log, logParameters, logEvaluatedParameters, logStates, logEvents)
        if success 
            return m2.x_init
        else
            return nothing
        end
    end    
=#



#    function p_derivatives!(der_x, x, p::Vector{FloatType1}, t::FloatType1)::Nothing
    function p_derivatives!(der_x, x, p, t::FloatType1)::Nothing 
        m1.evaluatedParameters.p .= p
        Base.invokelatest(m1.getDerivatives!, der_x, x, m1, t)
        return nothing
    end

    function p_derivatives!(der_x, x, p, t::FloatType2)::Nothing
        m2.evaluatedParameters.p .= p
        Base.invokelatest(m2.getDerivatives!, der_x, x, m2, t)
        return nothing
    end

    #function affect_outputs!(integrator)::Nothing
    #    if eltype(integrator.p) == FloatType1
    #        m1.storeResult = true
    #        Base.invokelatest(m1.getDerivatives!, m1.der_x, integrator.u, m1, integrator.t)
    #        m1.storeResult = false
    #    end
    #    return nothing
    #end
    #callback = DifferentialEquations.PresetTimeCallback((startTime2:interval:stopTime2), affect_outputs!)

    # Define problem and callbacks based on algorithm and model type
    abstol = 0.1*options1.tolerance
    return DifferentialEquations.ODEProblem(p_derivatives!, m1.x_init, tspan, p;   # x_init
             reltol=options1.tolerance, abstol=abstol,      # callback = callback, 
             modia_interval = options1.interval,
             modia_instantiatedModel = m1)
end


function ModiaSolve(problem, algorithm=missing; p, adaptive::Bool = true)

    startTime = problem.tspan[1]
    stopTime  = problem.tspan[2]
    interval  = problem.kwargs.data[:modia_interval]
    m         = problem.kwargs.data[:modia_instantiatedModel]
    m.algorithmType = typeof(algorithm)
    
    # Define problem and callbacks based on algorithm and model type
    if abs(interval) < abs(stopTime-startTime)
        tspan = startTime:interval:stopTime
    else
        tspan = [startTime, stopTime]
    end

    # Initial step size (the default of DifferentialEquations is too large) + step-size of fixed-step algorithm
    dt = adaptive ? interval/10 : interval    # initial step-size

    # Compute solution
    solution = ismissing(algorithm) ? DifferentialEquations.solve(problem, p=p,
                                        saveat = tspan, adaptive=adaptive, dt=dt) :
                                    DifferentialEquations.solve(problem, algorithm, p=p,
                                        saveat = tspan, adaptive=adaptive, dt=dt)
    m.solution = solution
    return solution
end

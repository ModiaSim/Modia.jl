export ModiaProblem, ModiaSolve

function ModiaProblem(m1::Modia.InstantiatedModel{FloatType1,FloatType1}, 
                      m2::Modia.InstantiatedModel{FloatType2,FloatType2};
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
    tspan_outputs = options1.startTime:options1.interval:options1.stopTime
    
    # Initialize/re-initialize InstantiatedModel
    if options1.log || options1.logParameters || options1.logEvaluatedParameters || options1.logStates
        println("... Construct simulation problem of DifferentialEquations.jl for ", m.modelName)
    end
     
    success1 = Modia.init!(m1)
    if !success1|| m1.eventHandler.restart == Modia.Terminate
        return nothing
    end
    
    success2 = Modia.init!(m2)
    if !success2|| m2.eventHandler.restart == Modia.Terminate
        return nothing
    end
    
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

    #= This call back gives an error:
    # type TrackedAffect has no field funciter
    function outputs!(x, t, integrator)::Nothing
       if eltype(integrator.p) == FloatType1
           m1.storeResult = true
           p_derivatives!(m1.der_x, x, integrator.p, t)
           m1.storeResult = false
       end
       return nothing
    end
    callback = DifferentialEquations.FunctionCallingCallback(outputs!, funcat=tspan_outputs)
    =#
    
    function affect_outputs!(integrator)::Nothing
       if eltype(integrator.p) == FloatType1
           m1.storeResult = true
           p_derivatives!(m1.der_x, integrator.u, integrator.p, integrator.t)
           m1.storeResult = false
       end
       return nothing
    end    
    callback = DifferentialEquations.PresetTimeCallback(tspan_outputs, affect_outputs!)

    # Define problem and callbacks based on algorithm and model type
    abstol = 0.1*options1.tolerance
    return DifferentialEquations.ODEProblem(p_derivatives!, m1.x_init, tspan, p;   # x_init
             reltol=options1.tolerance, abstol=abstol, callback = callback, 
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

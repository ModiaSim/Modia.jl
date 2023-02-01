# License for this file: MIT (expat)
# Copyright 2022, DLR Institute of System Dynamics and Control

newStatistics(startTime_s, algorithm::String, odeIntegrator::Bool, FloatType) = SignalTables.Map(
        initCpuTime_s      = 0.0,
        simCpuTime_s       = 0.0,
        initAlloc_MB       = 0.0,
        simAlloc_MB        = 0.0,
        FloatType          = string(FloatType),
        algorithm          = algorithm,
        odeIntegrator      = odeIntegrator,
        startTime_s        = startTime_s,
        terminationTime_s  = missing,
        nStates            = 0,
        nResults           = 0,
        nf_total           = 0,
        nf_integrator      = 0,
        nf_zeroCrossings   = 0,
        nJac               = 0,
        nAcceptedSteps     = 0,
        nRejectedSteps     = 0,
        nTimeEvents        = 0,
        nStateEvents       = 0,
        nRestartEvents     = 0,
        nFullRestartEvents = 0,
        linearSystemsSizes   = missing,
        linearSystemsRecFac  = missing,
        linearSystemsOdeMode = missing
    )

function updateStatistics_initStat!(statistics::OrderedDict{Symbol,Any}, initCpuTime_s, initAlloc_MB)::Nothing
    statistics[:initCpuTime_s] = Float64(initCpuTime_s)
    statistics[:initAlloc_MB]  = Float64(initAlloc_MB)
    return nothing
end

function updateStatistics_nStates!(statistics::OrderedDict{Symbol,Any}, nStates::Int)::Nothing
    statistics[:nStates] = nStates
    return nothing
end

function updateStatistics_nStates_segmented!(statistics::OrderedDict{Symbol,Any}, nStates::Int)::Nothing
    if typeof(statistics[:nStates]) == Int
        statistics[:nStates] = Int[statistics[:nStates], nStates]
    else
        push!(statistics[:nStates], nStates)
    end
    return nothing
end

function updateStatistics_final!(statistics::OrderedDict{Symbol,Any}, simCpuTime_s, simAlloc_MB)::Nothing
    statistics[:simCpuTime_s] = Float64(simCpuTime_s)
    statistics[:simAlloc_MB]  = Float64(simAlloc_MB)
    return nothing
end

function updateStatistics_performanceIndictors!(statistics::OrderedDict{Symbol,Any},
            terminationTime_s, nResults::Int,
            nf_total::Int, nf_integrator::Int, nf_zeroCrossings::Int,
            nJac, nAcceptedSteps, nRejectedSteps,
            nTimeEvents::Int, nStateEvents::Int, nRestartEvents::Int, nFullRestartEvents::Int,
            linearEquations)::Nothing
    if nFullRestartEvents == 0
        statistics[:terminationTime_s]    = Float64(terminationTime_s)
        statistics[:nResults]             = nResults
        statistics[:nf_total]             = nf_total
        statistics[:nf_integrator]        = nf_integrator
        statistics[:nf_zeroCrossings]     = nf_zeroCrossings
        statistics[:nJac]                 = nJac
        statistics[:nAcceptedSteps]       = nAcceptedSteps
        statistics[:nRejectedSteps]       = nRejectedSteps
        statistics[:nTimeEvents]          = nTimeEvents
        statistics[:nStateEvents]         = nStateEvents
        statistics[:nRestartEvents]       = nRestartEvents
        statistics[:linearSystemsSizes]   = Int[length(leq.b) for leq in linearEquations]
        statistics[:linearSystemsRecFac]  = Bool[leq.useRecursiveFactorization for leq in linearEquations]
        statistics[:linearSystemsOdeMode] = Bool[leq.odeMode for leq in linearEquations]

    elseif typeof(statistics[:nResults]) == Int
        statistics[:terminationTime_s]    = Float64[statistics[:terminationTime_s], Float64(terminationTime_s)]
        statistics[:nResults]             = Int[statistics[:nResults]         , nResults]
        statistics[:nf_total]             = Int[statistics[:nf_total]         , nf_total]
        statistics[:nf_integrator]        = Int[statistics[:nf_integrator]    , nf_integrator]
        statistics[:nf_zeroCrossings]     = Int[statistics[:nf_zeroCrossings] , nf_zeroCrossings]
        statistics[:nJac]                 = Union{Missing,Int}[statistics[:nJac]          , nJac]
        statistics[:nAcceptedSteps]       = Union{Missing,Int}[statistics[:nAcceptedSteps], nAcceptedSteps]
        statistics[:nRejectedSteps]       = Union{Missing,Int}[statistics[:nRejectedSteps], nRejectedSteps]
        statistics[:nTimeEvents]          = Int[statistics[:nTimeEvents]   , nTimeEvents]
        statistics[:nStateEvents]         = Int[statistics[:nStateEvents]  , nStateEvents]
        statistics[:nRestartEvents]       = Int[statistics[:nRestartEvents], nRestartEvents]
        statistics[:linearSystemsSizes]   = Vector{Int}[statistics[:linearSystemsSizes], Int[length(leq.b) for leq in linearEquations]]
        statistics[:linearSystemsRecFac]  = Vector{Int}[statistics[:linearSystemsRecFac], Bool[leq.useRecursiveFactorization for leq in linearEquations]]
        statistics[:linearSystemsOdeMode] = Vector{Int}[statistics[:linearSystemsOdeMode], Bool[leq.odeMode for leq in linearEquations]]

    else
        push!(statistics[:terminationTime_s]   , Float64(terminationTime_s))
        push!(statistics[:nResults]            , nResults)
        push!(statistics[:nf_total]            , nf_total)
        push!(statistics[:nf_integrator]       , nf_integrator)
        push!(statistics[:nf_zeroCrossings]    , nf_zeroCrossings)
        push!(statistics[:nJac]                , nJac)
        push!(statistics[:nAcceptedSteps]      , nAcceptedSteps)
        push!(statistics[:nRejectedSteps]      , nRejectedSteps)
        push!(statistics[:nTimeEvents]         , nTimeEvents)
        push!(statistics[:nStateEvents]        , nStateEvents)
        push!(statistics[:nRestartEvents]      , nRestartEvents)
        push!(statistics[:linearSystemsSizes]  , Int[length(leq.b) for leq in linearEquations])
        push!(statistics[:linearSystemsRecFac] , Bool[leq.useRecursiveFactorization for leq in linearEquations])
        push!(statistics[:linearSystemsOdeMode], Bool[leq.odeMode for leq in linearEquations])
    end
    statistics[:nFullRestartEvents] = nFullRestartEvents
    return nothing
end

function printStatistics(statistics::OrderedDict{Symbol,Any}, interval, tolerance):: Nothing
    if statistics[:nFullRestartEvents] == 0
        str_nstates = string(statistics[:nStates])
        lsSizes = statistics[:linearSystemsSizes]
        lsRec   = statistics[:linearSystemsRecFac]
        lsOde   = statistics[:linearSystemsOdeMode]
        str_sizes = "["
        str_rec   = "["
        str_ode   = "["
        for i in 1:length(lsSizes)
            if i == 1
                str_sizes = "[" * string(lsSizes[i])
                str_rec   = "[" * string(lsRec[i])
                str_ode   = "[" * string(lsOde[i])
            else
                str_sizes = str_sizes * ", " * string(lsSizes[i])
                str_rec   = str_rec   * ", " * string(lsRec[i])
                str_ode   = str_ode   * ", " * string(lsOde[i])
            end
        end
        str_sizes = str_sizes * "]"
        str_rec   = str_rec   * "]"
        str_ode   = str_ode   * "]"
    else
        nStates    = statistics[:nStates]
        nStatesMin = minimum(nStates)
        nStatesMax = maximum(nStates)
        if nStatesMin == nStatesMax
            str_nstates = string(nStatesMax)
        else
            str_nstates = string(nStatesMin) * ".." * string(nStatesMax)
        end
        lsSizes = statistics[:linearSystemsSizes]
        lsRec   = statistics[:linearSystemsRecFac]
        lsOde   = statistics[:linearSystemsOdeMode]
        str_sizes = "["
        str_ode   = "["
        str_rec   = "["
        for i in 1:length(lsSizes[1])
            smin = typemax(Int)
            smax = 0
            recmin = true
            recmax = false
            odemin = true
            odemax = false
            for j in 1:length(lsSizes)
                smin   = min(smin  , lsSizes[j][i])
                smax   = max(smax  , lsSizes[j][i])
                recmin = min(recmin, lsRec[j][i])
                recmax = max(recmax, lsRec[j][i])
                odemin = min(odemin, lsOde[j][i])
                odemax = max(odemax, lsOde[j][i])
            end
            if i > 1
                str_sizes = str_sizes * ", "
                str_rec   = str_rec   * ", "
                str_ode   = str_ode   * ", "
            end
            str_sizes = str_sizes * string(smin)
            str_rec   = str_rec   * string(recmin)
            str_ode   = str_ode   * string(odemin)
            if smax > smin
                str_sizes = str_sizes * ".." * string(smax)
            end
            if recmax > recmin
                str_rec = str_rec * ".." * string(recmax)
            end
            if odemax > odemin
                str_ode = str_ode * ".." * string(odemax)
            end
        end
        str_sizes = str_sizes * "]"
        str_rec   = str_rec   * "]"
        str_ode   = str_ode   * "]"
    end

    odeIntegrator::Bool     = statistics[:odeIntegrator]
    nFullRestartEvents::Int = statistics[:nFullRestartEvents]
    println("        initCpuTime          = ", round(statistics[:initCpuTime_s], sigdigits=3), " s")
    println("        simCpuTime           = ", round(statistics[:simCpuTime_s] , sigdigits=3), " s")
    println("        initAlloc            = ", round(statistics[:initAlloc_MB] , sigdigits=3), " MB")
    println("        simAlloc             = ", round(statistics[:simAlloc_MB]  , sigdigits=3), " MB")
    println("        FloatType            = ", statistics[:FloatType])
    println("        algorithm            = ", statistics[:algorithm], odeIntegrator ? " (ODE integrator)" : " (DAE integrator)")
    println("        startTime            = ", statistics[:startTime_s], " s")
    println("        terminationTime      = ", statistics[:terminationTime_s][end], " s")
    println("        interval             = ", interval, " s")
    println("        tolerance            = ", tolerance, " (relative tolerance)")
    println("        nStates              = ", str_nstates)
    println("        nResults             = ", sum(statistics[:nResults]))
    println("        nf_total             = ", sum(statistics[:nf_total])        , " (total number of getDerivatives! calls)")
    println("        nf_integrator        = ", sum(statistics[:nf_integrator])   , " (number of getDerivatives! calls from integrator)")  # solution.destats.nf
    println("        nf_zeroCrossings     = ", sum(statistics[:nf_zeroCrossings]), " (number of getDerivatives! calls for zero crossing detection)")
    println("        nJac                 = ", nFullRestartEvents==0 ? statistics[:nJac] : sum(statistics[:nJac]), " (number of Jacobian computations)")
    println("        nAcceptedSteps       = ", nFullRestartEvents==0 ? statistics[:nAcceptedSteps] : sum(statistics[:nAcceptedSteps]))
    println("        nRejectedSteps       = ", nFullRestartEvents==0 ? statistics[:nRejectedSteps] : sum(statistics[:nRejectedSteps]))
    println("        nTimeEvents          = ", sum(statistics[:nTimeEvents]))
    println("        nStateEvents         = ", sum(statistics[:nStateEvents]))
    println("        nRestartEvents       = ", sum(statistics[:nRestartEvents]))
    if nFullRestartEvents > 0
        println("        nFullRestartEvents   = ", nFullRestartEvents)
    end
    println("        linearSystemsSizes   = ", str_sizes)
    println("        linearSystemsRecFac  = ", str_rec, " (= true, if LU with RecursiveFactorization.jl)")
    #if !odeIntegrator
        println("        linearSystemsOdeMode = ", str_ode)
    #end
    return nothing
end

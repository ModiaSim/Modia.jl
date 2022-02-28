# License for this file: MIT (expat)
# Copyright 2017-2021, DLR Institute of System Dynamics and Control

# Provide functions to handle time and state events


"""
    @enum EventRestart NoRestart Restart FullRestart Terminate

Define how to continue or restart integration after an event. Usually, `Restart`
should be used. Only in special cases, the other flags are useful.

- `NoRestart`, continue integration without restarting the integrator
- `Restart`, restart integrator
- `FullRestart`, restart integrator and simulation state (so dimensions may change)
- `Terminate`, terminate integration
"""
@enum EventRestart NoRestart Restart FullRestart Terminate


using  ForwardDiff
getValue(v) = v
getValue(v::ForwardDiff.Dual) = v.value
getValue(v::Measurements.Measurement) = Measurements.value(v)


const nClock = 100
const nSample = 100
const nAfterDefault = 100

mutable struct EventHandler{FloatType,TimeType}
    stopTime::TimeType              # stopTime of simulation

    # Logging
    logEvents::Bool                 # = true, if events shall be logged
    nZeroCrossings::Int             # Number of zero crossing calls
    nRestartEvents::Int             # Number of Restart events
    nStateEvents::Int               # Number of state events
    nTimeEvents::Int                # Number of time events

    # Input values for the event functions
    time::TimeType                  # Current simulation time
    initial::Bool                   # = true, if model is called at initialization
                                    #         (if initial, event=true)
    terminal::Bool                  # = true, if model is called for termination (close files, streams, visualization, ...)
    event::Bool                     # = true, if model is called at an event
    afterSimulationStart::Bool      # = true, if model is called after simulation start
    crossing::Bool                  # = true, if model is called to compute crossing function
    firstEventIterationDirectlyAfterInitial::Bool # = true, at the first event iteration directly after initial at model.options.startTime.
    triggerEventDirectlyAfterInitial::Bool        # = true, if time event shall be triggered directly after initial at model.options.startTime

    # Values especially for simulations with several segments
    firstInitialOfAllSegments::Bool # = true, if model is called at initialization of the first simulation segment.
    terminalOfAllSegments::Bool     # = true, if model is called for termination at the last simulation segment.

    # Computed by the event functions
    # For time events:
    maxTime::TimeType         # Integrate at most up to maxTime
                              # (if integrateToEvent == true, maxTime = nextEventTime,
                              #  otherwise maxTime > nextEventTime)
    nextEventTime::TimeType   # Next time event instant (= realmax(Time) if no event)
    integrateToEvent::Bool    # = true, if integrator shall integrate to nextEventTime
                              # = false, if integrator can integrate beyond nextEventTime
    restart::EventRestart     # Type of integrator restart at current event
    newEventIteration::Bool   # = true, if another event iteration; = false, if no event iteration anymore
    firstEventIteration::Bool # = true, if first iteration at an event.

    # For state events:
    zEps::FloatType         # Epsilon for zero-crossing hysteresis
    nz::Int                 # Number of event indicators
    z::Vector{FloatType}    # Vector of event indicators (zero crossings). If one of z[i] passes
                            # zero, that is beforeEvent(z[i])*z[i] < 0, an event is triggered
                            # (allocated during instanciation according to nz).
    zPositive::Vector{Bool} # = true if z > 0 at the last event instant, otherwise false

    # For after time events:
    nafter::Int             # Number of after time events
    after::Vector{Bool}     # after[i] = time >= after(<variable i>)

    # For clocks (time events):
    clock::Vector{TimeType}
    sample::Vector{Any}

    function EventHandler{FloatType,TimeType}(; nz::Int=0, nAfter::Int=0, logEvents::Bool=false) where {FloatType,TimeType}
        @assert(nz >= 0)
        @assert(nAfter >= 0)
        nAfter = nAfter > 0 ? nAfter : nAfterDefault
        zEps   = FloatType(1e-10)
        new(floatmax(TimeType), logEvents, 0, 0, 0, 0, convert(TimeType,0),
            false, false, false, false, false, false, false, false, false, floatmax(TimeType), floatmax(TimeType),
            true, NoRestart, false, false, zEps, nz, ones(FloatType,nz), fill(false, nz), nAfter, fill(false,nAfter),
            fill(convert(TimeType,0),nClock), Vector{Any}(undef, nSample))
    end
end

# Default constructors
#EventHandler(           ; kwargs...)                   = EventHandler{Float64  ,Float64}(; kwargs...)
#EventHandler{FloatType}(; kwargs...) where {FloatType} = EventHandler{FloatType,Float64}(; kwargs...)


function reinitEventHandler(eh::EventHandler{FloatType,TimeType}, stopTime::TimeType, logEvents::Bool)::Nothing where {FloatType,TimeType}
    eh.logEvents      = logEvents
    eh.nZeroCrossings = 0
    eh.nRestartEvents = 0
    eh.nStateEvents   = 0
    eh.nTimeEvents    = 0

    eh.time     = convert(TimeType, 0)
    eh.initial  = false
    eh.terminal = false
    eh.event    = false
    eh.crossing = false
    eh.firstEventIterationDirectlyAfterInitial = false
    eh.triggerEventDirectlyAfterInitial = false

    eh.afterSimulationStart = false
    eh.firstInitialOfAllSegments = false
    eh.terminalOfAllSegments     = false
    eh.stopTime      = stopTime
    eh.maxTime       = floatmax(TimeType)
    eh.nextEventTime = floatmax(TimeType)
    eh.integrateToEvent = false
    eh.restart          = Restart
    eh.newEventIteration   = false
    eh.firstEventIteration = true
    eh.z .= convert(FloatType, 0)
    eh.after .= false

    return nothing
end


positiveCrossingAsString(positive::Bool) = positive ? " became > 0" : " became <= 0"
negativeCrossingAsString(negative::Bool) = negative ? " became < 0" : " became >= 0"


function initEventIteration!(h::EventHandler{FloatType,TimeType}, t::TimeType)::Nothing where {FloatType,TimeType}
    h.time          = t
    h.restart       = NoRestart
    h.maxTime       = floatmax(TimeType)
    h.nextEventTime = floatmax(TimeType)
    h.newEventIteration   = false
    h.firstEventIteration = true
    return nothing
end


function setNextEvent!(h::EventHandler{FloatType,TimeType}, nextEventTime::TimeType;
                       triggerEventDirectlyAfterInitial=false,
                       integrateToEvent::Bool=true,
                       restart::EventRestart=Restart)::Nothing where {FloatType,TimeType}
    #println("... setNextEvent!: time = ", h.time, ", nextEventTime = ", nextEventTime, ", restart = ", restart)
    if (h.event && nextEventTime > h.time) || (h.initial && nextEventTime >= h.time)
        if h.initial && triggerEventDirectlyAfterInitial && nextEventTime == h.time
            h.triggerEventDirectlyAfterInitial = true

            if nextEventTime < h.nextEventTime
                h.nextEventTime = nextEventTime
                if h.logEvents
                    println("        nextEventTime = ", round(nextEventTime, sigdigits=6), " s")
                end
            end
        else
            if integrateToEvent
                h.maxTime = min(h.maxTime, nextEventTime)
            end
            if nextEventTime < h.nextEventTime
                h.nextEventTime    = nextEventTime
                h.integrateToEvent = integrateToEvent
                if h.logEvents
                    println("        nextEventTime = ", round(nextEventTime, sigdigits=6),
                            " s, integrateToEvent = ", integrateToEvent ? "true" : "false")
                end
            end
        end
        h.restart = max(h.restart, restart)
    end
    return nothing
end


function after!(h::EventHandler{FloatType,TimeType}, nr::Int, t::Number, tAsString::String,
                leq::Union{Nothing,ModiaBase.LinearEquations{FloatType}};
                restart::EventRestart=Restart)::Bool where {FloatType,TimeType}
    # time >= t  (it is required that t is a discrete-time expression, but this cannot be checked)
    t2 = convert(TimeType,t)
    if h.initial
        h.after[nr] = h.time >= t2
        if h.logEvents
            println("        after(", tAsString, " (= ", t2, ")) became ", h.after[nr] ? "true" : "false")
        end
        if t2 > h.time
            setNextEvent!(h, t2, restart = restart)
        end

    elseif h.event
        if abs(h.time - t2) < 1E-10
            if h.logEvents && !h.after[nr]
                println("        after(", tAsString, " (= ", t2, ")) became true")
            end
            h.after[nr] = true
        elseif t2 > h.time
            setNextEvent!(h, t2, restart = restart)
            h.after[nr] = false
        else
            if h.logEvents && !h.after[nr]
                println("        after(", tAsString, " (= ", t2, ")) became true")
            end
            h.after[nr] = true
        end
    end
    return h.after[nr]
end


#positive!(h, nr, crossing, crossingAsString, leq_mode; restart=Restart) = positive!(h, nr, getValue(crossing), crossingAsString, leq_mode; restart=restart)

function positive!(h::EventHandler{FloatType,TimeType}, nr::Int, crossing, crossingAsString::String,
                   leq::Union{Nothing,ModiaBase.LinearEquations{FloatType}};
                   restart::EventRestart=Restart)::Bool where {FloatType,TimeType}
    crossing = getValue(crossing)

    if h.initial
        if !isnothing(leq) && leq.mode >= 0
            # crossing has no meaningful value (called in algebraic loop with zero or unit vectors of unknowns).
            h.zPositive[nr] = false    # set to an arbitrary value
            h.z[nr]         = -h.zEps  # set to an arbitrary value that is consistent to zPositive
            return h.zPositive[nr]
        end

        h.zPositive[nr] = crossing > convert(FloatType,0)
        if h.logEvents
            println("        ", crossingAsString, " (= ", crossing, ")", positiveCrossingAsString(h.zPositive[nr]))
        end

    elseif h.event
        if !isnothing(leq) && leq.mode >= 0
            # crossing has no meaningful value (called in algebraic loop with zero or unit vectors of unknowns).
            return h.zPositive[nr]
        end
        new_zPositive = crossing > convert(FloatType,0)
        change = (h.zPositive[nr] && !new_zPositive) || (!h.zPositive[nr] && new_zPositive)
        h.zPositive[nr] = new_zPositive

        if change
            h.restart = max(h.restart, restart)
            if h.logEvents
                println("        ", crossingAsString, " (= ", crossing, ")", positiveCrossingAsString(h.zPositive[nr]))
            end
            if !isnothing(leq) && leq.mode == -1
                # Solution of mixed linear equation system is not consistent to positve(..)
                leq.success = false
                if leq.niter > leq.niter_max
                    push!(leq.inconsistentPositive, crossingAsString)
                end
            end
        end
    end
    h.z[nr] = crossing + (h.zPositive[nr] ? h.zEps : -h.zEps)

    return h.zPositive[nr]
end


function negative!(h::EventHandler{FloatType,TimeType}, nr::Int, crossing, crossingAsString::String,
                   leq::Union{Nothing,ModiaBase.LinearEquations{FloatType}};
                   restart::EventRestart=Restart)::Bool where {FloatType,TimeType}
    crossing = getValue(crossing)

    if h.initial
        if !isnothing(leq) && leq.mode >= 0
            # crossing has no meaningful value (called in algebraic loop with zero or unit vectors of unknowns).
            h.zPositive[nr] = true    # set to an arbitrary value
            h.z[nr]         = h.zEps  # set to an arbitrary value that is consistent to zPositive
            return !h.zPositive[nr]
        end

        h.zPositive[nr] = !(crossing < convert(FloatType,0))
        if h.logEvents
            println("        ", crossingAsStringg, " (= ", crossing, ")", negativeCrossingAsString(!h.zPositive[nr]))
        end

    elseif h.event
        if !isnothing(leq) && leq.mode >= 0
            # crossing has no meaningful value (called in algebraic loop with zero or unit vectors of unknowns).
            return !h.zPositive[nr]
        end
        new_zPositive = !(crossing < convert(FloatType,0))
        change = (h.zPositive[nr] && !new_zPositive) || (!h.zPositive[nr] && new_zPositive)
        h.zPositive[nr] = new_zPositive

        if change
            h.restart = max(h.restart, restart)
            if h.logEvents
                println("        ", crossingAsStringg, " (= ", crossing, ")", negativeCrossingAsString(!h.zPositive[nr]))
            end
            if !isnothing(leq) && leq.mode == -1
                # Solution of mixed linear equation system is not consistent to negative(..)
                leq.success = false
                if leq.niter > leq.niter_max
                    push!(leq.inconsistentNegative, crossingAsString)
                end
            end
        end
    end
    h.z[nr] = crossing + (h.zPositive[nr] ? h.zEps : -h.zEps)

    return !h.zPositive[nr]
end


#=
function change!(h::EventHandler{FloatType,TimeType}, nr::Int, crossing::FloatType, crossingAsString::String,
                 leq::Union{Nothing,ModiaBase.LinearEquations{FloatType}};
                 restart::EventRestart=Restart)::Bool where {FloatType,TimeType}
    if leq_mode >= 0
        return crossing > convert(FloatType,0)
    end

    h.z[nr] = crossing
    if h.initial
        h.zPositive[nr] = crossing > convert(FloatType,0)
        if h.logEvents
            println("        ", crossingAsString, " (= ", crossing, ")", positiveCrossingAsString(h.zPositive[nr]))
        end
        return false

    elseif h.event
        new_zPositive = crossing > convert(FloatType,0)
        change = (h.zPositive[nr] && !new_zPositive) || (!h.zPositive[nr] && new_zPositive)
        h.zPositive[nr] = new_zPositive

        if change
            h.restart = max(h.restart, restart)
            if h.logEvents
                println("        ", crossingAsString, " (= ", crossing, ")",  positiveCrossingAsString(h.zPositive[nr]))
            end
            h.newEventIteration = true
            return true
        end
    end

    return false
end
=#


function edge!(h::EventHandler{FloatType,TimeType}, nr::Int, crossing, crossingAsString::String,
               leq::Union{Nothing,ModiaBase.LinearEquations{FloatType}};
               restart::EventRestart=Restart)::Bool where {FloatType,TimeType}

    if !isnothing(leq)
        @error "edge(" * crossingAsString * ") is called in a linear system of equations with\n" *
               "iteration variables $(leq.vTear_names). This is not supported."
    end

    crossing = getValue(crossing)

    if h.initial
        h.zPositive[nr] = crossing > convert(FloatType,0)

        if h.logEvents
            println("        ", crossingAsString, " (= ", crossing, ")", positiveCrossingAsString(h.zPositive[nr]))
        end

    elseif h.event
        new_zPositive = crossing > convert(FloatType,0)
        edge     = !h.zPositive[nr] && new_zPositive
        downEdge = h.zPositive[nr] && !new_zPositive
        h.zPositive[nr] = new_zPositive

        if edge
            h.restart = max(h.restart, restart)
            if h.logEvents
                println("        ", crossingAsString, " (= ", crossing, ") became > 0")
            end
            h.newEventIteration = true
            h.z[nr] = crossing + (h.zPositive[nr] ? h.zEps : -h.zEps)
            return true
        elseif downEdge
            h.restart = max(h.restart, NoRestart)
            if h.logEvents
                println("        ", crossingAsString, " (= ", crossing, ") became <= 0")
            end
        end
    end

    h.z[nr] = crossing + (h.zPositive[nr] ? h.zEps : -h.zEps)
    return false
end


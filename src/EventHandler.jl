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


mutable struct EventHandler{FloatType,TimeType}
    # Logging
    logEvents::Bool                 # = true, if events shall be logged
    nZeroCrossings::Int             # Number of zero crossing calls
    nRestartEvents::Int             # Number of Restart events
    nStateEvents::Int               # Number of state events
    
    # Input values for the event functions
    time::TimeType                  # Current simulation time                                
    initial::Bool                   # = true, if model is called at initialization
                                    #         (if initial, event=true)
    terminal::Bool                  # = true, if model is called for termination (close files, streams, visualization, ...)
    event::Bool                     # = true, if model is called at an event  
    afterSimulationStart::Bool      # = true, if model is called after simulation start
    crossing::Bool                  # = true, if model is called to compute crossing function
    
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
    nz::Int                 # Number of event indicators
    z::Vector{FloatType}    # Vector of event indicators (zero crossings). If one of z[i] passes
                            # zero, that is beforeEvent(z[i])*z[i] < 0, an event is triggered
                            # (allocated during instanciation according to nz).
    zPositive::Vector{Bool} # = true if z > 0 at the last event instant, otherwise false
     
    function EventHandler{FloatType,TimeType}(; nz::Int=0, logEvents::Bool=false) where {FloatType,TimeType} 
        @assert(nz >= 0)
        new(logEvents, 0, 0, 0, convert(TimeType,0), false, false, false, false, false, false, false, floatmax(TimeType), floatmax(TimeType),
            true, NoRestart, false, false, nz, ones(FloatType,nz), fill(false, nz))
    end
end

# Default constructors
EventHandler(           ; kwargs...)                   = EventHandler{Float64  ,Float64}(; kwargs...)
EventHandler{FloatType}(; kwargs...) where {FloatType} = EventHandler{FloatType,Float64}(; kwargs...)


function reinitEventHandler(eh::EventHandler{FloatType,TimeType}, logEvents::Bool)::Nothing where {FloatType,TimeType} 
    eh.logEvents      = logEvents
    eh.nZeroCrossings = 0
    eh.nRestartEvents = 0
    eh.nStateEvents   = 0  
    
    eh.time     = convert(TimeType, 0)
    eh.initial  = false
    eh.terminal = false
    eh.event    = false 
    eh.crossing = false
    eh.afterSimulationStart = false
    eh.firstInitialOfAllSegments = false
    eh.terminalOfAllSegments     = false 
    eh.maxTime       = floatmax(TimeType)
    eh.nextEventTime = floatmax(TimeType)
    eh.integrateToEvent = false
    eh.restart          = Restart
    eh.newEventIteration   = false
    eh.firstEventIteration = true
    eh.z .= convert(FloatType, 0)
    
    return nothing
end


positiveCrossingAsString(positive::Bool) = positive ? " became > 0" : " became <= 0"
negativeCrossingAsString(negative::Bool) = negative ? " became < 0" : " became >= 0"

   
function initEventIteration!(h::EventHandler{FloatType,TimeType}, t)::Nothing where {FloatType,TimeType}
    h.time          = t
    h.restart       = NoRestart
    h.maxTime       = floatmax(TimeType)
    h.nextEventTime = floatmax(TimeType)
    h.newEventIteration   = false
    h.firstEventIteration = true
    return nothing
end

function terminateEventIteration!(h::EventHandler{FloatType,TimeType}) where {FloatType,TimeType}
    if h.restart == Terminate || h.restart == FullRestart
        result = true
    elseif h.initial
        result = false
    else
        result = !h.newEventIteration
    end
    h.initial = false  
    h.firstInitialOfAllSegments = false          
    h.firstEventIteration = false
    h.newEventIteration   = false
    return result
end


const zEps = 1.e-10


function setNextEvent!(h::EventHandler{FloatType,TimeType}, nextEventTime::TimeType; 
                       integrateToEvent::Bool=true, 
                       restart::EventRestart=Restart) where {FloatType,TimeType}
    if (h.event && nextEventTime > h.time) || (h.initial && nextEventTime >= h.time)
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
        h.restart = max(h.restart, restart)
    end
    return nothing
end


function positive!(h::EventHandler{FloatType,TimeType}, nr::Int, crossing::FloatType, crossingAsString::String, leq_mode::Int; 
                   restart::EventRestart=Restart)::Bool where {FloatType,TimeType}
    leq_mode = -1
    if leq_mode >= 0
        return crossing > convert(FloatType,0)
        
    elseif h.initial
        h.zPositive[nr] = crossing > convert(FloatType,0)
        if h.logEvents
            println("        ", crossingAsString, " (= ", crossing, ")", positiveCrossingAsString(h.zPositive[nr]))
        end

    elseif h.event
        # println("... nr = ", nr, ", crossing = ", crossing)
        new_zPositive = crossing > convert(FloatType,0)
        change = (h.zPositive[nr] && !new_zPositive) || (!h.zPositive[nr] && new_zPositive)
        h.zPositive[nr] = new_zPositive
        
        if change
            h.restart = max(h.restart, restart)
            if h.logEvents
                println("        ", crossingAsString, " (= ", crossing, ")", positiveCrossingAsString(h.zPositive[nr]))
            end
            h.newEventIteration = true
        end
    end
    h.z[nr] = crossing + (h.zPositive[nr] ? zEps : -zEps)
   
    return h.zPositive[nr]
end


function negative!(h::EventHandler{FloatType,TimeType}, nr::Int, crossing::FloatType, crossingAsString::String, leq_mode::Int; 
                   restart::EventRestart=Restart)::Bool where {FloatType,TimeType}
                   
    if leq_mode >= 0
        return crossing >= convert(FloatType,0)
                           
    elseif h.initial
        h.zPositive[nr] = crossing >= convert(FloatType,0)
        if h.logEvents
            println("        ", crossingAsStringg, " (= ", crossing, ")", negativeCrossingAsString(!h.zPositive[nr]))
        end

    elseif h.event
        new_zPositive = crossing >= convert(FloatType,0)
        change = (h.zPositive[nr] && !new_zPositive) || (!h.zPositive[nr] && new_zPositive)
        h.zPositive[nr] = new_zPositive
        
        if change
            h.restart = max(h.restart, restart)
            if h.logEvents
                println("        ", crossingAsStringg, " (= ", crossing, ")", negativeCrossingAsString(!h.zPositive[nr]))
            end
            h.newEventIteration = true
        end
    end
    h.z[nr] = crossing + (h.zPositive[nr] ? zEps : -zEps)
  
    return !h.zPositive[nr]
end


function change!(h::EventHandler{FloatType,TimeType}, nr::Int, crossing::FloatType, crossingAsString::String, leq_mode::Int; 
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


function edge!(h::EventHandler{FloatType,TimeType}, nr::Int, crossing::FloatType, crossingAsString::String, leq_mode::Int; 
               restart::EventRestart=Restart)::Bool where {FloatType,TimeType}
    
    if leq_mode >= 0
        return false
    end    
    
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
            h.z[nr] = crossing + (h.zPositive[nr] ? zEps : -zEps)                 
            return true
        elseif downEdge 
            h.restart = max(h.restart, NoRestart)
            if h.logEvents
                println("        ", crossingAsString, " (= ", crossing, ") became <= 0")
            end
        end        
    end
        
    h.z[nr] = crossing + (h.zPositive[nr] ? zEps : -zEps) 
    return false
end


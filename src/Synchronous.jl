
#= Synchronous Modelica primitives

Clock() 
Clock(intervalCounter::Int, resolution::Int=1)  # Add version with rational number
Clock(interval::Float64) 
Clock(condition::Bool, startInterval::Float64=0.0) 
Clock(c::Clock, solverMethod) 

previous(x)

sample(u, c::Clock)
sample(u) 

hold(u) = u

subSample(u, factor::Int)
superSample(u, factor::Int)
shiftSample(u, shiftCounter::Int, resolution::Int=1)  # Add version with rational number
backSample(u, backCounter::Int, resolution::Int=1)  # Add version with rational number

noClock() 

interval(u)

=#

export Clock, sample, hold, previous


Clock(interval, m::InstantiatedModel, nr::Int) = Clock(m.options.startTime, interval, m, nr)


function Clock(startTime, interval, m::InstantiatedModel{FloatType,TimeType}, nr::Int)::Bool where {FloatType,TimeType}
# Clock ticks at startTime, startTime + interval, startTime + 2*interval, ...
# This is independent of the starTime of the simulation. Example:
#   if startTime = -3, interval = 1.5, m.options.startTime = 1.0 then
#   the clock ticks at 1.5, 3.0, 4.5, ...
    eh = m.eventHandler
    if isEvent(m)
        if isInitial(m)
            startTime2 = convert(TimeType, startTime)
            interval2  = convert(TimeType, interval)
            
            tFirst = startTime2 >= eh.time ? startTime2 : startTime2 + div(eh.time - startTime2, interval2, RoundUp)*interval2
            if abs(eh.time - tFirst) < 1e-10
                eh.clock[nr] = eh.time     
                setNextEvent!(eh, eh.clock[nr], triggerEventDirectlyAfterInitial=true)
            else    
                eh.clock[nr] = tFirst 
                setNextEvent!(eh, eh.clock[nr])
            end
         
        elseif isFirstEventIterationDirectlyAfterInitial(m)
            eh.clock[nr] = eh.time + convert(TimeType,interval)            
            setNextEvent!(eh, eh.clock[nr])
            return true
            
        elseif isFirstEventIteration(m) && isAfterSimulationStart(m)
            tick = abs(eh.time - eh.clock[nr]) < 1E-10
            if tick
                eh.clock[nr] = eh.time + convert(TimeType,interval)            
                setNextEvent!(eh, eh.clock[nr])
                return true
            end
            setNextEvent!(eh, eh.clock[nr])
        end
    end

    return false
end


@inline function sample(v, clock::Bool, m::InstantiatedModel{FloatType,TimeType}, nr::Int) where {FloatType,TimeType}
    eh = m.eventHandler
    if isInitial(m) || clock
        eh.sample[nr] = v
    end
    return eh.sample[nr] 
end


@inline function previous(clock::Bool, m::InstantiatedModel, nr::Int)
    # m.previous[nr] is initialized with the start/init value of v before the first model evaluation
    if clock 
        m.previous[nr] = m.nextPrevious[nr]
    end
    return m.previous[nr]
end


hold(v) = v


@inline function hold(v, clock::Bool, m::InstantiatedModel, nr::Int)
    # m.hold[nr] is initialized with the start/init value of v before the first model evaluation
    if clock 
        m.hold[nr] = v
    end
    return m.hold[nr]
end


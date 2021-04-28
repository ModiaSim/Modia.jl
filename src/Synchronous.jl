
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

# The argument list of function F is extended with the argument simulationModel.

export Clock, @Clock, sample, hold, previous, updatePrevious

# -------------------------------------------------------------------

macro Clock(interval, nr)
    code = :(Clock($interval, instantiatedModel, $nr))
    return esc(code)
end


function Clock(interval, m::SimulationModel{FloatType,TimeType}, nr::Int) where {FloatType,TimeType}
    eh = m.eventHandler
    if isInitial(m)
        #if eh.logEvents
        #    println("        in Clock, nr = $nr (isInitial)")
        #end

        eh.clock[nr] = eh.time + convert(TimeType,interval)     
        setNextEvent!(eh, eh.clock[nr])
    
    elseif isEvent(m) && isAfterSimulationStart(m)
        tick = abs(eh.time - eh.clock[nr]) < 1E-10
        if tick
            #if eh.logEvents
            #    println("        in Clock, nr = $nr (isEvent; clock is active)")
            #end
            eh.clock[nr] = eh.time + convert(TimeType,interval)            
            setNextEvent!(eh, eh.clock[nr])
            return true
        end
        setNextEvent!(eh, eh.clock[nr])
    end

    return false
end


function sample(v, clock::Bool, m::SimulationModel{FloatType,TimeType}, nr::Int) where {FloatType,TimeType}
    eh = m.eventHandler
    if isInitial(m)
        #if eh.logEvents
        #    println("        in sample, nr = $nr (isInitial)")
        #end

        eh.sample[nr] = v
        return v
    end
  
    if clock 
        #if eh.logEvents
        #    println("        in sample, nr = $nr (clock is active)")
        #end
        eh.sample[nr] = v
        return v
    else
        return eh.sample[nr] 
    end
end

# initPrevious not used
#=
function initPrevious(v, m::SimulationModel{FloatType,TimeType}, nr::Int) where {FloatType,TimeType}
    eh = m.eventHandler
    if isInitial(m)
        #if eh.logEvents
        #    println("        in initPrevious, nr = $nr (initialize previous store)")
        #end

        eh.previous[nr] = v
        return v
    else
        return eh.previous[nr]   
    end
end
=#

function previous(v, clock::Bool, m::SimulationModel{FloatType,TimeType}, nr::Int) where {FloatType,TimeType}
    eh = m.eventHandler
    if isInitial(m)
        #if eh.logEvents
        #    println("        in previous, nr = $nr (initialize previous store)")
        #end 

        eh.previous[nr]     = v
        eh.nextPrevious[nr] = v
        return v
    end
  
    if clock 
        eh.previous[nr] = eh.nextPrevious[nr];
        #if eh.logEvents
        #    println("        in previous, nr = $nr (clock is active)")
        #end
    end
    return eh.previous[nr]
end

#=
function previous(v, clock::Bool)
  return v
end
=#

function updatePrevious(v, m::SimulationModel{FloatType,TimeType}, nr::Int) where {FloatType,TimeType}
    eh = m.eventHandler
    if isEvent(m) && isAfterSimulationStart(m)
    # if eh.logEvents
    #   println("        in updatePrevious, nr = $nr, v = $v, time = ", eh.time)
    # end
        eh.nextPrevious[nr] = v
    end
end

hold(v) = v

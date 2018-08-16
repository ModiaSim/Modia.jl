module Synchronous

using ..Instantiation

import ModiaMath 
import ModiaMath.ModiaToModiaMath  # ModiaSimulationModel, EventRestart, Restart, positive!, ...
import ModiaMath.ModiaToModiaMath.ModiaSimulationModel
import ModiaMath.Restart

@static if VERSION < v"0.7.0-DEV.2005"
    Nothing = Void 
end

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

const nClock = 100
const nSample = 100
const nPrevious = 100
const nCrossings = 100

mutable struct Store
    clock::Vector{ModiaMath.Time}
    sample::Vector{Any}
    previous::Vector{Any}
    nextPrevious::Vector{Any}
    crossings::Vector{Any}
  
    @static if VERSION < v"0.7.0-DEV.2005"
      Store() = new(zeros(nClock), Array{Any}(nSample), zeros(nPrevious), zeros(nPrevious), zeros(nCrossings))
    else
      Store() = new(zeros(nClock), Array{Any}(undef, nSample), zeros(nPrevious), zeros(nPrevious), zeros(nCrossings))
    end    
end

isLog(m) = true

# -------------------------------------------------------------------

# Temporary:
allInstances(s) = return s


function Clock(interval::ModiaMath.Time, m::ModiaSimulationModel, nr::Int) 
    if ModiaToModiaMath.isPreInitial(m)
        if ModiaToModiaMath.isLogEvents(m)
            println("        in Clock, nr = $nr (isPreInitial)")
        end
        return false
    end
  
    if ModiaToModiaMath.isInitial(m)
        if ModiaToModiaMath.isLogEvents(m)
            println("        in Clock, nr = $nr (isInitial)")
        end

        if typeof(m.store) == Nothing
            m.store = Store()
        end

        m.store.clock[nr] = ModiaMath.getTime(m)
        ModiaToModiaMath.setNextEvent!(m, m.store.clock[nr])
    
    elseif ModiaToModiaMath.isEvent(m) && ModiaToModiaMath.isAfterSimulationStart(m)
        tick = abs(ModiaMath.getTime(m) - m.store.clock[nr]) < 1E-10
        if tick
            if ModiaToModiaMath.isLogEvents(m)
                println("        in Clock, nr = $nr (isEvent; clock is active)")
            end
        
            m.store.clock[nr] = ModiaMath.getTime(m) + interval
            ModiaToModiaMath.setNextEvent!(m, m.store.clock[nr])
            return true
        end
        ModiaToModiaMath.setNextEvent!(m, m.store.clock[nr])
    end

    return false
end

function positive(crossing::Float64, m::ModiaSimulationModel, nr::Int; restart::ModiaToModiaMath.EventRestart=Restart)
    if typeof(m.store) == Nothing
        m.store = Store()
    end
        println("positive:"); @show nr m; ModiaToModiaMath.positive!(crossing, m, nr; restart=restart)
end

positiveChange(crossing::Float64, m::ModiaSimulationModel, nr::Int; restart::ModiaToModiaMath.EventRestart=Restart) =
         ModiaToModiaMath.positiveChange!(crossing, m, nr; restart=restart)
positiveEdge(crossing::Float64, m::ModiaSimulationModel, nr::Int; restart::ModiaToModiaMath.EventRestart=Restart) =
         ModiaToModiaMath.positiveEdge!(crossing, m, nr; restart=restart)

function sample(v, clock::Bool, m::ModiaSimulationModel, nr::Int) 
    if ModiaToModiaMath.isPreInitial(m)
        if ModiaToModiaMath.isLogEvents(m)
            println("        in sample, nr = $nr (isPreInitial)")
        end
        return v
    elseif ModiaToModiaMath.isInitial(m)
        if ModiaToModiaMath.isLogEvents(m)
            println("        in sample, nr = $nr (initialize sample store)")
        end

        if typeof(m.store) == Nothing
            m.store = Store()
        end 

        m.store.sample[nr] = v
        return v
    end
  
    if clock 
        if ModiaToModiaMath.isLogEvents(m)
            println("        in sample, nr = $nr (clock is active)")
        end
        m.store.sample[nr] = v;
        return v
    else
        return m.store.sample[nr] 
    end
end

# initPrevious not used
function initPrevious(v, m::ModiaSimulationModel, nr::Int) 
    if ModiaToModiaMath.isPreInitial(m)
        if ModiaToModiaMath.isLogEvents(m)
            println("        in initPrevious, nr = $nr (isPreInitial)")
        end
        return v

    elseif ModiaToModiaMath.isInitial(m)
        if ModiaToModiaMath.isLogEvents(m)
            println("        in initPrevious, nr = $nr (initialize previous store)")
        end

        if typeof(m.store) == Nothing
            m.store = Store()
        end

        m.store.previous[nr] = v
        return v
    else
        return m.store.previous[nr]   
    end
end

function previous(v, clock::Bool, m::ModiaSimulationModel, nr::Int) 
    if ModiaToModiaMath.isPreInitial(m)
        if ModiaToModiaMath.isLogEvents(m)
            println("        in previous, nr = $nr (isPreInitial)")
        end
        return v

    elseif ModiaToModiaMath.isInitial(m)
        if ModiaToModiaMath.isLogEvents(m)
            println("        in previous, nr = $nr (initialize previous store)")
        end

        if typeof(m.store) == Nothing
            m.store = Store()
        end      

        m.store.previous[nr]     = v
        m.store.nextPrevious[nr] = v
        return v
    end
  
    if clock 
        m.store.previous[nr] = m.store.nextPrevious[nr];
        if ModiaToModiaMath.isLogEvents(m)
            println("        in previous, nr = $nr (clock is active)")
        end
    end
    return m.store.previous[nr]
end

#=
function previous(v, clock::Bool)
  return v
end
=#

function updatePrevious(v, m::ModiaSimulationModel, nr::Int) 
    if ModiaToModiaMath.isInitial(m)
        if typeof(m.store) == Nothing
            m.store = Store()
        end 
    end

    if ModiaToModiaMath.isEvent(m) && ModiaToModiaMath.isAfterSimulationStart(m)
    # if ModiaToModiaMath.isLogEvents(m)
    #   println("        in updatePrevious, nr = $nr, v = $v, time = ", ModiaMath.getTime(m))
    # end
        m.store.nextPrevious[nr] = v
    end
end

hold(v) = v

end
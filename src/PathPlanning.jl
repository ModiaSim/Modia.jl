# License for this file: MIT (expat)
# Copyright 2017-2021, DLR Institute of System Dynamics and Control

export PTP_path, pathEndTime, getPosition!, getPosition, getIndex, getPath, plotPath

using OrderedCollections
using Unitful


"""
    path = PTP_path(names; positions, v_max, a_max, startTime=0.0)

Generate a new path object to move as fast as possible from
positions[i,:] to positions[i+1,:]. The `positions[i,:]` can be a set of translational
positions in [m], that is absolute distances, and/or rotational positions in 
[rad] that is angles. In robotics such a movement is called PTP (Point-To-Point).
The signals are constructed in such a way that it is not possible
to move faster, given the maximally allowed velocity `v_max[j]` and
the maximally allowed acceleration `a_max[j]` for signal `names[j]`
and have a velocity of zero at the given `positions`.

If there are two or more signals (that is length(names) > 1) then
the path is constructed such that all signals
are in the same periods in the acceleration, constant velocity
and deceleration phase. This means that only one of the signals
is at its limits whereas the others are synchronized in such a way
that the end point is reached at the same time instant.

For example, this means that the signals have a velocity of zero at positions[1,:],
one of the signals is accelerated with its maximally allowed acceleration until
one of the signals reaches its maximally allowed velocity. At a proper time instant,
one of the signals is decelerated with the negative value of its maximally allowed
acceleration, so that all signals reach positions[2,:] with velocity zero.

This element is useful to generate a reference signal for a controller
which controls, e.g., a drive train, or to drive
a flange according to a given acceleration.


# Example
```julia
using Modia
@usingModiaPlot

const ptp_path = PTP_path(["angle1", "angle2", "angle3"],
                          positions = [0.0 2.0 3.0;  # angle1=0.0, angle2=2.0, angle3=3.0
                                       0.5 3.0 4.0;
                                       0.8 1.5 0.3;
                                       0.2 1.5 0.8],
                          startTime = 0.1,
                          v_max = 2*ones(3),
                          a_max = 3*ones(3))
angles = zeros(3)
getPosition!(ptp_path, 0.5, angles)   # angles = [0.12, 2.24, 3.24]
path = getPath(ptp_path)
plot(path, ["angle1", "angle2", "angle3"])
```
"""
mutable struct PTP_path{FloatType}
   names::Vector{AbstractString}
   startTime::FloatType
   v_max::Vector{FloatType}
   a_max::Vector{FloatType}
   positions::Matrix{FloatType}

   delta::Matrix{FloatType}
   hasPath::Vector{Bool}
   sd_max::Vector{FloatType}
   sdd_max::Vector{FloatType}
   Ta1::Vector{FloatType}
   Ta2::Vector{FloatType}
   noWphase::Vector{Bool}
   Tv::Vector{FloatType}
   Te::Vector{FloatType}
   Ta1s::Vector{FloatType}
   Ta2s::Vector{FloatType}
   Tvs::Vector{FloatType}
   Tes::Vector{FloatType}
   sd_max2::Vector{FloatType}
   s1::Vector{FloatType}
   s2::Vector{FloatType}
   s3::Vector{FloatType}
   Tend::FloatType
   posTemp::Vector{FloatType}   # Temporary storage that can be used by the functions operation on PTP_path

   function PTP_path{FloatType}(
                     names::AbstractVector;
                     positions::Matrix{FloatType}, 
                     v_max::Vector{FloatType},
                     a_max::Vector{FloatType},
                     startTime =0.0) where {FloatType}
        startTime = convertTimeVariable(FloatType, startTime) 
        
        #@assert(size(positions,1) > 1)
        #@assert(size(positions,2) == Base.length(names))
        #@assert(Base.length(v_max) == Base.length(names))
        #@assert(Base.length(a_max) == Base.length(names))
        np    = Base.length(names)
        npath = size(positions,1) - 1   # number of path points

        for i in eachindex(v_max)
            @assert(v_max[i] > convert(FloatType,0))
            @assert(a_max[i] > convert(FloatType,0))
        end

        delta    = zeros(FloatType, npath, np)
        hasPath  = fill(false,npath)
        sd_max   = zeros(FloatType,npath)
        sdd_max  = zeros(FloatType,npath)
        Ta1      = zeros(FloatType,npath)
        Ta2      = zeros(FloatType,npath)
        noWphase = fill(false,npath)
        Tv       = zeros(FloatType,npath)
        Te       = zeros(FloatType,npath)
        Ta1s     = zeros(FloatType,npath)
        Ta2s     = zeros(FloatType,npath)
        Tvs      = zeros(FloatType,npath)
        Tes      = zeros(FloatType,npath)
        sd_max2  = zeros(FloatType,npath)
        s1       = zeros(FloatType,npath)
        s2       = zeros(FloatType,npath)
        s3       = zeros(FloatType,npath)
        aux1     = zeros(FloatType,np)
        aux2     = zeros(FloatType,np)
        small    = 1000*eps(FloatType)

        for i in 1:npath
            delta[i,:] = positions[i+1,:] - positions[i,:]

            for j in 1:np
                aux1[j] = delta[i,j]/v_max[j]
                aux2[j] = delta[i,j]/a_max[j]
            end

            sd_max_inv  = maximum(abs.(aux1))
            sdd_max_inv = maximum(abs.(aux2))

            if sd_max_inv <= small || sdd_max_inv <= small
                hasPath[i] = false
            else
                hasPath[i] = true
                sd_max[i]  = 1/sd_max_inv
                sdd_max[i] = 1/sdd_max_inv

                Ta1[i] = sqrt(1/sdd_max[i])
                Ta2[i] = sd_max[i]/sdd_max[i]
                noWphase[i] = Ta2[i] >= Ta1[i]

                Tv[i]   = noWphase[i] ? Ta1[i] : 1/sd_max[i]
                Te[i]   = noWphase[i] ? Ta1[i] + Ta1[i] : Tv[i] + Ta2[i]

                Tbegin  = i==1 ? startTime : Tes[i-1]
                Ta1s[i] = Ta1[i] + Tbegin
                Ta2s[i] = Ta2[i] + Tbegin
                Tvs[i]  = Tv[i]  + Tbegin
                Tes[i]  = Te[i]  + Tbegin
                sd_max2[i] = sdd_max[i]*Ta1[i]
                s1[i] = sdd_max[i]*(noWphase[i] ? Ta1[i]*Ta1[i] : Ta2[i]*Ta2[i])/2
                s2[i] = s1[i] + (noWphase[i] ? sd_max2[i]*(Te[i] - Ta1[i]) - (sdd_max[i]/2)*(Te[i] - Ta1[i])^2 : sd_max[i]*(Tv[i] - Ta2[i]))
                s3[i] = s2[i] + sd_max[i]*(Te[i] - Tv[i]) - (sdd_max[i]/2)*(Te[i] - Tv[i])^2
            end
        end

        Tend = Tes[end]

        new(names, startTime, v_max, a_max, positions, delta,
            hasPath, sd_max, sdd_max, Ta1, Ta2, noWphase, Tv, Te, Ta1s, Ta2s, Tvs, Tes, sd_max2, s1, s2, s3, Tend, zeros(FloatType,np))
    end
end

PTP_path(args...; kwargs...) = PTP_path{Float64}(args...; kwargs...)

                     

"""
    Tend = pathEndTime(path)

Given a `path::PTP_path` return the end time `Tend` of the path.
"""
pathEndTime(path::PTP_path) = path.Tend



"""
    getPosition!(path, time, position)

Given a `path::PTP_path` and a time instant `time`, return the actual
position at time `time` in vector `position`.
"""
function getPosition!(path::PTP_path{FloatType}, time, position::Vector{FloatType}) where {FloatType}
    time = convertTimeVariable(FloatType, time)
    @assert(length(position) == size(path.positions,2))
    npath = length(path.hasPath)
    np    = length(position)

    # Search correct time interval
    i = 0
    if time <= path.startTime
        i = 1
        s = 0

    else
        while i < npath
            i = i+1
            if time <= path.Tes[i]
               break
            end
        end
        if time >= path.Tes[end]
            i = npath
            s = path.noWphase[i] ? path.s2[end] : path.s3[end]
            #println("... time=$time i=$i s=$s qbegin=", path.positions[i,1], ", qdelta = ", path.delta[i,1])
        else
            Tbegin = i==1 ? path.startTime : path.Tes[i-1]
            if path.noWphase[i]
                if time < path.Ta1s[i]
                    s = (path.sdd_max[i]/2)*(time - Tbegin)^2
                elseif time < path.Tes[i]
                    s = path.s1[i] + path.sd_max2[i]*(time - path.Ta1s[i]) - (path.sdd_max[i]/2)*(time - path.Ta1s[i])^2
                else
                    s = path.s2[i]
                end
            elseif time < path.Ta2s[i]
                s = (path.sdd_max[i]/2)*(time - Tbegin)^2
            elseif time < path.Tvs[i]
                s = path.s1[i] + path.sd_max[i]*(time - path.Ta2s[i])
            elseif time < path.Tes[i]
                s = path.s2[i] + path.sd_max[i]*(time - path.Tvs[i]) - (path.sdd_max[i]/2)*(time - path.Tvs[i])*(time - path.Tvs[i])
            else
                s = path.s3[i]
            end
        end
    end

    for j in 1:np
        # println("... i = $i, j= $j, np = $np, s=$s, positions[i,1]=",path.positions[i,1])
        position[j]= path.positions[i,j] + path.delta[i,j]*s
    end
end


"""
    getPosition!(path, time, position, velocity, acceleration)

Given a `path::PTP_path` and a time instant `time`, return the actual
position, velocity and acceleration at time `time` in vectors
`position, velocity, acceleration`.
"""
function getPosition!(path::PTP_path{FloatType}, time,
                      position::Vector{FloatType},
                      velocity::Vector{FloatType},
                      acceleration::Vector{FloatType})::Nothing where {FloatType}
    time = convertTimeVariable(FloatType, time)                    
    @assert(length(position)     == size(path.positions,2))
    @assert(length(velocity)     == length(position))
    @assert(length(acceleration) == length(acceleration))
    npath = length(path.hasPath)
    np    = length(position)
    s   = convert(FloatType, 0)
    sd  = convert(FloatType, 0)
    sdd = convert(FloatType, 0)

    # Search correct time interval
    i = 0
    if time <= path.startTime
        i   = 1
        s   = convert(FloatType, 0)
        sd  = convert(FloatType, 0)
        sdd = convert(FloatType, 0)
    else
        while i < npath
            i = i+1
            if time <= path.Tes[i]
               break
            end
        end
        if time >= path.Tes[end]
            i   = npath
            s   = path.noWphase[i] ? path.s2[end] : path.s3[end]
            sd  = convert(FloatType, 0)
            sdd = convert(FloatType, 0)
            #println("... time=$time i=$i s=$s qbegin=", path.positions[i,1], ", qdelta = ", path.delta[i,1])
        else
            Tbegin = i==1 ? path.startTime : path.Tes[i-1]
            if path.noWphase[i]
                if time < path.Ta1s[i]
                    s   = (path.sdd_max[i]/2)*(time - Tbegin)^2
                    sd  = path.sdd_max[i]*(time - Tbegin)
                    sdd = path.sdd_max[i]
                elseif time < path.Tes[i]
                    s   = path.s1[i] + path.sd_max2[i]*(time - path.Ta1s[i]) - (path.sdd_max[i]/2)*(time - path.Ta1s[i])^2
                    sd  = path.sd_max2[i] - path.sdd_max[i]*(time - path.Ta1s[i])
                    sdd = -path.sdd_max[i]
                else
                    s   = path.s2[i]
                    sd  = convert(FloatType, 0)
                    sdd = convert(FloatType, 0)
                end
            elseif time < path.Ta2s[i]
                s   = (path.sdd_max[i]/2)*(time - Tbegin)^2
                sd  = path.sdd_max[i]*(time - Tbegin)
                sdd = path.sdd_max[i]
            elseif time < path.Tvs[i]
                s   = path.s1[i] + path.sd_max[i]*(time - path.Ta2s[i])
                sd  = path.sd_max[i]
                sdd = convert(FloatType, 0)
            elseif time < path.Tes[i]
                s   = path.s2[i] + path.sd_max[i]*(time - path.Tvs[i]) - (path.sdd_max[i]/2)*(time - path.Tvs[i])^2
                sd  = path.sd_max[i] - path.sdd_max[i]*(time - path.Tvs[i])
                sdd = -path.sdd_max[i]
            else
                s   = path.s3[i]
                sd  = convert(FloatType, 0)
                sdd = convert(FloatType, 0)
            end
        end
    end

    for j in 1:np
        # println("... i = $i, j= $j, np = $np, s=$s, positions[i,1]=",path.positions[i,1])
        position[j]     = path.positions[i,j] + path.delta[i,j]*s
        velocity[j]     = path.delta[i,j]*sd
        acceleration[j] = path.delta[i,j]*sdd
    end

    return nothing
end


"""
    pos = getPosition(path, index, time)

Given a `path::PTP_path`, the `index` of a signal, and a time instant `time`, return the actual
position at time `time`.
"""
function getPosition(path::PTP_path, index, time)
    getPosition!(path, time, path.posTemp)
    return path.posTemp[index]
end



"""
    index = getIndex(path, name)

Return the index of `name` in `path` or trigger an error, if not present.
"""
function getIndex(path, name)
    index = findfirst(x -> x==name, path.names)
    if typeof(index) == Nothing
        error("getIndex(path,name): \"", name, "\" not in path")
    end
    return index
end


"""
    getPath(path;
            names=path.names, tend=1.1*path.Tend, ntime=101, onlyPositions=false)

Given a `path::PTP_path`, return a SignalTables.SignalTable with the time series
of the path over `time` up to `tend` for all `ntime` time points.
"""
function getPath(path::PTP_path{FloatType}; names=path.names, 
                 ntime=101, tend = 1.1*path.Tend, onlyPositions=false) where {FloatType}
    tend = convertTimeVariable(FloatType, tend)                
    tvec = range(0.0, tend, length=ntime)
    @show tvec
    indices = indexin(names, path.names)
    names2  = deepcopy(names)
    for i in eachindex(indices)
        if typeof(i) == Nothing
            @warn "getPath(path, ...): \""*names[i]*"\" is ignored, because not in path"
            deleteat!(indices,i)
            deleteat!(names2,i)
        end
    end

    np   = length(indices)
    q    = zeros(FloatType,length(tvec), np)
    qt   = zeros(FloatType,length(path.names))

    series = SignalTable("time" => Var(values=tvec, unit="s", independent=true))

    if onlyPositions
        for i in eachindex(tvec)
            getPosition!(path, tvec[i], qt)
            q[i,:] = qt[indices]
        end

        for i in eachindex(names2)
            series[names2[i]] = Var(values = q[:,i])
        end
        
    else
        der_names2  = "der(" .* names2 .* ")"
        der2_names2 = "der2(" .* names2 .* ")"
        qd   = zeros(FloatType,length(tvec), np)
        qdd  = zeros(FloatType,length(tvec), np)
        qtd  = zeros(FloatType,length(path.names))
        qtdd = zeros(FloatType,length(path.names))
        for i in eachindex(tvec)
            getPosition!(path, stripUnit(tvec[i]), qt, qtd, qtdd)
            q[i,:]   = qt[indices]
            qd[i,:]  = qtd[indices]
            qdd[i,:] = qtdd[indices]
        end

        for i in eachindex(names2)
            series[names2[i]]      = Var(values = q[:,i])
            series[der_names2[i]]  = Var(values = qd[:,i])
            series[der2_names2[i]] = Var(values = qdd[:,i])
        end
    end

    return series
end


"""
    plotPath(path, plot::Function;
             names=path.names, heading="PTP plots",
             tend=1.1*path.Tend, figure=1, ntime=101, onlyPositions=true)

Given a `path::PTP_path`, plot the path over `time` up to `tend` for all points
identified by the vector or tuple `names` to figure `figure`
using `ntime` time points. 

# Example

```julia
using Modia
@usingModiaPlot

const ptp_path = PTP_path(["angle1", "angle2", "angle3"],
                          positions = [0.0 2.0 3.0;  # angle1=0.0, angle2=2.0, angle3=3.0
                                       0.5 3.0 4.0;
                                       0.8 1.5 0.3;
                                       0.2 1.5 0.8],
                          startTime = 0.1,
                          v_max = 2*ones(3),
                          a_max = 3*ones(3))
angles = zeros(3)
getPosition!(ptp_path, 0.5, angles)   # angles = [0.12, 2.24, 3.24]
plotPath(ptp_path, plot)   # use plot(..) defined with @usingModiaPlot
```
"""
function plotPath(path::PTP_path{FloatType}, plot::Function; names=path.names, heading="PTP plots", figure=1,
                  ntime=101, tend = 1.1*path.Tend, onlyPositions=true)::Nothing where {FloatType}
    time = range(0u"s",(tend)u"s",length=ntime)
    indices = indexin(names, path.names)
    names2  = deepcopy(names)
    for i in eachindex(indices)
        if isnothing(i)
            @warn "plotPath(path, ...): \""*names[i]*"\" is ignored, because not in path"
            deleteat!(indices,i)
            deleteat!(names2,i)
        end
    end

    pathSignalTable = getPath(path; names=names, tend=tend, ntime=ntime, onlyPositions=onlyPositions)
    if onlyPositions
        plot(pathSignalTable, Tuple(names2), heading=heading, figure=figure)
    else
        der_names2  = "der(" .* names2 .* ")"
        der2_names2 = "der2(" .* names2 .* ")"
        plot(pathSignalTable, [Tuple(names2), Tuple(der_names2), Tuple(der2_names2)],
             heading=heading, figure=figure)
    end

    return nothing
end
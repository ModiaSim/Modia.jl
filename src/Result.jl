# License for this file: MIT (expat)
# Copyright 2022, DLR Institute of System Dynamics and Control
#
# Modia result datastrucure and functions operating on results.

using  OrderedCollections: OrderedDict, OrderedSet
import ModiaResult

"""
    @enum ResultKind RESULT_ELIMINATED RESULT_CONSTANT_INVARIANT RESULT_CONSTANT_SEGMENT RESULT_T RESULT_X RESULT_DER_X RESULT_W_INVARIANT RESULT_W_SEGMENT

Kind of result variable.

| ResultKind value   | Description                                                             |
|:-------------------|:------------------------------------------------------------------------|
| RESULT_ELIMINATED  | Variable is eliminated - alias info stored in result.info               |
| RESULT_CONSTANT    | Variable is constant all the time - value stored in result.info         |
| RESULT_T           | Independent variable (usually time) - stored in result.t                |
| RESULT_X           | State (invariant/segment variable) - stored in result.x                 |
| RESULT_DER_X       | State derivative (invariant/segment variable) - stored in result.der_x  |
| RESULT_W_INVARIANT | Invariant algebraic variable - stored in result.w_invariant             |
| RESULT_W_SEGMENT   | Segment algebraic variable - stored in result.w_segment                 |
"""
@enum ResultKind RESULT_ELIMINATED RESULT_CONSTANT RESULT_T RESULT_X RESULT_DER_X RESULT_W_INVARIANT RESULT_W_SEGMENT


"""
    locationID = LocationID(index, size)            # Location id of invariant variable
    locationID = LocationID(segment, index, size)   # Location id of segment variable

Return a new location id (defining the location of the values of a variable with respect to ResultKind).
"""
struct LocationID
    segment::Int  # Index of simulation segment (= -1, if invariant variable)
    index::Int    # Index or start index with respect to simulation segment (= -1, if time (RESULT_T))
    size::Tuple   # Dimensions with respect to simulation segment. size=(): scalar variable or size now known

    LocationID(index, size) = new(-1, index, size)
    LocationID(segment, index, size) = new(segment, index, size)
end


"""
    info = ResultInfo(kind, aliasName; negate=false)                    # Define alias ResultInfo
    info = ResultInfo(kind, value, unit)                                # Define constant ResultInfo
    info = ResultInfo(kind, default, unit, invariant)                   # Define partial ResultInfo of x or der_x
    info = ResultInfo(kind, default, unit, signalKind, index)           # Define ResultInfo for invariant variable
    info = ResultInfo(kind, default, unit, signalKind, segment, index)  # Define ResultInfo for segment variable

Return info how to access a result variable.
"""
struct ResultInfo
    kind::ResultKind                   # Kind of result variable in simulation segment sk at time instant ti and (invariant or varying) index:
                                       # = RESULT_ELIMINATED : Variable is eliminated. Alias info is stored in result.info
                                       # = RESULT_CONSTANT   : Variable is constant all the time. Value is stored in result.info
                                       # = RESULT_T          : result.t[          sk][ti]
                                       # = RESULT_X          : result.x[          sk][ti][index:index+prod(dims)-1]
                                       # = RESULT_DER_X      : result.der_x[      sk][ti][index:index+prod(dims)-1]
                                       # = RESULT_W_INVARIANT: result.w_invariant[sk][ti][index]
                                       # = RESULT_W_SEGMENT  : result.w_segment[  sk][ti][index]

    # Only if kind = RESULT_ELIMINATED
    aliasName::String                  # Name of non-eliminated variable
    aliasNegate::Bool                  # = true, if info[aliasName] signal must be negated

    # Only if kind = RESULT_CONSTANT
    value::Union{Any,Nothing}          # Value of constant variable (without unit)

    # For all kinds with exception of RESULT_ELIMINATED and RESULT_W_INVARIANT
    type::Union{DataType,Nothing}      # Type of variable (to make sure that the type is not changing)
    unit::String                       # Unit of variable as a parseable string.
    ndims::Int                         # Number of dimensions (= length(dims); to make sure that ndims is not changing)
                                       # ndims = -1, if type=nothing

    # For all kinds with exception of RESULT_ELIMINATED and RESULT_CONSTANT
    signalKind::ModiaResult.SignalType # Kind of signal
    invariant::Bool                    # = true, invariant variable
                                       # = false, segment variable

    # For all kinds with exception of RESULT_ELIMINATED, RESULT_CONSTANT and RESULT_T
    locationID::Vector{LocationID}    # if invariant = true : locationID[1] is location for all segments
                                      # if invariant = false: locationID[i] is segment specific location

    ResultInfo(kind::ResultKind, aliasName::String; negate::Bool=false)             = begin
                                                                                        @assert(kind == RESULT_ELIMINATED)
                                                                                        (kind, aliasName, negate)
                                                                                      end
    ResultInfo(kind::ResultKind, value, unit::String)                               = begin
                                                                                        @assert(kind == RESULT_CONSTANT)
                                                                                        (kind, "", false, value, typeof(value), unit, ndims(value))
                                                                                      end
    ResultInfo(kind::ResultKind, default, unit::String, invariant::Bool)            = begin
                                                                                        @assert(kind == RESULT_X || kind == RESULT_DER_X )
                                                                                        new(kind, "", false, nothing, typeof(default), unit, ndims(default), ModiaResult.Continuous, invariant, LocationID[])
                                                                                      end
    ResultInfo(kind::ResultKind, default, unit::String, signalKind, index)          = begin
                                                                                        @assert(kind != RESULT_ELIMINATED && kind != RESULT_CONSTANT)
                                                                                        new(kind, "", false, nothing, typeof(default), unit, isnothing(default) ? -1 : ndims(default),
                                                                                            signalKind, true, LocationID[LocationID(index  , isnothing(default) ? () : size(default))])
                                                                                      end
    ResultInfo(kind::ResultKind, default, unit::String, signalKind, segment, index) = begin
                                                                                        @assert(kind != RESULT_ELIMINATED && kind != RESULT_CONSTANT)
                                                                                        new(kind, "", false, nothing, typeof(default), unit, ndims(default), signalKind, false, LocationID[LocationID(segment,index,size(default))])
                                                                                      end
end


"""
    result = Result{FloatType,TimeType}(timeNameAsString, equationInfo, w_invariant_names, vEliminated, vProperty, var_name)

Return a new result data structure filled with invariant variable definitions.
"""
mutable struct Result{FloatType,TimeType}
    # Result access information
    info::OrderedDict{String,ResultInfo}        # Dictionary with key = variable name and value = info to access variable value
    timeName::String                            # Name of independent variable (usually time)
    timeResultInfo::ResultInfo                  # Reference to result.info[timeName]
    n_w_invariant::Int                          # Number of w_invariant variables
    alias_segment_names::OrderedSet{String}     # Names of alias_segment results in current segment
    w_segment_names::OrderedSet{String}         # Names of w_segment results in current segment
    w_segment_temp::Vector{Any}                 # Memory to store temporarily references to w_segment results in current segment;
                                                # length(w_segment_temp) = length(w_segment_names)
    t::Vector{Vector{TimeType}}                 # t[sk][ti] is time instant ti in segment sk

    # A variable v[sk][ti][j] is variable with index j at time instant ti in segment sk
    x::Vector{Vector{Vector{FloatType}}}        # x[sk][ti][j] - invariant and segment states
    der_x::Vector{Vector{Vector{FloatType}}}    # der_x[sk][ti][j] - invariant and segment state derivatives
    w_invariant::Vector{Vector{Tuple}}          # w_invariant[sk][ti][j] - invariant algebraic variables
    w_segment::Vector{Vector{Vector{Any}}}      # w_segment[sk][ti][j] - segment algebraic variables

    function Result{FloatType,TimeType}(timeNameAsString::String, eqInfo::EquationInfo, w_invariant_names, vEliminated, vProperty, var_name) where {FloatType,TimeType}
        info                = OrderedDict{String, ResultInfo}()
        n_w_invariant       = length(w_invariant_names)
        alias_segment_names = OrderedSet{String}()
        w_segment_names     = OrderedSet{String}()
        w_segment_temp      = Any[]
        t                   = fill(TimeType[],0)
        x                   = fill(Vector{FloatType}[], 0)
        der_x               = fill(Vector{FloatType}[], 0)
        w_invariant         = fill(Tuple[],0)
        w_segment           = fill(Vector{Any}[], 0)

        # Fill info with time
        timeResultInfo = ResultInfo(RESULT_T, TimeType(0), "s", ModiaResult.Independent, -1)
        info[timeNameAsString] = timeResultInfo

        # Fill info with x_invariant, der_x_invariant (but with dummy locationID, since not yet known)
        for i = 1:eqInfo.nx_info_invariant
            xi_info = eqInfo.x_info[i]
            x_unit     = xi_info.unit
            der_x_unit = x_unit == "" ? "1/s" : unitAsString(unit(uparse(x_unit)/u"s"))
            @assert(!haskey(info, xi_info.x_name))
            @assert(!haskey(info, xi_info.der_x_name))
            if isnothing(xi_info.startOrInit)
                xi_info.startOrInit = FloatType(0)
            end
            info[xi_info.x_name]     = ResultInfo(RESULT_X    , xi_info.startOrInit, x_unit    , true)
            info[xi_info.der_x_name] = ResultInfo(RESULT_DER_X, xi_info.startOrInit, der_x_unit, true)
        end

        # Fill info with w_invariant
        for (w_invariant_index, w_invariant_name) in enumerate(w_invariant_names)
            name = string(w_invariant_name)
            @assert(!haskey(info, name))
            info[name] = ResultInfo(RESULT_W_INVARIANT, nothing, "", ModiaResult.Continuous, w_invariant_index)
        end

        # Fill info with eliminated variables
        for v in vEliminated
            name = var_name(v)
            @assert(!haskey(info, name))
            if ModiaBase.isZero(vProperty, v)
                info[name] = ResultInfo(RESULT_CONSTANT, FloatType(0), "")
            elseif ModiaBase.isAlias(vProperty, v)
                aliasName = var_name( ModiaBase.alias(vProperty, v) )
                @assert(haskey(info, aliasName))
                info[name] = ResultInfo(RESULT_ELIMINATED, aliasName)
            else # negated alias
                negatedAliasName = var_name( ModiaBase.negAlias(vProperty, v) )
                @assert(haskey(info, negatedAliasName))
                info[name] = ResultInfo(RESULT_ELIMINATED, negatedAliasName, negate=true)
            end
        end

        new(info, timeNameAsString, timeResultInfo, n_w_invariant, alias_segment_names, w_segment_names, w_segment_temp, t, x, der_x, w_invariant, w_segment)
    end
end


function newResultSegment!(result::Result{FloatType,TimeType})::Nothing where {FloatType,TimeType}
    push!(result.t          , TimeType[])
    push!(result.x          , Vector{FloatType}[])
    push!(result.der_x      , Vector{FloatType}[])
    push!(result.w_invariant, Tuple[])
    push!(result.w_segment  , Vector{Any}[])
    return nothing
end


"""
    iterator = SegmentAndIndex(resultInfo::ResultInfo, nSegments::Int)

Iterate over the segments of result/resultInfo.
"""
struct SegmentAndIndex
    resInfo::ResultInfo
    nSegments::Int       # Number of segments
end

@inline function getInvariantLocation(resultInfo::ResultInfo, iter::Int)::NTuple{3,Int}
    id = resultInfo.locationID[1]
    return (iter, id.index, id.index+prod(id.size)-1)
end

@inline function getSegmentLocation(resultInfo::ResultInfo, iter::Int)::NTuple{3,Int}
    id = resultInfo.locationID[iter]
    return (id.segment, id.index, id.index+prod(id.size)-1)
end

Base.iterate(s::SegmentAndIndex, iter=1) = s.resInfo.invariant ? (iter > s.nSegments                  ? nothing : (getInvariantLocation(s.resInfo,iter), iter+1)) :
                                                                 (iter > length(s.resInfo.locationID) ? nothing : (getSegmentLocation(  s.resInfo,iter), iter+1))
Base.length(s::SegmentAndIndex) = s.resInfo.invariant ? s.nSegments : length(s.resInfo.locationID)


struct Segment
    resInfo::ResultInfo
    nSegments::Int       # Number of segments
end

Base.iterate(s::Segment, iter=1) = s.resInfo.invariant ? (iter > s.nSegments                  ? nothing : (iter, iter+1)) :
                                                         (iter > length(s.resInfo.locationID) ? nothing : (s.resInfo.locationID[iter].segment, iter+1))
Base.length(s::Segment) = s.resInfo.invariant ? s.nSegments : length(s.resInfo.locationID)


function rawSignal(result::Result, name::AbstractString; unitless=true)::Tuple{AbstractVector, AbstractVector, ModiaResult.SignalType}
    nSegments = length(result.t)
    resInfo   = result.info[name]
    negate    = false

    if resInfo.kind == RESULT_ELIMINATED
        resInfo = result.info[resInfo.aliasName]
        negate  = result.info[resInfo.aliasNegate]
    end

    tSig = [ result.t[sk] for sk in Segment(resInfo,nSegments) ]

    if resInfo.kind == RESULT_T
        sig = result.t
    elseif resInfo.kind == RESULT_X
        sig = [ [ibeg==iend ? x_ti[ibeg] : x_ti[ibeg:iend] for x_ti in result.x[sk]] for (sk,ibeg,iend) in SegmentAndIndex(resInfo,nSegments) ]
    elseif resInfo.kind == RESULT_DER_X
        sig = [ [ibeg==iend ? der_x_ti[ibeg] : der_x_ti[ibeg:iend] for der_x_ti in result.der_x[sk]] for (sk,ibeg,iend) in SegmentAndIndex(resInfo,nSegments) ]
    elseif resInfo.kind == RESULT_W_INVARIANT
        sig = [ [w_ti[index] for w_ti in result.w_invariant[sk]] for (sk,index,dummy) in SegmentAndIndex(resInfo,nSegments) ]
    elseif resInfo.kind == RESULT_W_SEGMENT
        sig = [ [w_ti[index] for w_ti in result.w_segment[sk]]   for (sk,index,dummy) in SegmentAndIndex(resInfo,nSegments) ]
    elseif resInfo.kind == RESULT_CONSTANT
        sig = [ ModiaResult.OneValueVector(resInfo.value, length(result.t[sk])) for sk = 1:length(result.t) ]
    else
        error("Bug in Modia.rawSignal: name=\"$name\" has ResultInfo=$resInfo, but ResultInof.kind = $(resInfo.kind) is not known.")
    end

    if negate
        sig *= -1
    end
    if resInfo.kind == RESULT_W_INVARIANT
        # Result has already unit, if compiled with unitless=false
        if unitless && !m.unitless
            sig = stripUnit(sig)
        end
    elseif !unitless && resInfo.unit != ""
        sig *= uparse(resInfo.unit)
    end
    return (tSig, sig, resInfo.signalKind)
end



# License for this file: MIT (expat)
# Copyright 2022, DLR Institute of System Dynamics and Control
#
# Modia result datastrucure and functions operating on results.

using  OrderedCollections: OrderedDict, OrderedSet

"""
    @enum ResultKind RESULT_ELIMINATED RESULT_CONSTANT RESULT_T RESULT_X RESULT_DER_X RESULT_W_INVARIANT RESULT_W_SEGMENTED

Kind of result variable.

| ResultKind value   | Description                                                              |
|:-------------------|:-------------------------------------------------------------------------|
| RESULT_ELIMINATED  | Variable is eliminated - alias info stored in result.info                |
| RESULT_CONSTANT    | Variable is constant all the time - value stored in result.info          |
| RESULT_T           | Independent variable (usually time) - stored in result.t                 |
| RESULT_X           | State (invariant/segmented variable) - stored in result.x                |
| RESULT_DER_X       | State derivative (invariant/segmented variable) - stored in result.der_x |
| RESULT_W_INVARIANT | Invariant algebraic variable - stored in result.w_invariant              |
| RESULT_W_SEGMENTED | Segmented algebraic variable - stored in result.w_segmented              |
"""
@enum ResultKind RESULT_ELIMINATED RESULT_CONSTANT RESULT_T RESULT_X RESULT_DER_X RESULT_W_INVARIANT RESULT_W_SEGMENTED


"""
    id = ValuesID(index, dims)            # ValuesID of invariant variable
    id = ValuesID(segment, index, dims)   # ValuesID of segmented variable

Return a new id that defines where the values of a variable are stored.
"""
struct ValuesID
    segment::Int               # Index of simulation segment (= -1, if index holds for every segment)
    index::Int                 # Index or start index with respect to simulation segment
    dims::Union{Dims,Nothing}  # Dimensions with respect to simulation segment; = nothing, if not yet defined (for RESULt_W_INVARIANT, before value is inquired)

    ValuesID(index, dims)          = new(     -1, index, dims)
    ValuesID(segment, index, dims) = new(segment, index, dims)
end


hasInlineValues(kind::ResultKind) = kind == RESULT_X || RESULT_DER_X
hasDims(        kind::ResultKind) = kind != RESULT_W_INVARIANT

isInvariant(id::Vector{ValuesID})                    = length(id) == 1
isSegmented(id::Vector{ValuesID}, t::AbstractVector) = !(length(id) == 1 || length(id) == length(t))

index_i(id::Vector{ValuesID}, i::Int)                                      = isInvariant(id) ? id[1].index : id[i].index
dims_i( id::Vector{ValuesID}, i::Int, kind::ResultKind, v::AbstractVector) = hasDims(kind) ? (isInvariant(id) ? id[1].dims : id[i].dims) :
                                                                                             (isInvariant(id) ? size(v[1][1][id[1].index]) : size(v[i][1][id[i].index]))


"""
    info = ResultInfo(kind, signal)                    # x, der_x (id is initially not known)
    info = ResultInfo(kind, signal, id)                # t, w_invariant, w_segmented
    info = ResultInfo(signal, value)                   # constant
    info = ResultInfo(signal, aliasName, aliasNegate)  # alias and negative alias

Return info how to access a result variable.
"""
mutable struct ResultInfo
    kind::ResultKind                 # Kind of result variable in simulation segment sk at time instant ti with index = index_i(..):
                                     # = RESULT_ELIMINATED : Variable is eliminated. Alias info is stored in result.info
                                     # = RESULT_CONSTANT   : Variable is constant all the time. Value is stored in result.info
                                     # = RESULT_T          : result.t[          sk][ti]
                                     # = RESULT_X          : result.x[          sk][ti][index:index+prod(dims_i(..))-1]
                                     # = RESULT_DER_X      : result.der_x[      sk][ti][index:index+prod(dims_i(..))-1]
                                     # = RESULT_W_INVARIANT: result.w_invariant[sk][ti][index]
                                     # = RESULT_W_SEGMENTED: result.w_segmented[sk][ti][index]
    signal::SignalTables.SymbolDictType # = SignalTables.Var() or SignalTables.Par()

    aliasName::String                # Name of non-eliminated variable
    aliasNegate::Bool                # = true, if info[aliasName] signal must be negated
    id::Vector{ValuesID}             # Location of the variable values with respect to ResultKind and Result
    _eltypeOrType                        # If known, eltypeOrType(signal.values/.values); if not known: Nothing
    value::Any                       # Value of constant variable (without unit)

    ResultInfo(kind::ResultKind, signal, _eltypeOrType)                                       = new(kind             , signal, ""       , false      , ValuesID[]  , _eltypeOrType, nothing)
    ResultInfo(kind::ResultKind, signal, id::ValuesID, _eltypeOrType)                         = new(kind             , signal, ""       , false      , ValuesID[id], _eltypeOrType, nothing)
    ResultInfo(signal::SignalTables.SymbolDictType, value)                                = new(RESULT_CONSTANT  , signal, ""       , false      , ValuesID[]  , eltypeOrType(value), value)
    ResultInfo(signal::SignalTables.SymbolDictType, aliasName::String, aliasNegate::Bool) = new(RESULT_ELIMINATED, signal, aliasName, aliasNegate, ValuesID[], Nothing, nothing)
end


"""
    result = Result{FloatType,TimeType}(timeNameAsString, equationInfo, w_invariant_names, w_invariant_initial, vEliminated, vProperty, var_name)

Return a new result data structure.
"""
mutable struct Result{FloatType,TimeType}
    # Result access information
    info::OrderedDict{String,ResultInfo}        # Dictionary with key = variable name and value = info to access variable value
    timeName::String                            # Name of independent variable (usually time)
    timeResultInfo::ResultInfo                  # Reference to result.info[timeName]
    n_w_invariant::Int                          # Number of w_invariant variables
    alias_segmented_names::OrderedSet{String}   # Names of alias_segmented results in current segment
    w_segmented_names::OrderedSet{String}       # Names of w_segmented results in current segment
    w_segmented_temp::Vector{Any}               # Memory to store temporarily references to w_segmented results in current segment
                                                # length(w_segmented_temp) = length(w_segmented_names)
    firstIndexOfSegment::Vector{Int}            # firstIndexOfSegment[sk] is the first index of the "expanded" time vector of segment sk
    t::Vector{Vector{TimeType}}                 # t[sk][ti] is time instant ti in segment sk

    # A variable v[sk][ti][j] is variable with index j at time instant ti in segment sk
    x::Vector{Vector{Vector{FloatType}}}        # x[sk][ti][j] - invariant and segmented states
    der_x::Vector{Vector{Vector{FloatType}}}    # der_x[sk][ti][j] - invariant and segmented state derivatives
    w_invariant::Vector{Vector{Tuple}}          # w_invariant[sk][ti][j] - invariant algebraic variables
    w_segmented::Vector{Vector{Vector{Any}}}    # w_segmented[sk][ti][j] - segmented algebraic variables

    function Result{FloatType,TimeType}(eqInfo::EquationInfo, timeNameAsString::String, w_invariant_names, vEliminated, vProperty, var_name) where {FloatType,TimeType}
        info                  = OrderedDict{String, ResultInfo}()
        n_w_invariant         = length(w_invariant_names)
        alias_segmented_names = OrderedSet{String}()
        w_segmented_names     = OrderedSet{String}()
        w_segmented_temp      = Any[]
        t                     = fill(TimeType[],1)
        x                     = fill(Vector{FloatType}[], 1)
        der_x                 = fill(Vector{FloatType}[], 1)
        w_invariant           = fill(Tuple[], 1)
        w_segmented           = fill(Vector{Any}[], 1)

        # Fill info with time
        firstIndexOfSegment = Int[1]
        timeResultInfo = ResultInfo(RESULT_T, SignalTables.Var(unit="s", independent=true), ValuesID(1,()), TimeType)
        info[timeNameAsString] = timeResultInfo

        # Fill info with x, der_x (note: id is not yet known, because init/start value might be changed in evaluatedParameters(..), which is called after Result(...)
        for i in 1:length(eqInfo.x_info)
            xi_info = eqInfo.x_info[i]
            @assert(!haskey(info, xi_info.x_name))
            @assert(!haskey(info, xi_info.der_x_name))
            x_unit     = xi_info.unit
            der_x_unit = x_unit == "" ? "1/s" : unitAsParseableString(uparse(x_unit)/u"s")
            if x_unit == ""
                x_var = SignalTables.Var(start=xi_info.startOrInit, fixed=xi_info.fixed, state=true, der=xi_info.der_x_name)
            else
                x_var = SignalTables.Var(unit=x_unit, start=xi_info.startOrInit, fixed=xi_info.fixed, state=true, der=xi_info.der_x_name)
            end
            if !isnan(xi_info.nominal)
                x_var[:nominal] = xi_info.nominal
            end
            if xi_info.unbounded
                x_var[:unbounded] = true
            end
            info[xi_info.x_name]     = ResultInfo(RESULT_X, x_var, FloatType)
            if der_x_unit == ""
                info[xi_info.der_x_name] = ResultInfo(RESULT_DER_X, SignalTables.Var(), FloatType)
            else
                info[xi_info.der_x_name] = ResultInfo(RESULT_DER_X, SignalTables.Var(unit=der_x_unit), FloatType)
            end
        end

        # Fill info with w_invariant
        for (w_invariant_index, w_invariant_name) in enumerate(w_invariant_names)
            name = string(w_invariant_name)
            @assert(!haskey(info, name))
            info[name] = ResultInfo(RESULT_W_INVARIANT, SignalTables.Var(), ValuesID(w_invariant_index, nothing), Nothing)
        end

        # Fill info with eliminated variables
        for v in vEliminated
            name = var_name(v)
            @assert(!haskey(info, name))
            if ModiaBase.isZero(vProperty, v)
                info[name] = ResultInfo(SignalTables.Var(), FloatType(0))
            elseif ModiaBase.isAlias(vProperty, v)
                aliasName = var_name( ModiaBase.alias(vProperty, v) )
                @assert(haskey(info, aliasName))
                info[name] = ResultInfo(SignalTables.Var(), aliasName, false)
            else # negated alias
                negatedAliasName = var_name( ModiaBase.negAlias(vProperty, v) )
                @assert(haskey(info, negatedAliasName))
                info[name] = ResultInfo(SignalTables.Var(), negatedAliasName, true)
            end
        end

        new(info, timeNameAsString, timeResultInfo, n_w_invariant, alias_segmented_names, w_segmented_names, w_segmented_temp,
            firstIndexOfSegment, t, x, der_x, w_invariant, w_segmented)
    end
end


"""
    nResults(result::Result)

Returns the number of result points (= number of values of time vector).
"""
nResults(result::Result) = result.firstIndexOfSegment[end] + length(result.t[end]) - 1


"""
    nResultsCurrentSegment(result::Result)

Returns the number of result points of the current segment (= number of values of time vector).
"""
nResultsCurrentSegment(result::Result) = length(result.t[end])


"""
    newResultSegment!(result, equationInfo, nsegments)

Start a new result segment (nsegments > 1)
"""
function newResultSegment!(result::Result{FloatType,TimeType}, equationInfo::EquationInfo, nsegments::Int)::Nothing where {FloatType,TimeType}
    @assert(nsegments > 1)
    empty!(result.alias_segmented_names)
    empty!(result.w_segmented_names)
    empty!(result.w_segmented_temp)

    # Start new segment
    push!(result.firstIndexOfSegment, result.firstIndexOfSegment[end] + length(result.t[end]))
    push!(result.t          , TimeType[])
    push!(result.x          , Vector{FloatType}[])
    push!(result.der_x      , Vector{FloatType}[])
    push!(result.w_invariant, Tuple[])
    push!(result.w_segmented, Vector{Any}[])

    return nothing
end


dims_range(dims::Dims) = Tuple([1:i for i in dims])


"""
    signalResultValues(t, s, resultInfo::ResultInfo; log=false, name="")

Return a SignalTables.Var() values vector from independent values t, dependent values s, and resultInfo.
"""
function signalResultValues(t::AbstractVector, s::AbstractVector, resultInfo::ResultInfo, result::Result; log=false, name::AbstractString="")
    id = resultInfo.id
    @assert(length(id) > 0)
    inlineValues = resultInfo.kind == RESULT_X || resultInfo.kind == RESULT_DER_X
    _eltypeOrType    = resultInfo._eltypeOrType
    ndims_s      = length(id[1].dims)
    if id[1].segment == -1 && ndims_s == 0
        # Scalar signal that is defined in every segment
        index = id[1].index
        sc = _eltypeOrType[ustrip.(ti[index]) for sk in s for ti in sk]

    else
        # Find largest dims = dimsMax in all segments
        dimsMax::Dims{ndims_s} = id[1].dims
        invariant::Bool  = id[1].segment==-1   # = true, if variable is defined at all time instants and variable index does not change
        hasMissing::Bool = !(invariant || length(id) == length(t))   # = true, if at least one variable value is not defined at all time instants
        if length(id) > 1
            dMax::Vector{Int} = [dimsMax...]
            for i = 2:length(id)
                dimi::Dims{ndims_s} = id[i].dims
                @assert(length(dMax) == length(dimi))
                if dimi != dMax
                    hasMissing = true
                    for j = 1:length(dimi)
                        if dimi[j] > dMax[j]
                            dMax[j] = dimi[j]
                        end
                    end
                end
            end
            dimsMax = Tuple(dMax)
        end

        # Allocate memory for target signal
        dims1 = sum(size(tk,1) for tk in t)
        dims  = (dims1, dimsMax...)
        if hasMissing
            # Allocate target memory with missing values
            sc = Array{Union{_eltypeOrType,Missing}, length(dims)}(missing, dims)
        else
            # Allocate target memory with undef values
            sc = Array{_eltypeOrType, length(dims)}(undef, dims)
        end

        # Copy subset of s-values to target sc
        j = 1
        firstIndexOfSegment = result.firstIndexOfSegment
        if length(dimsMax) == 0
            # Target is a scalar signal
            @assert(id[1].segment != -1)
            for id_k in resultInfo.id
                index   = id_k.index
                segment = id_k.segment
                j       = firstIndexOfSegment[segment]
                for s_ti in s[segment]
                    sc[j] = s_ti[index]
                    j += 1
                end
            end
        else
            # Target is not a scalar signal  setindex!(A,x,(2,2:4)...)
            j = 1
            if inlineValues
                if invariant
                    dims = id[1].dims
                    dimr = dims_range(dims)
                    ibeg = id[1].index
                    iend = ibeg + prod(dims) - 1
                    for sk in s
                        for s_ti in sk
                            setindex!(sc, reshape(view(s_ti,ibeg:iend),dims), (j,dimr...)...)
                            j += 1
                        end
                    end
                else
                    for id_k in resultInfo.id
                        dims    = id_k.dims
                        dimr    = dims_range(dims)
                        ibeg    = id_k.index
                        iend    = ibeg + prod(dims) - 1
                        segment = id_k.segment
                        j       = firstIndexOfSegment[segment]
                        for s_ti in s[segment]
                            setindex!(sc, reshape(view(s_ti,ibeg:iend),dims), (j,dimr...)...)
                            j += 1
                        end
                    end
                end
            else
                if invariant
                    index = id[1].index
                    dimr  = dims_range(id[1].dims)
                    for sk in s
                        for s_ti in sk
                            setindex!(sc, ustrip.(s_ti[index]), (j,dimr...)...)
                            j += 1
                        end
                    end
                else
                    for id_k in resultInfo.id
                        index   = id_k.index
                        dimr    = dims_range(id_k.dims)
                        segment = id_k.segment
                        j       = firstIndexOfSegment[segment]
                        for s_ti in s[segment]
                            setindex!(sc, ustrip.(s_ti[index]), (j,dimr...)...)
                            j += 1
                        end
                    end
                end
            end
        end
    end

    if log
        println("$name[id] = $sc")
        println("typeof($name[id]) = ", typeof(sc))
    end
    return sc
end
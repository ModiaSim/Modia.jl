# License for this file: MIT (expat)
# Copyright 2022, DLR Institute of System Dynamics and Control
#
# Modia result datastrucure and functions operating on results.

using  OrderedCollections: OrderedDict, OrderedSet
import ModiaResult

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
    segment::Int               # Index of simulation segment (= -1, if invariant variable)
    index::Int                 # Index or start index with respect to simulation segment
    dims::Union{Dims,Missing}  # Dimensions with respect to simulation segment (dims=missing, if dimensions not known)

    ValuesID(index, dims) = new(-1, index, dims)
    ValuesID(segmented, index, dims) = new(segmented, index, dims)
end


hasIinlineValues(kind::ResultKind) = kind == RESULT_X || RESULT_DER_X
hasDims(          kind::ResultKind) = kind != RESULT_W_INVARIANT

isInvariant(id::Vector{ValuesID})                    = length(id) == 1
isSegmented(id::Vector{ValuesID}, t::AbstractVector) = !(length(id) == 1 || length(id) == length(t))

index_i(id::Vector{ValuesID}, i::Int)                                      = isInvariant(id) ? id[1].index : id[i].index 
dims_i( id::Vector{ValuesID}, i::Int, kind::ResultKind, v::AbstractVector) = hasDims(kind) ? (isInvariant(id) ? id[1].dims : id[i].dims) :
                                                                                             (isInvariant(id) ? size(v[1][1][id[1].index]) : size(v[i][1][id[i].index]))                      
                                                                                                                                  

"""
    info = ResultInfo(kind, aliasName; negate=false)               # Define alias ResultInfo
    info = ResultInfo(kind, defaultOrValue, unit)                  # Define constant ResultInfo or partial ResultInfo of x or der_x
    info = ResultInfo(kind, default       , unit, index)           # Define ResultInfo for invariant variable
    info = ResultInfo(kind, default       , unit, segment, index)  # Define ResultInfo for segmented variable

Return info how to access a result variable.
"""
struct ResultInfo
    kind::ResultKind  # Kind of result variable in simulation segment sk at time instant ti with index = index_i(..):
                      # = RESULT_ELIMINATED : Variable is eliminated. Alias info is stored in result.info
                      # = RESULT_CONSTANT   : Variable is constant all the time. Value is stored in result.info
                      # = RESULT_T          : result.t[          sk][ti]
                      # = RESULT_X          : result.x[          sk][ti][index:index+prod(dims_i(..))-1]
                      # = RESULT_DER_X      : result.der_x[      sk][ti][index:index+prod(dims_i(..))-1]
                      # = RESULT_W_INVARIANT: result.w_invariant[sk][ti][index]
                      # = RESULT_W_SEGMENTED: result.w_segmented[sk][ti][index]

    aliasName::String                     # Name of non-eliminated variable
    aliasNegate::Bool                     # = true, if info[aliasName] signal must be negated
    value::Union{Any,Missing}             # Value of constant variable (without unit)
    
    VariableType::Union{DataType,Missing} # Type of variable, if known (to make sure that the VariableType is not changing)
    unit::String                          # Unit of variable as a parseable string. If not known, unit="".
    id::Vector{ValuesID}                  # Location of the variable values with respect to ResultKind and Result

    ResultInfo(kind::ResultKind, aliasName::String; negate::Bool=false) = begin
                                                                            @assert(kind == RESULT_ELIMINATED)
                                                                            (kind, aliasName, negate, missing, missing, "", ValuesID[])
                                                                          end
    ResultInfo(kind::ResultKind, defaultOrValue, unit::String)          = begin
                                                                            @assert(kind == RESULT_CONSTANT || kind == RESULT_X || kind == RESULT_DER_X)
                                                                            if kind == RESULT_CONSTANT
                                                                                (kind, "", false, defaultOrValue, typeof(defaultOrValue), unit, ValuesID[])
                                                                            else
                                                                                (kind, "", false, missing, typeof(defaultOrValue), unit, ValuesID[])
                                                                            end
                                                                          end
    ResultInfo(kind::ResultKind, default, unit::String, index)          = ResultInfo(kind,default,unit,-1,index)
    ResultInfo(kind::ResultKind, default, unit::String, segment, index) = begin
                                                                            @assert(kind != RESULT_ELIMINATED && kind != RESULT_CONSTANT)
                                                                            new(kind, "", false, default, unit, ValuesID[ValuesID(segment,index,ismissing(default) ? missing : size(default))])
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
    alias_segmented_names::OrderedSet{String}   # Names of alias_segmented results in current segment
    w_segmented_names::OrderedSet{String}       # Names of w_segmented results in current segment
    w_segmented_temp::Vector{Any}               # Memory to store temporarily references to w_segmented results in current segment
                                                # length(w_segmented_temp) = length(w_segmented_names)
    t::Vector{Vector{TimeType}}                 # t[sk][ti] is time instant ti in segment sk

    # A variable v[sk][ti][j] is variable with index j at time instant ti in segmented sk
    x::Vector{Vector{Vector{FloatType}}}        # x[sk][ti][j] - invariant and segmented states
    der_x::Vector{Vector{Vector{FloatType}}}    # der_x[sk][ti][j] - invariant and segmented state derivatives
    w_invariant::Vector{Vector{Tuple}}          # w_invariant[sk][ti][j] - invariant algebraic variables
    w_segmented::Vector{Vector{Vector{Any}}}    # w_segmented[sk][ti][j] - segmented algebraic variables

    function Result{FloatType,TimeType}(timeNameAsString::String, eqInfo::EquationInfo, w_invariant_names, vEliminated, vProperty, var_name) where {FloatType,TimeType}
        info                  = OrderedDict{String, ResultInfo}()
        n_w_invariant         = length(w_invariant_names)
        alias_segmented_names = OrderedSet{String}()
        w_segmented_names     = OrderedSet{String}()
        w_segmented_temp      = Any[]
        t                     = fill(TimeType[],0)
        x                     = fill(Vector{FloatType}[], 0)
        der_x                 = fill(Vector{FloatType}[], 0)
        w_invariant           = fill(Tuple[],0)
        w_segmented           = fill(Vector{Any}[], 0)

        # Fill info with time
        timeResultInfo = ResultInfo(RESULT_T, TimeType(0), "s", 1)
        info[timeNameAsString] = timeResultInfo

        # Fill info with x_invariant, der_x_invariant (but with dummy id, since not yet known)
        for i = 1:eqInfo.nx_info_invariant
            xi_info = eqInfo.x_info[i]
            x_unit  = xi_info.unit
            der_x_unit = x_unit == "" ? "1/s" : unitAsString(unit(uparse(x_unit)/u"s"))
            @assert(!haskey(info, xi_info.x_name))
            @assert(!haskey(info, xi_info.der_x_name))
            if isnothing(xi_info.startOrInit)
                xi_info.startOrInit = FloatType(0)
            end
            info[xi_info.x_name]     = ResultInfo(RESULT_X    , xi_info.startOrInit, x_unit    )
            info[xi_info.der_x_name] = ResultInfo(RESULT_DER_X, xi_info.startOrInit, der_x_unit)
        end

        # Fill info with w_invariant
        for (w_invariant_index, w_invariant_name) in enumerate(w_invariant_names)
            name = string(w_invariant_name)
            @assert(!haskey(info, name))
            info[name] = ResultInfo(RESULT_W_INVARIANT, missing, "", w_invariant_index)
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

        new(info, timeNameAsString, timeResultInfo, n_w_invariant, alias_segmented_names, w_segmented_names, w_segmented_temp, t, x, der_x, w_invariant, w_segmented)
    end
end


function newResultSegment!(result::Result{FloatType,TimeType})::Nothing where {FloatType,TimeType}
    push!(result.t          , TimeType[])
    push!(result.x          , Vector{FloatType}[])
    push!(result.der_x      , Vector{FloatType}[])
    push!(result.w_invariant, Tuple[])
    push!(result.w_segmented  , Vector{Any}[])
    return nothing
end


dims_range(dims::Dims) = Tuple([1:i for i in dims])
dims_i( inlineValues::Bool, id::Vector{ValuesID}, k::Int, s::AbstractVector) = inlineValues ? id[k].dims         : size( s[k][1][id[k].index])  
ndims_i(inlineValues::Bool, id::Vector{ValuesID}, k::Int, s::AbstractVector) = inlineValues ? length(id[k].dims) : ndims(s[k][1][id[k].index]) 
   
signalResultValues(t::AbstractVector, s::AbstractVector, resultInfo::ResultInfo; log=false, name::AbstractString="") =
    signalResultValues(t, s, resultInfo.id, kind == RESULT_X || kind == RESULT_DER_X, resultInfo.VariableType, log=log, name=name)

function signalResultValues(t::AbstractVector, s::AbstractVector, id::Vector{ValuesID}, inlineValues::Bool, VariableType=Float64; log=false, name::AbstractString="")  
    @assert(length(id) > 0)
    
    if length(id) == 1 && ndims_i(inlineValues,id,1,s) == 0
        # Scalar signal that is defined in every segment
        index = id[1].index
        sc = VariableType[ti[index] for sk in s for ti in sk]
        
    else
        # Find largest dims = dimsMax in all segments
        dimsMax::Dims{ndims_i(inlineValues,id,1,s)} = dims_i(inlineValues,id,1,s)
        hasMissing::Bool = !(length(id) == 1 || length(id) == length(t))
        if length(id) > 1
            dMax::Vector{Int} = [dimsMax...]
            for i = 2:length(id)
                dimi = dims_i(inlineValues,id,i,s)
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
            sc = Array{Union{VariableType,Missing}, length(dims)}(missing, dims)
        else
            # Allocate target memory with undef values
            @assert(length(dimsMax) > 0)
            sc = Array{VariableType, length(dims)}(undef, dims)
        end
        
        # Copy subset of s-values to target sc
        invariant  = length(id) == 1
        segmentOld = 1        
        j = 1
        if length(dimsMax) == 0
            # Target is a scalar signal
            if invariant
                index = id[1].index
            end
            for k = 1:length(s)
                if !invariant
                    index   = id[k].index
                    segment = id[k].segment  
                    if segmentOld < segment-1
                        j = sum(length(t[i]) for i = segmentOld:segment-1) + 1
                        segmentOld = segment
                    end
                end
                for s_ti in s[k]
                    sc[j] = s_ti[index]
                    j += 1
                end
            end            
        else
            # Target is not a scalar signal  setindex!(A,x,(2,2:4)...)
            if inlineValues
                if invariant
                    dims = id[1].dims
                    dimr = dims_range(dims)
                    ibeg = id[1].index
                    iend = ibeg + prod(dims) - 1
                end
                for k = 1:length(s)
                    if !invariant
                        dims = id[k].dims
                        dimr = dims_range(dims)
                        ibeg = id[k].index
                        iend = ibeg + prod(dims) - 1
                        segment = id[k].segment
                        if segmentOld < segment-1
                            j = sum(length(t[i]) for i = segmentOld:segment-1) + 1
                            segmentOld = segment
                        end
                    end                  
                    for s_ti in s[k]
                        setindex!(sc, reshape(view(s_ti,ibeg:iend),dims), (j,dimr...)...)
                        j += 1
                    end
                end
            else
                if invariant
                    index = id[1].index
                    dimr  = dims_range( size(s[1][1][index]) )
                end            
                for k = 1:length(s)
                    if !invariant
                        index   = id[k].index
                        dimr    = dims_range( size(s[k][1][index]) )
                        segment = id[k].segment                        
                        if segmentOld < segment-1
                            j = sum(length(t[i]) for i = segmentOld:segment-1) + 1
                            segmentOld = segment
                        end                        
                    end
                    for s_ti in s[k]
                        setindex!(sc, s_ti[index], (j,dimr...)...)
                        j += 1
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


function signalResultValues(result::Result, name::AbstractString; unitless=true)
    resInfo = result.info[name]
    negate  = false

    if resInfo.kind == RESULT_ELIMINATED
        resInfo = result.info[resInfo.aliasName]
        negate  = result.info[resInfo.aliasNegate]
    end

    if resInfo.kind == RESULT_T
        sig = signal(result.t, result.t, resInfo)
    elseif resInfo.kind == RESULT_X
        sig = signal(result.t, result.x, resInfo)
    elseif resInfo.kind == RESULT_DER_X
        sig = signal(result.t, result.der_x, resInfo)
    elseif resInfo.kind == RESULT_W_INVARIANT
        sig = signal(result.t, result.w_invariant, resInfo)
    elseif resInfo.kind == RESULT_W_SEGMENTED
        sig = signal(result.t, result.w_segmented, resInfo)
    elseif resInfo.kind == RESULT_CONSTANT
        sig = ModiaResult.OneValueVector(resInf.value, sum(length(tk) for tk in result.t))
    else
        error("Bug in Modia.signal: name=\"$name\" has ResultInfo=$resInfo, but ResultInfo.kind = $(resInfo.kind) is not known.")
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
    
    return sig
end
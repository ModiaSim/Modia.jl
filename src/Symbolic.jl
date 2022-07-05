# Moved from ModiaBase/src/Symbolic.jl, because dependent on FloatType that is only available in ModiaLang

prependPar(ex, prefix, parameters=[], inputs=[]) = ex

function prependPar(ex::Symbol, prefix, parameters=[], inputs=[]) 
    if prefix == nothing; ex elseif ex in [:time, :instantiatedModel, :_leq_mode, :_x]; ex else Expr(:ref, prefix, QuoteNode(ex)) end
end

function prependPar(ex::Expr, prefix, parameters=[], inputs=[])
    if isexpr(ex, :.)
        Expr(:ref, prependPar(ex.args[1], prefix, parameters, inputs), QuoteNode(ex.args[2].value))
    else
        Expr(ex.head, [prependPar(arg, prefix, parameters, inputs) for arg in ex.args]...)
    end
end



function addCastAndTypeCheckBasicValueType(ex, value, FloatType)
    # typeof(FloatType) is neither Measurements nor MonteCarloMeasurements
    T = baseType(FloatType)
    if isQuantity(typeof(value))
        ustr = unitAsParseableString(unit(value))
        :( Modia.quantity($T, @u_str($ustr))($ex)::Modia.quantity($T, @u_str($ustr)) )                 
    else
        :( $T($ex)::$T )
    end
end

                
"""
    e = addCastAndTypeCheck(ex,value,FloatType)

If `eltype(value) <: Number` add either 

1. RequiredType(ex)::RequiredType or
2. typeof(value)(ex)::typeof(value)

If it is manageable, variant (1) is used and the type cast is made with respect to _FloatType.
This is done for scalars that are no Measurements and no MonteCarloMeasurements.
Otherwise, variant (2) is used.

As a result, the generated code is type stable and does not work with Any types 
(although the values are stored in an Any dictionary).
So the code will be more efficient and memory allocated at run-time will be reduced.
Furthermore, if value is a scalar that has a unit, the conversion includes also the conversion
to the unit used during code generation.

Note, this function should only be used on parameter values.
"""
function addCastAndTypeCheck(ex,value,FloatType)
    if eltype(value) <: Number
        # Numeric scalar or array
        valueType = typeof(value)
        
        if valueType <: Number
            # Numeric scalar
            T = baseType(FloatType)            
            if FloatType <: Measurements.Measurement
                if isMeasurements(valueType)               
                    if isQuantity(valueType)
                        ustr = unitAsParseableString(unit(value))
                        :( Modia.quantity(Modia.Measurements.Measurement{$T}, @u_str($ustr))($ex)::Modia.quantity(Modia.Measurements.Measurement{$T}, @u_str($ustr)) )                         
                    else
                        :( Modia.Measurements.Measurement{$T}($ex)::Modia.Measurements.Measurement{$T} )
                    end
                else
                    addCastAndTypeCheckBasicValueType(ex,value,T)
                end
                
            elseif FloatType <: MonteCarloMeasurements.AbstractParticles
                if isMonteCarloMeasurements(valueType)
                    if isQuantity(valueType)
                        error("Error in Modia (MonteCarloMeasurements with Quantities appears - which is not supported.)")
                    else
                        :( ($ex)::$valueType )
                    end
                else
                    addCastAndTypeCheckBasicValueType(ex,value,T)
                end
                
            elseif valueType <: AbstractFloat || isQuantity(valueType)
                # FloatType is neither Measurements nor MonteCarloMeasurements and
                # is either an AbstractFloat or has a unit (is treated as Float)
                addCastAndTypeCheckBasicValueType(ex,value,T)
                
            else
                # For example Bool or Int values: No conversion, just define the type in the code
                :( ($ex)::$valueType )
            end            
        else
            # Numeric array
            :( ($ex)::$valueType )
        end
    else
        # No numeric scalar or array
        ex
    end
end

"""
    e = makeDerVar(ex, parameters, inputs, FloatType, evaluateParameters=false)

Recursively converts der(x) to Symbol(:(der(x))) in expression `ex`

* `ex`: Expression or array of expressions
* `return `e`: ex with der(x) converted 
"""
function makeDerVar(ex, parameters, inputs, FloatType, evaluateParameters=false) 
    if typeof(ex) in [Symbol, Expr]
        if ex in keys(parameters)
            addCastAndTypeCheck( prependPar(ex, :(_p), parameters, inputs),  parameters[ex], FloatType )
        elseif ex in keys(inputs)
            prependPar(ex, :(_p), parameters, inputs)
        else 
            ex
        end
    else
        ex
    end
end

function makeDerVar(ex::Expr, parameters, inputs, FloatType, evaluateParameters=false)
    if ex.head == :call && ex.args[1] == :der
        Symbol(ex)
	elseif isexpr(ex, :.) && ex in keys(parameters)
        if evaluateParameters
            parameters[ex]
        else
            addCastAndTypeCheck( prependPar(ex, :(_p), parameters, inputs), parameters[ex], FloatType )
        end
	elseif isexpr(ex, :.) && ex in keys(inputs)
        if evaluateParameters
            inputs[ex]
        else
            prependPar(ex, :(_p), parameters, inputs)
        end
    elseif ex.head == :.
        Symbol(ex)
    elseif ex.head == :call # Don't change dot-notation for function calls
        Expr(ex.head, ex.args[1], [makeDerVar(arg, parameters, inputs, FloatType, evaluateParameters) for arg in ex.args[2:end]]...)
    else
        Expr(ex.head, [makeDerVar(arg, parameters, inputs, FloatType, evaluateParameters) for arg in ex.args]...)
    end
end
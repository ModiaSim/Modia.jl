"""
Encoding and decoding Modia models as JSON.

* Developer: Hilding Elmqvist, Mogram AB
* First version: August 2021
* License: MIT (expat)

"""
module JSONModel

export modelToJSON, JSONToModel, cloneModel, writeModel, readModel

import JSON
using Modia
using Unitful
using OrderedCollections
using ModiaBase.Symbolic: removeBlock
using Measurements
import MonteCarloMeasurements

"""
    jModel = encodeModel(model, expressionsAsStrings=true)

Encodes a model suitable to convert to a JSON string.

* `model`: model (declarations and equations)
* `expressionsAsStrings`: Julia expressions stored as strings or ast (Expr struct, etc) stored as dictionary
* `return transformed model`
"""
function encodeModel(model, expressionsAsStrings=true)
    Base.remove_linenums!(model)
    model = removeBlock(model)
    if typeof(model) <: Quantity
        return Dict(:_value => ustrip.(model), :_unit=>SignalTables.unitAsParseableString(unit(model)))
    elseif expressionsAsStrings && typeof(model) == Expr
        return Dict(:_Expr=>string(model))
    elseif typeof(model) == Expr
        return Dict("head"=>model.head, "args"=>[encodeModel(a, expressionsAsStrings) for a in model.args])
    elseif typeof(model) == Symbol
        return Dict(:_Symbol=>string(model))
    elseif typeof(model) <: AbstractDict 
        jModel = OrderedDict()
        for (n,v) in model
            jModel[n] = encodeModel(v, expressionsAsStrings)
        end
        return jModel
    elseif typeof(model) <: QuoteNode || typeof(model) <: Measurement || typeof(model) <: MonteCarloMeasurements.StaticParticles || typeof(model) <: MonteCarloMeasurements.Particles
        return Dict(:_type=>"$(typeof(model))", :_value => string(model))
    elseif typeof(model) <: AbstractArray
        return [encodeModel(m, expressionsAsStrings) for m in model]
    else
        return model
    end
end

"""
    json = modelToJSON(model; expressionsAsStrings=true)

Encodes a model as JSON. Expressions can be stored as strings or as ASTs.

* `model`: model (declarations and equations)
* `expressionsAsStrings`: Julia expressions stored as strings or Expr struct stored as dictionary
* `return JSON string for the model`
"""
function modelToJSON(model; expressionsAsStrings=true)
    jModel = encodeModel(model, expressionsAsStrings)
    return JSON.json(jModel, 2)
end



"""
    model = decodeModel(jModel)

Encodes a model suitable to convert to a JSON string.

* `jmodel`: encoded model
* `return model`
"""
function decodeModel(jModel)
    if typeof(jModel) <: AbstractDict && length(keys(jModel)) == 2 && "_value" in keys(jModel) && "_unit" in keys(jModel)
        unit = jModel["_unit"]
        #unit = replace(unit, ' ' => '*')  # Workaround for https://github.com/PainterQubits/Unitful.jl/issues/391  # Remove, since already correctly done in encodeModel(...)
        return jModel["_value"]*uparse(unit)
    elseif typeof(jModel) <: AbstractDict && length(keys(jModel)) == 1 && "_Expr" in keys(jModel) 
        return Base.remove_linenums!( Meta.parse(jModel["_Expr"]) )  # without remove_linenums!, there are errors later with if-expressions.
    elseif typeof(jModel) <: AbstractDict && length(keys(jModel)) == 1 && "_Symbol" in keys(jModel) 
        return Symbol(jModel["_Symbol"])
    elseif typeof(jModel) <: AbstractDict && length(keys(jModel)) == 2 && "head" in keys(jModel) && "args" in keys(jModel)
        ast = Expr(Symbol(jModel["head"]), decodeModel.(jModel["args"])...)
        return ast
    elseif typeof(jModel) <: AbstractDict && length(keys(jModel)) == 2 && "_type" in keys(jModel) && "_value" in keys(jModel)
        if jModel["_type"] == "QuoteNode"
            return QuoteNode(eval(Meta.parse(jModel["_value"])))
        else
            return eval(Meta.parse(jModel["_value"]))
        end
    elseif typeof(jModel) <: AbstractDict 
        model = OrderedDict()
        for (n,v) in jModel
            model[Symbol(n)] = decodeModel(v)
        end
        return model
    elseif typeof(jModel) <: AbstractArray
        return [decodeModel(m) for m in jModel]
    else
        return jModel
    end
end

"""
    jModel = JSONToModel(json)

Decodes JSON string as model.

* `json`: JSON string for the model
* `return model`
"""
function JSONToModel(json)
    jModel = JSON.parse(json, dicttype=OrderedDict)
    model = decodeModel(jModel)
end


"""
    writeModel(filename::String, model; log=true)
    
Write model in JSON format on file `filename`.
"""
function writeModel(filename::String, model::AbstractDict; log=true) 
    if log
        println("  Write model in JSON format on file \"$filename\"")
    end
    write(filename, modelToJSON(model))
end


"""
    model = readModel(filename::AbstractString; log=true)
    
Read a model that is stored in JSON format on file `filename` and return it.
"""
function readModel(filename::String; log=true)
    if log
        println("  Read model from JSON file \"$filename\"")
    end
    JSONToModel( read(filename,String) )
end



"""
    modelClone = cloneModel(model; expressionsAsStrings=true)

Clones a model by encoding/decoding as JSON.

* `model`: model
* `expressionsAsStrings`: Julia expressions stored as strings or Expr struct stored as dictionary
* `return modelClone`
"""
function cloneModel(model; expressionsAsStrings=true)
#    @showModel(model)
    jsonModel = modelToJSON(model, expressionsAsStrings=expressionsAsStrings)
#    @show jsonModel
    modelClone = JSONToModel(jsonModel)
#    @showModel(modelClone)
    return modelClone
end

# --------------------------------------------------------------------

#=
function test()

    model = Model(
        p = 3u"m/s",
        q = 4.5,
        r = :(p+q),
        v = Var(min=0),
        s = "asdf qwerty",  
        t = 1 ± 0.1,
#        u = 1 ∓ 0.1,
        equations = :[
            2*der(x) = -p*x
            y = 2*x
        ]
    )

    model2 = Model(
#        t = 1 ± 0.1,
#        u = 1 ∓ 0.1,
    )

    println("Original model:")
    @showModel(model)
    println()
    @instantiateModel(model, log=true)

    json = modelToJSON(model)

    println("JSON for model:")
    println(json)

    println()

    modelClone = JSONToModel(json)
    @show modelClone
    println("Cloned model:")
    @showModel(modelClone)
    @instantiateModel(modelClone, log=true)

    println()

    json = modelToJSON(modelClone)
    println("JSON for cloned model:")
    println(json)

    # ------------------------

    Pin = Model( v = potential, i = flow )
    Ground = Model( p = Pin, equations = :[ p.v = 0.0u"V" ] ) 
    Gro = cloneModel(Ground)
    @showModel(Gro)
    @instantiateModel(Gro)
end

test()
=#

end

"""
Modia module with utility functions.

* Developer: Hilding Elmqvist, Mogram AB
* First version: July 2016
* Copyright (c) 2016-2018: Hilding Elmqvist, Toivo Henningsson, Martin Otter
* License: MIT (expat)

"""
module Utilities

using ..Instantiation: GetField, Der, InitVariable, Extends, Instance, Variable, Variability, constant, parameter, discrete, continuous, Connect, vars_of
using ..Execution: split_variables
using Base.Meta: quot, isexpr
using Unitful
using ..Instantiation

using ..ModiaLogging

export @show_io, printSymbolList, showModel, showInstance, checkSizes, printJSON

using Base: show_unquoted

# Version of @show for any stream
macro show_io(io, exs...)
    blk = Expr(:block)
    for ex in exs
        push!(blk.args, :(loglnModia($io, $(sprint(show_unquoted, ex) * " = "),
                                  repr(begin value = $(esc(ex)) end))))
    end
    if !isempty(exs); push!(blk.args, :value); end
    return blk
end

const maxSymbols = 1000

function printSymbolList(label, symbols, numbering=false, vertical=false, A=[])
    logModia(label, ": ")

    if vertical
        loglnModia()
    end

    for i in 1:min(maxSymbols, length(symbols))
        if i > 1
            if vertical
                loglnModia()
            else
                logModia(", ") 
            end
        end

        if numbering
            if vertical
                logModia(lpad(i, 5, " "), ": ")
            else
                logModia(i, ": ")      
            end
        end

        if vertical
            logModia(rpad(prettyfy(symbols[i]), 30, " "))
        else
            logModia(prettyfy(symbols[i]))
        end

        if A != [] 
            if A[i] != 0
                logModia("A[$i] = $(A[i])")
            end        
        end
    end

    if length(symbols) > maxSymbols
        logModia(", ...")
    end

    loglnModia()
end

function showModel(mod, indent="")
    loglnModia("showModel:")
    #=
    loglnModia("------")
    @show mod
    loglnModia("------")
    =#
    newIndent = string(indent, "  ")
    loglnModia("Instance(model ", mod.name)
    
    for i in 1:length(mod.initializers)
        if isa(mod.initializers[i], InitVariable)
            logModia(indent, mod.initializers[i].name, " = ")
            loglnModia(mod.initializers[i].init)
        else
            loglnModia(indent, mod.initializers[i])
            if isa(mod.initializers[i], Extends)
                # println("Extends($(mod.initializers[i].fdef.args[2].args[2].args[2].args[1]))")
                # Check that we have an Expr
                if isa(mod.initializers[i].fdef.args[2].args[2].args[2], Expr)
                    partial = mod.initializers[i].fdef.args[2].args[2].args[2].args[1] == :PARTIAL
                    if partial
                        # mod.partial = true
                    end
                else
                    error("Missing parentheses in extends clause: @extends $(mod.initializers[i].fdef.args[2].args[2].args[2])")
                end
            end
        end
    end
    loglnModia(indent, "end)\n")
end

function showVariable(v)
    if isa(v, Variable)
        logModia("Variable(")
        first = true

        if v.variability != continuous
            if !first; logModia(", ") end
            logModia("variability = ", v.variability)
            first = false
        end
        
        if v.typ != Any
            if !first; logModia(", ") end
            logModia("T = ", v.typ)
            first = false
        end
        
        if v.value != nothing
            if !first; logModia(", ") end
            logModia("value = ", v.value)
            first = false
        end
        
        if v.size != nothing
            if !first; logModia(", ") end
            logModia("size = ", v.size)
            first = false
        end
        
        if v.start != nothing
            if !first; logModia(", ") end
            logModia("start = ", v.start)
            first = false
        end
        
        if v.unit != NoUnits
            if !first; logModia(", ") end
            logModia("unit = ", v.unit)
            first = false
        end
        
        if v.min != nothing
            if !first; logModia(", ") end
            logModia("min = ", v.min)
            first = false
        end
        
        if v.max != nothing
            if !first; logModia(", ") end
            logModia("max = ", v.max)
            first = false
        end
        
        if v.flow != nothing && v.flow
            if !first; logModia(", ") end
            logModia("flow = ", v.flow)
            first = false
        end
        
        if v.state != true
            if !first; logModia(", ") end
            logModia("state = ", v.state)
            first = false
        end

        if v.info != ""
            if !first; logModia(", ") end
            logModia("info = \"", v.info, "\"")
            first = false
        end

        loglnModia(")")
    else
        loglnModia(v)
    end
end

function showInstance(inst, indent="")
    # loglnModia("showInstance:")
    #=
    loglnModia("------")
    @show inst
    loglnModia("------")
    =# 
    newIndent = string(indent, "  ")
    loglnModia("@model ", inst.model_name, " begin")
    for key in keys(inst.variables) 
        if isa(key, Instance)
            logModia(indent, "  ", key, " = ")
            showInstance(inst.variables[key], newIndent)
        else
            @static if VERSION < v"0.7.0-DEV.2005"
                keyname = replace(string(key), ".", "_")
            else
                keyname = replace(string(key), "." => "_")
            end
            logModia(indent, "  ", keyname, " = ") 
            showVariable(inst.variables[key])
        end
    end
    
    if length(inst.equations) > 0
        loglnModia(indent, "@equations begin")
    end
    
    for e in inst.equations
        loglnModia(indent, "  ", prettyPrint(e))
    end
    
    loglnModia(indent, "  end")
    loglnModia(indent, "end")
end

function checkSizes(VSizes, ESizes)
    # Check that system matrix is square
    ok = true
    if length(VSizes) != length(ESizes)
        ModiaLogging.increaseLogCategory(:UnbalancedModel)
        # For TestArray5
        #    error("The number of unknowns ($(length(VSizes))) is not equal to the number of equations ($(length(ESizes)))")
    end

    loglnModia()
    # Should check that arrays are not empty
    scalarV = sum(length(zeros(v)) for v in VSizes)
    scalarE = sum(length(zeros(e)) for e in ESizes)
    
    if scalarV != scalarE  
        # error("Scalarized system matrix is not square: $scalarE x $scalarV")
        ModiaLogging.closeLogModia()
        error("The number of scalarized unknowns (= $scalarV) is not equal to the number of scalarized equations (= $scalarE).\n",
              "If option `simulate(<model>, ...; logTranslation=true)` is set, inspect <user>/ModiaResults/<model>.txt for more info.")
        ok = false
    else
        # loglnModia("Scalarized system matrix is square: $scalarE x $scalarV")
        loglnModia("The number of scalarized unknowns (= $scalarV) is equal to the number of scalarized equations (= $scalarE).")
    end

    if false && (sort(VSizes) != sort(ESizes) )  ### Testing
        @show sort(VSizes) sort(ESizes)
        # Show the number of variables (and equations) of different sizes
        differentVSizes = unique(sort(VSizes))
        VSizeCount = [length([s for s in VSizes if s == d ]) for d in differentVSizes]
        differentESizes = unique(sort(ESizes))
        ESizeCount = [length([s for s in ESizes if s == d ]) for d in differentESizes]
        @show differentVSizes VSizeCount differentESizes ESizeCount
        ok = false
        error("The sizes of variables and equations don't match.")
    end
end


# ---------------------------------------------------------------

# Exprimental code for outputting connection structure to be used by Kieler automatic layout algorithm.

function isPort(inst)
    for (n, v) in inst.variables
        if isa(v, Variable) && v.flow
            return true
        end
    end
    return false
end

const size = 30
const sizePort = 3

function printJSON(file, inst::Instance, fullName, name, parent="", indent="", level=1)
    loglnModia(file, indent, "{")
    loglnModia(file, indent, "\"id\": \"$(fullName)\",")
    loglnModia(file, indent, "\"class\": \"$(inst.model_name)\",")
  
    params, unknowns = split_variables(vars_of(inst))
    if length(params) > 0
        loglnModia(file, indent, "\"parameters\": {")
        indent1 = "  " * indent
        i = 0
        for (n, p) in params
            i += 1
            logModia(file, indent1, "\"", n, "\"", " : ", p)
            if i < length(params)
                loglnModia(file, ",") 
            else
                loglnModia(file) 
            end
        end    
        loglnModia(file, indent, "},")
    end
  
    loglnModia(file, indent, "\"labels\": [{\"text\": \"$name\"}],")
    
    if isPort(inst)
        orientation = if length(search(string(parent), "G")) > 0; "\"NORTH\"" elseif length(search(string(name), "p")) > 0 || length(search(string(name), "in")) > 0; "\"WEST\"" else "\"EAST\"" end
        loglnModia(file, indent, "\"properties\": {\"de.cau.cs.kieler.portSide\": $orientation},")
        loglnModia(file, indent, "\"width\": $sizePort,")
        loglnModia(file, indent, "\"height\": $sizePort")
    else
        loglnModia(file, indent, "\"properties\": {\"de.cau.cs.kieler.portConstraints\": \"FIXED_SIDE\"},")
        loglnModia(file, indent, "\"width\": $size,")
        loglnModia(file, indent, "\"height\": $size")
    end
  
    first = true
    for (n, v) in inst.variables
        if isa(v, Instance) && isPort(v) && level > 1 # length(search(string(n), "n")) == 0
            if first
                loglnModia(file, indent, ", \"ports\": [")        
            else
                loglnModia(file, indent, ", ")
            end
            first = false
            printJSON(file, v, string(name) * string(n), n, name, "  " * indent, level + 1)
        end
    end

    if !first
        loglnModia(file, indent, "]")
    end

    first = true
    for (n, v) in inst.variables
        if isa(v, Instance) && (!isPort(v) || level == 1) #  || (isa(v, Instance) && isPort(v) && length(search(string(n), "n")) > 0)
            if first
                loglnModia(file, indent, ", \"children\": [")        
            else
                loglnModia(file, indent, ", ")
            end
            first = false
            printJSON(file, v, n, n, name, "  " * indent, level + 1)
        end
    end
    
    if !first
        loglnModia(file, indent, "]")
    end

    first = true
    id = 0
    for e in inst.equations
        if isa(e, Connect)
            if first
                loglnModia(file, indent, ", \"edges\": [") 
            else
                loglnModia(file, indent, ", ")
            end
    
            first = false
            id += 1
            loglnModia(file, indent, "  {\"id\": \"id$id\",")
    
            if isa(e.a.base, GetField)
                loglnModia(file, indent, "   \"source\": \"$(e.a.base.name)\",")
                loglnModia(file, indent, "   \"sourcePort\": \"$(string(e.a.base.name) * string(e.a.name))\",")
            else
                loglnModia(file, indent, "   \"source\": \"$(e.a.name)\",")        
            end
    
            if isa(e.b.base, GetField)
                loglnModia(file, indent, "   \"target\": \"$(e.b.base.name)\",")
                loglnModia(file, indent, "   \"targetPort\": \"$(string(e.b.base.name) * string(e.b.name))\"}")
            else
                loglnModia(file, indent, "   \"target\": \"$(e.b.name)\"}")        
            end
        end
    end
    
    if !first
        loglnModia(file, indent, "]")
    end

    loglnModia(file, indent, "}")
end


end

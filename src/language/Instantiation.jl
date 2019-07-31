"""
Modia module for instantiation and flattening of models.

* Original developer: Toivo Henningsson, Lund
* Developer: Hilding Elmqvist, Mogram AB
* Copyright (c) 2016-2019: Hilding Elmqvist, Toivo Henningsson, Martin Otter
* License: MIT (expat)

"""
module Instantiation

export Variability, constant, parameter, discrete, continuous
export Property, general, symmetric, orthogonal, rotationGroup3D
export Variable, RealType, Par, FLOAT
export @model, @equations, @equation
export VariableDict, Model
export simulationModel_global
export mark_solved_equations
export instantiate, flatten
export Connect, addEquation!, deleteEquation!
export prettyPrint, prettyfy, operator_table
export setOptions

global logMacros = false
global logInstantiation = false
global logFlattening = false
#const distributed = false

function setOptions(options) 
    global logMacros = false
    if haskey(options, :logMacros)
        global logMacros = options[:logMacros]
        @show logMacros
        delete!(options, :logMacros)
    end

    global logInstantiation = false
    if haskey(options, :logInstantiation)
        global logInstantiation = options[:logInstantiation]
        @show logInstantiation
        delete!(options, :logInstantiation)
    end

    global logFlattening = false
    if haskey(options, :logFlattening)
        global logFlattening = options[:logFlattening]
        @show logFlattening
        delete!(options, :logFlattening)
    end

end

# Desired:
#   using  DataStructures: OrderedDict
#   using  LinearAlgebra
#   import ModiaMath
#   using  Unitful
#   import Markdown
#
# Since package Instantiation is included directly in a test (in Modia/test/symbolic/BLTandPantelides/setup.jl)
# all these packages above are included via Modia, in order that these packages need not to be defined in the 
# user environment.
#
import Modia
using  Modia.DataStructures: OrderedDict
#0.7 using SparseArrays
@static if VERSION < v"0.7.0-DEV.2005"
    Nothing = Void 
    AbstractDict = Associative
    macro __MODULE__()
        return current_module()
    end
else
    using Modia.LinearAlgebra
end

import Modia.ModiaMath 
using  Modia.Unitful
using ..ModiaLogging
#using ..Synchronous


@static if VERSION >= v"0.7.0-DEV.2005"
    import Modia.Markdown
else
    import Base.Markdown
end
using Base.Meta:quot, isexpr
import Base.Docs

const shortSyntax = true

# Name that should stand for the current instance inside of model declaration
const this_symbol = :this

"Enum for variability of a variable"
@enum Variability constant = 1 parameter discrete continuous
"Indicates a constant valued variable"
constant
"Indicates an input value that stays constant with time"
parameter
"Indicates a continuous variable"
continuous
"Indicates a discrete variable"
discrete

@enum Property general = 1 symmetric orthogonal rotationGroup3D


FLOAT = Union{Float64,AbstractArray{Float64,1},AbstractArray{Float64,2}}


#"""
#The main object for a model variable with attributes
#"""
mutable struct Variable  # {T,n}
    variability::Variability
    typ
    size
    value #  ::T
    unit::Unitful.Unitlike
    displayUnit
    min
    max
    start
    fixed::Bool
    nominal
    info::AbstractString
    flow::Bool
    state::Bool
    property::Property
end

"""
A constructor for a `Variable`, the main object for a model variable with attributes

## Keyword arguments 

  * `value = undefined`: the value of the Variable
  * `info = ""`: documentation string
  * `unit`: unit of measure: defaults to no units unless provided with the `value`
  * `displayUnit = unit`: unit used for display
  * `min = undefined`: minimum value
  * `max = undefined`: maximum value
  * `start = undefined`: starting value
  * `fixed::Bool = false`: fixed value
  * `nominal = undefined`: nominal value
  * `variability = continuous`: other options include `parameter`, `constant`, or `discrete`
  * `T`: type of the value; taken from `value`
  * `size::Tuple = ()`: size of the variable
  * `flow::Bool = false`: indicates a flow variable for connectors
  * `state::Bool = true`: indicates a state Variable
  * `property = general`: other options include `symmetric`, `orthogonal`, and `rotationGroup3D`
"""
Variable(;
    # The variability, type and info are added as attributes in the type for uniform treatment.  
    # Input/output, etc should also be added.

    value=nothing, 
    info="", 
    unit=if typeof(value) <: Unitful.Quantity; Unitful.unit(value) else Unitful.NoUnits end, 
    displayUnit=unit, 
    min=nothing, 
    max=nothing, 
    start=nothing, 
    fixed::Bool=false, 
    nominal=nothing,
    variability=continuous, 
    T=if value == nothing; Any else typeof(value) end, 
    size=if value == nothing; nothing else Base.size(value) end, 
    flow::Bool=false, 
    state::Bool=true, 
    property=general) = 
    Variable(variability, T, size, value, 
    unit, displayUnit, min, max, start, fixed, nominal, info, flow, state, property)
    
function Base.show(io::IO, v::Variable)
    print(io, "Variable(")
    first = true
    
    if v.value != nothing
        if !first; print(io, ", ") end
        print(io, "value = ", v.value)
        first = false
    end

    if v.variability != continuous
        if !first; print(io, ", ") end
        print(io, "variability = ", v.variability)
        first = false
    end
    
    if v.typ != Any
        if !first; print(io, ", ") end
        print(io, "T = ", v.typ)
        first = false
    end
    
    if v.size != nothing
        if !first; print(io, ", ") end
        print(io, "size = ", v.size)
        first = false
    end
    
    if v.start != nothing
        if !first; print(io, ", ") end
        print(io, "start = ", v.start)
        first = false
    end
    
    if v.fixed != nothing
        if !first; print(io, ", ") end
        print(io, "fixed = ", v.fixed)
        first = false
    end

    if v.nominal != nothing
        if !first; print(io, ", ") end
        print(io, "nominal = ", v.nominal)
        first = false
    end
    
    if v.unit != NoUnits
        if !first; print(io, ", ") end
        print(io, "unit = ", v.unit)
        first = false
    end
    
    if v.min != nothing
        if !first; print(io, ", ") end
        print(io, "min = ", v.min)
        first = false
    end
    
    if v.max != nothing
        if !first; print(io, ", ") end
        print(io, "max = ", v.max)
        first = false
    end
    
    if v.flow != false
        if !first; print(io, ", ") end
        print(io, "flow = ", v.flow)
        first = false
    end
    
    if v.state != true
        if !first; print(io, ", ") end
        print(io, "state = ", v.state)
        first = false
    end
    
    println(io, ")")
end


"Check that a start value (possibly default) exists for the var, or give an error."
function check_start(var::Variable, name)
    if (var.start === nothing) && (var.typ !== Nothing) && !applicable(zero, var.typ)
        error("Variable ", name, " has no start value and no default exists for type ", var.typ)
    end
end

"Get a default value for a given type `T`."
default(::Type{Nothing}) = 0
default(T::Type) = zero(T)

"Get the start value of `var`. Returns a default value if no start value is given."
get_start(var::Variable) = var.start === nothing ? default(var.typ) : var.start

"Get the array dimensions of `var`."
get_dims(var::Variable) = size(get_start(var))


# ---------------------------- @model, @equations -----------------------------

"Base type for model specific AST node types."
abstract type Symbolic end

"AST node for a field access to a field of a `Symbolic` node."
struct GetField <: Symbolic
    base::Symbolic
    name::Symbol
end

struct Ref <: Symbolic
    base::Symbolic
    inds::Vector
end

function Base.show(io::IO, g::GetField)
    if shortSyntax
        print(io, g.base, ".", g.name)  
    else
        print(io, "GetField(", g.base, ", ", g.name, ")")
    end
end

"AST node for the derivative of a `Symbolic` node."
struct Der <: Symbolic
    base::Symbolic
end

function Base.show(io::IO, d::Der)
    if shortSyntax
        print(io, "der(", d.base, ")")
    else
        print(io, "Der(", d.base, ")")  
    end
end

"AST node for access to the current instance."
struct This <: Symbolic
end

function Base.show(io::IO, this::This)
    if shortSyntax
        print(io, "this")
    else
        print(io, "This()")  
    end
end

"Representation of a connect equation, used in the equations list."
struct Connect
    a::Symbolic
    b::Symbolic
end

function Base.show(io::IO, c::Connect)
    print(io, "connect(", c.a, ", ", c.b, ")")
end

"""
AST node for special variables that are global to the models.

Unlike instance variables, global variables refer to the same variable
no matter in which instance the reference occurs.
"""
struct Global <: Symbolic
    name::Symbol
end

const time_global = Global(:time)
const simulationModel_global = Global(:simulationModel)

function Base.show(io::IO, g::Global)
    print(io, g.name)
end

is_linenumber(ex) = isa(ex, LineNumberNode) || isexpr(ex, :line)

"""
Field access function used for equations in `@model`.

Creates an AST node for the field access:
* For field access to a `Symbolic` node, create a `GetField` node (deferred access).
* For field access to other values, look up now and quote the result.
"""
function model_getfield(x, name)
    @show x
    @show name
    error()
end
model_getfield(x, name::Symbol) = selquote(getfield(x, name))
model_getfield(x::Symbolic, name::Symbol) = GetField(x, name)

#model_ref(x, inds::Vector) = selquote(x[inds...])
#model_ref(x::Symbolic, inds::Vector) = Ref(x, inds)
#= Multi body no units does not work:
model_ref(x, ind) = selquote(x[ind])
model_ref(x::Symbolic, ind) = Ref(x, [ind])
=#

"Quote an expression; don't quote if the meaning in the AST is the same anyway."
selquote(x) = x
selquote(sym::Symbol) = quot(sym)
selquote(ex::Base.ExprNode) = quot(ex)

"""
Recode an AST node in an equation.

Return code that evaluates to a new AST node:
* Preserve expression structure.
* Look up identifiers and replace them with their values.
* Recode field access using `model_getfield`.
* Create `Der` and `Connect` nodes from calls to `der()` and `connect()`.
"""
recode(ex) = :( $(quot(selquote))($ex) )
function recode(ex::Expr)
    if is_linenumber(ex); return ex; end
    head, args = ex.head, ex.args
    if head == :block
        args = collect(filter(ex -> (!is_linenumber(ex)), args))
#=
        println("recode:")
        dump(args)

        rex = Expr(:block)
        for a in args
            ra = recode(a)
            println("ra:")
            dump(ra)
            push!(rex.args, ra)
        end
        println("rex:")
        dump(rex)
        return rex
=#
        if length(args) == 1;  return recode(args[1]);  end
    elseif head == :.
        return :( $(quot(model_getfield))($(recode(args[1])), $(args[2])) )
    #= Problem with multi body no units
    elseif head == :ref
        return :( $(quot(model_ref))($(recode(args[1])), $(recode(args[2])) ))
    =#
    elseif head == :kw
        @assert length(args) == 2
        arg1, arg2 = quot(args[1]::Symbol), recode(args[2])
        return :( $(quot(Expr))($(quot(:kw)), $arg1, $arg2) )
    elseif head == :quote
        return quot(ex)
    elseif head == :call && length(args) == 2 && args[1] == :der
        return :( $(quot(Der))($(recode(args[2]))) )
    elseif head == :$
        return args[1]
    elseif head == :call && args[1] == :connect
        if length(args) == 3
            return :( $(quot(Connect))($(recode(args[2])), $(recode(args[3]))) )
        else
            print("The connect statement takes two arguments: $ex")
        end
    end

    args = map(recode, ex.args)
    :( $(quot(Expr))($(quot(head)), $(args...)) )
end

"Recode and return the list of equations contained in the AST node `eqs`."
function recode_equations(eqs)
    equations = []
    @assert isexpr(eqs, :block)
    for eq in eqs.args
        if is_linenumber(eq)
            continue
        elseif isexpr(eq, :for)
            println("recode_equation: for found")
            dump(eq)
            return
#=
        elseif isexpr(eq, :if, 3)  # need to handle missing else and elseif
            println("if equation found")
            dump(eq)
            cond = eval(eq.args[1])
            @show cond
            if cond
                push!(equations, recode(eq.args[2]))
            else
                push!(equations, recode(eq.args[3]))
            end       
=#            
        else
            push!(equations, recode(eq))
        end
    end
    equations
end

"Implementation of the standalone @equations macro. (no loner used)"
function code_equations(eqs)
    equations = recode_equations(eqs)
    quote
        let $(this_symbol) = $(quot(This()))
            [$(equations...)]
        end
    end
end

"Return the list of preprocessed ASTs for the supplied equations."
macro equations(eqs)
    esc(code_equations(eqs))
end

"Implementation of the standalone @equations macro. (no loner used)"
function code_equation(eq)
    equation = recode(eq)
    # @show equation
    quote
        let $(this_symbol) = $(quot(This()))
            $(equation)
        end
    end
end

"Return the list of preprocessed ASTs for the supplied equations."
macro equation(eq)
    esc(code_equation(eq))
end

"""
Recode an AST node to be used in a variable initializer/`@extends`.

Return a new AST node. Field accesses are replaced with calls to `initializer_getfield`.
"""
recode_initializer(ex) = ex
function recode_initializer(ex::Expr)
    head, args = ex.head, ex.args

    if head === :.
        @assert length(args) == 2
        @assert isa(args[2], QuoteNode) || isexpr(args[2], :quote)
        :( $(quot(initializer_getfield))($(recode_initializer(args[1])), $(args[2])) )
    else
        Expr(head, [recode_initializer(arg) for arg in args]...)
    end
end

"Get the variable name declared in a variable declaration AST node."
get_declared_var_name(name::Symbol) = name
get_declared_var_name(ex::Expr) = (if ! isexpr(ex, :(=), 2); dump(ex); end; @assert isexpr(ex, :(=), 2); ex.args[1])

"""
Create code that initializes all variables that are known as local to the instance.

The list of variables should be gathered from the declarations local to the model
and any `@inherits` declarations. For each name, `varnames` should indicate whether
the corresponding variable has been initialized yet:
* Initialized variables are read from the instance through `this` (`$this_symbol`).
* Not yet initialized variables are bound as uninitialized locals to trigger
  an error if the code tries to use them.
"""
function code_init_locals(varnames::AbstractDict{Symbol,Bool})
    if isempty(varnames);  return quote end;  end
    Expr(:local,
        [(init ? :($name = $(quot(initializer_getfield))($this_symbol, $(quot(name)))) : name)
            for (name, init) in varnames]...)
end

"Recode a variable declaration based on the corresponding AST node from a `@model` invocation."
function code_variable(name::Symbol, varnames)
    # Variable without initializer. Use the `no_initializer` functions as initializer,
    # which will throw an error at initialization unless a modifier is used to provide
    # another initializer.
    fdef = :( (this, time) -> $(quot(no_initializer))(this, $(quot(name))) )
    :( $(quot(InitVariable))($(quot(name)), $fdef, $(QuoteNode(fdef))) )
end

function code_variable(ex::Expr, varnames)
    @assert isexpr(ex, :(=), 2)
    lhs, rhs = ex.args
#    if typeof(rhs) == Expr && (rhs.head == call || rhs.head == :call)
    if typeof(rhs) == Expr && rhs.head == :call
        args = rhs.args
        if length(args) > 0 && typeof(args[1]) == Expr && args[1].head == :curly
            first = args[1].args[2]
#        if isa(eval(typ),DataType)
            if !isa(first, Int64)
                typ = first
                siz = args[1].args[3:end]
            else
                typ = Any
                siz = args[1].args[2:end]
            end
#        @show typ
            baseTyp = typ

            if siz != []
                typ = :(Array{$typ,$(length(siz))})
                  
            end
            if baseTyp in [:Int64, :Float64, :Bool]
                sta = Expr(:kw, :start, zeros(eval(baseTyp), siz...))
                rhs = Expr(:call, :Variable, Expr(:kw, :typ, typ), sta, args[2:end]...)
            else
                rhs = Expr(:call, :Variable, Expr(:kw, :typ, typ), args[2:end]...)
            end
        end
    end

    rhs = recode_initializer(rhs)
    locals = code_init_locals(varnames)
    fdef = :( ($this_symbol, time) -> ($locals; $rhs) )
    :( $(quot(InitVariable))($(quot(lhs)), $fdef, $(QuoteNode(fdef))) )
end

"Recode an `@extends` declaration based on the corresponding AST node from a `@model` invocation."
function code_extends(ex, varnames)
    locals = code_init_locals(varnames)
    fdef = :( ($this_symbol, time) -> ($locals; $(recode_initializer(ex))) )
    :( $(quot(Extends))($fdef, $(QuoteNode(fdef))) )
end

"Recode an `@equations` declaration based on the corresponding AST node from a `@model` invocation."
code_eqs(ex) = :( $(quot(Equations))([$(recode_equations(ex)...)]) )

"Parse an `@inherits` declaration an update `varnames`"
parse_inherits!(varnames::AbstractDict, name::Symbol, init::Bool) = (varnames[name] = init; nothing)
function parse_inherits!(varnames::AbstractDict, ex::Expr, init::Bool)
    @assert ex.head === :tuple
    for name in ex.args;  varnames[name] = init;  end
end

# Symbols that we need to recognize and can't quote with :
const extends_symbol = Symbol("@extends")
const inherits_symbol = Symbol("@inherits")
const equations_symbol = Symbol("@equations")

# Symbols that should be used consistently for unique variables in the generated code
const time_symbol = gensym("time")
const simulationModel_symbol = gensym("simulationModel")

@static if VERSION < v"0.7.0-DEV.2005"
    const macroCallNumberOfArguments = 2
else
    const macroCallNumberOfArguments = 3
end

"Implementation of `@model`"
function code_model(head, top_ex)
    #=
    println("---------------------")
    println("code_model")
    @show name top_ex
    println("---------------------")
    @assert isexpr(top_ex, :block)
    =#
    # Extract name and any model arguments
    if isa(head, Symbol)
        name = head
        args = []
    else
        @assert isexpr(head, :call)
        name = head.args[1]::Symbol
        args = head.args[2:end];
    end

    # Gather local variable names - short version of loop below
    varnames = OrderedDict{Symbol,Bool}() # variable name => initialized?
    for ex in top_ex.args
            if is_linenumber(ex); continue
        elseif isexpr(ex, :macrocall, macroCallNumberOfArguments)
            if ex.args[1] == inherits_symbol
                parse_inherits!(varnames, ex.args[macroCallNumberOfArguments], false)
            end
        elseif isexpr(ex, :if, 2)
            println("if found")
            dump(ex)
            cond = eval(ex.args[1])
            @show cond
            if cond
                varnames[get_declared_var_name(ex.args[2].args[2])] = false
            end        
        else
            varnames[get_declared_var_name(ex)] = false
        end
    end

    # Create code for initializers
    initializers = []
    for ex in top_ex.args
        if is_linenumber(ex)
            continue
        elseif isexpr(ex, :macrocall, macroCallNumberOfArguments) && ex.args[1] == equations_symbol
            @assert length(ex.args) == macroCallNumberOfArguments
            push!(initializers, code_eqs(ex.args[macroCallNumberOfArguments]))
        elseif isexpr(ex, :macrocall, macroCallNumberOfArguments) && ex.args[1] == extends_symbol
            @assert length(ex.args) == macroCallNumberOfArguments
            push!(initializers, code_extends(ex.args[macroCallNumberOfArguments], varnames))
        elseif isexpr(ex, :macrocall, macroCallNumberOfArguments) && ex.args[1] == inherits_symbol
            @assert length(ex.args) == macroCallNumberOfArguments
            parse_inherits!(varnames, ex.args[macroCallNumberOfArguments], true)
        elseif isexpr(ex, :if, 2)
            println("if found 2")
            dump(ex)
            cond = eval(ex.args[1])
            @show cond
            if cond
                varnames[get_declared_var_name(ex.args[2].args[2])] = true
            end        
        else
            initvar_ex = code_variable(ex, varnames)
            push!(initializers, initvar_ex)
            varnames[get_declared_var_name(ex)] = true
        end
    end

    var_bindings = [:( $(name) = $(quot(GetField(This(), name))) ) for name in keys(varnames)]
    push!(var_bindings, :(time = $(quot(time_global))))
    push!(var_bindings, :(simulationModel = $(quot(simulationModel_global))))

    quote
        $(name) = (let $(this_symbol) = $(quot(This())), $(var_bindings...)
#            $(quot(Model))($(quot(name)), [$(initializers...)])
            $(quot(Model))($(quot(name)), Modia.Instantiation.@__MODULE__, [$(args...)], [$(initializers...)])
        end)
    end
end

"""
    @model <Name> begin
        <declarations>
        @equations begin
            <equations>
        end
    end

Fill in a `Model` instance with the given declarations and equations and assign it
to a constant named <Name>.
"""
macro model(head, ex)
    coded = code_model(head, ex)
    # Resolve code in macro environment
    esccoded = esc(coded)
    
    if logMacros
        println("@model:", head)
        @show Base.remove_linenums!(ex)
        println()
        dump(Base.remove_linenums!(ex), maxdepth=100)
        println("------------------")
        println("coded:")
        @show Base.remove_linenums!(coded)
        println()
        dump(Base.remove_linenums!(coded), maxdepth=100)
        println("------------------")
        println("Resolved (esc(coded)):")
        @show Base.remove_linenums!(esccoded)
        println()
        dump(Base.remove_linenums!(esccoded), maxdepth=100)
        println("------------------")
    end

    return esccoded
end


# ------------------------------- Initializer ---------------------------------

"""
Base class for initializers.

An initializer is a step for filling in an instance.
The contents of a Model is described as a list of initializers.
"""
abstract type Initializer end

"""
Initializer for a named variable with given default initializer function.

The `init` function should be called with arguments `(this,time)`, where `this`
is the partially filled instance and `time` is the current time.
The function should return the value to be bound to the variable.
"""
struct InitVariable <: Initializer
    name::Symbol
    init::Function
    fdef # AST of init, for display
end

function Base.show(io::IO, iv::InitVariable)
    print(io, "InitVariable(", repr(iv.name), ", ", iv.fdef, ")")
end

"""
Initializer for an `@extends` clause.

The `init` function should work the same as for `InitVariable`, but return an Instantiations
for an instance of a model that the current instance will inherit the contents from.
"""
struct Extends <: Initializer
    init::Function
    fdef # AST of init, for display
end

function Base.show(io::IO, e::Extends)
    print(io, "Extends(", e.fdef, ")")
end

"Initializer to add equations to an instance."
struct Equations <: Initializer
    equations::Vector
end


# ---------------------------------- Model ------------------------------------

const VariableDict = OrderedDict{Symbol,Any}

"""
    type Model

A `Model` object is a description on how to fill in an `Instance` object.
"""
mutable struct Model
    name::Symbol
    mod::Module
    arguments::Vector
    initializers::Vector{Initializer}
end

#=
function Base.show(io::IO, model::Model)
    #println(io, "model ", model.name)
    println(io, "Model(", repr(model.name), ", [")
    for init in model.initializers
        print(io, "    ", init, ",\n")
    end
    #println(io, "end ", model.name)
    println(io, "])")
end
=#
"Get the mode name of a `Model` or `Instance`."
model_name_of(model::Model) = model.name

"""
A model to be instantiated with given arguments.

Used as an intermediate representation to allow to apply modifiers before a model
is actually instantiated. When an `Instantiations` is stored as value for a variable of
an `Instance`, it will be instantiated.
"""
mutable struct Instantiations
    model::Model
    kwargs
end
Instantiations(model::Model) = Instantiations(model, [])

# Calling a `Model` (optionally with keyword arguments) creates an Instantiations.
# call(model::Model; kwargs...) = Instantiations(model, kwargs)
(model::Model)(;kwargs...) = Instantiations(model, kwargs)

"Instance of a model, with variable bindings and equations."
mutable struct Instance
    model_name::Symbol
    mod::Module
    info::String
    variables::VariableDict
    equations::Vector{Any}
    partial::Bool

    # Extra code to be inserted raw before and after the computations (experimental)
    initial_pre::Vector{Any}
    initial_post::Vector{Any}
    F_pre::Vector{Any}
    F_post::Vector{Any}
end

function Instance(model_name::Symbol, mod, info, variables, equations, partial)
    Instance(model_name, mod, info, VariableDict(variables), 
        collect(Any, equations), partial, [], [], [], [])
end

#=
function Base.show(io::IO, inst::Instance)
  println(io, "Instance(")
  println(io, "  name = $(inst.model_name),")
  println(io, "  variables = $(inst.variables)")
  if length(inst.equations) > 0
    println(io, "  equations = Any[")
    for e in inst.equations
      print(io, "    ")
      print(io, ":(", prettyPrint(e), ")")
#      print(io, ":(", e, ")")
      println(io)
    end
    println(io, "  ]")
  end

  println(io, ")")
end
=#


"Get the variables of an `Instance` as a `VariableDict`."
vars_of(instance::Instance) = instance.variables
"Get an iterable of the equations of an `Instance`."
eqs_of(instance::Instance)  = instance.equations
model_name_of(instance::Instance) = instance.model_name

"""
Field access used in initializer code.

Look up variables in `Instance`s, otherwise like regular `getfield`.
"""
function initializer_getfield(x, name)
    @show x
    @show name
    error()
end

initializer_getfield(x, name::Symbol) = getfield(x, name)

function initializer_getfield(instance::Instance, name::Symbol)
    if !haskey(instance.variables, name)
        ModiaLogging.increaseLogCategory(:NoFieldDefined)
        error("No field ", name, " defined in model ", instance.model_name, " (yet)")
    end
    instance.variables[name]
end

"Give an error message that no initializer has been defined for a varible."
function no_initializer(inst::Instance, name::Symbol)
    error("No initializer given for variable $(model_name_of(inst)).$name")
end

# -------------------------------- instantiate --------------------------------

# Pass through anything that's not an Instantiations
as_field_value(x, time::Float64) = x
as_field_value(inst::Instantiations, time::Float64) = instantiate(inst.model, time, inst.kwargs)
as_field_value(insts::Vector{Instantiations}, time::Float64) = Instance[instantiate(inst.model, time, inst.kwargs) for inst in insts]

function add_variable!(instance::Instance, name::Symbol, value)
    if haskey(instance.variables, name)
        # error("Multiple definitions of $(name) in $(model_name_of(instance))")
        println("Multiple definitions of $(name) in $(model_name_of(instance))")
    end
    instance.variables[name] = value
end

function initialize!(instance::Instance, iv::InitVariable, time::Float64, kwargs::AbstractDict)
    # @show iv.fdef
    add_variable!(instance, iv.name, 
        as_field_value(haskey(kwargs, iv.name) ? kwargs[iv.name] : iv.init(instance, time),
                       time))
end

function initialize!(instance::Instance, ext::Extends, time::Float64, kwargs::AbstractDict)
    base_instantiation = ext.init(instance, time)::Instantiations
    base_instantiation.kwargs = merge!(Dict{Symbol,Any}(base_instantiation.kwargs), kwargs)
    base = as_field_value(base_instantiation, time)

    for (name, value) in base.variables
        add_variable!(instance, name, as_field_value(value, time))
    end
    append!(instance.equations, base.equations)
end

function instantiate_equation!(instance::Instance, eq)
    if !isexpr(eq, :if) 
        push!(instance.equations, eq)
        return
    end
    #=
    if isexpr(eq, :for)
        dump(eq)
        return
    end
    =#
    
    println("Conditional equation:")
    println(prettyPrint(eq)) 
    if length(eq.args) > 3 
        error("elseif is presently not handled.")
    end
    cond = eq.args[1]
    if typeof(cond) == GetField # only handle name that resolves in the model for now
        cond_value = lookup(instance, cond)
        if typeof(cond_value) == Variable
            cond_value = cond_value.value
        end
    else
        op = cond.args[1]
        if typeof(cond.args[2]) == GetField 
            cond_value = lookup(instance, cond.args[2])
            if typeof(cond_value) == Variable
              cond_value = cond_value.value
            end
        else
            error("Too complex expression: $cond")
        end
        if op == !
            cond_value = ! cond_value
        else
            error("Not handled operator.")
        end
    end
    println("condition = ", cond_value)
    
    eq_index = cond_value ? 2 : 3
    if eq_index <= length(eq.args)
      branch_eq = eq.args[eq_index]
      if ! (branch_eq.head in [:(=), :block])
          error("At most one equation is currently allowed in conditional equation.")
      end
      
      if branch_eq.head == :(=) 
          instantiate_equation!(instance, branch_eq)
      end
    end
end

function initialize!(instance::Instance, eqs::Equations, time::Float64, kwargs::AbstractDict)
    #append!(instance.equations, eqs.equations)
    for eq in eqs.equations
        instantiate_equation!(instance, eq)
    end
end

# using Distributed
# using SharedArrays

function instantiate(model::Model, time::Float64, kwargs=[])
    kwargs = Dict(kwargs)
    if logInstantiation
        println("\ninstanciate:::::::")
        @show model
        println("dump(model):")
        dump(model, maxdepth=100)
        @show kwargs
    end
    instance = Instance(model_name_of(model), model.mod, get(kwargs, :info, ""), 
			VariableDict(), [], :partial in model.arguments)
    if logInstantiation
        println("dump(instance):")
        dump(instance, maxdepth=100)
    end
#    if ! distributed
        for initializer in model.initializers
            if logInstantiation
                println("dump(initializer):")
                dump(initializer, maxdepth=100)
            end
            initialize!(instance, initializer, time, kwargs)
        end
#=
    else # This is not working yet
        insts = SharedArray{Instance,1}(length(model.initializers))
        @distributed for i in 1:length(model.initializers)
            insts[i] = []
            initializer = model.initializers[i]
            initialize!(insts[i], initializer, time, kwargs)
            @show insts[i]
        end    
        @show insts
        instance = vcat(insts)
    end
=#
    instance
end

# --------------------------------- flatten ----------------------------------

flatten_this(ex, prefix::AbstractString) = ex

function flatten_this(ex::Expr, prefix::AbstractString)
    Expr(ex.head, [flatten_this(arg, prefix) for arg in ex.args]...)
end

flatten_this(ex::Der, prefix::AbstractString) = Der(flatten_this(ex.base, prefix))

function flatten_this(ex::GetField, prefix::AbstractString)
    if ex.base == This()
        GetField(This(), Symbol(prefix, ex.name))
    else
        base = flatten_this(ex.base, prefix)::GetField
        @assert base.base === This()
        GetField(This(), Symbol(base.name, ".", ex.name))
    end
end

lookup(instance::Instance, ::This) = instance

function lookup(instance::Instance, g::GetField)
    instance = lookup(instance, g.base)::Instance
    var = vars_of(instance)[g.name]
    var
end

function lookup(instance::Instance, r::Ref)
    #@show instance
    #@show r
#error()
    lookup(instance, r.base)[r.inds...]
end

get_flow(var) = false
get_flow(var::Variable) = var.flow

#get_connector_type(inst::Instance) = [name => get_flow(var) for (name, var) in vars_of(inst)]
get_connector_type(inst::Instance) = Dict(name => get_flow(var) for (name, var) in vars_of(inst))
get_connector_type(var::Variable) = Dict(nothing => get_flow(var)) # consider: better representation?

const Connection = Pair{Symbol,Bool}

struct Flat
    vars::VariableDict
    eqs::Vector{Any}
    connections::Dict{Connection,Connection}
    connection_types::Dict{Connection,Dict}
end

function get_representative(flat::Flat, c::Connection)
    if haskey(flat.connections, c)
        flat.connections[c] = get_representative(flat, flat.connections[c])
    else
        c
    end
end


# --------------------------------- connect ----------------------------------

function connect!(flat::Flat, a::Connection, b::Connection, ctype::Dict)
    a = get_representative(flat, a)
    b = get_representative(flat, b)
    if a != b
        flat.connections[b] = a
        flat.connection_types[a] = ctype
    end
end

get_this_fieldname(g::GetField) = (@assert g.base === This(); g.name)
get_connection_sign(ex::This) = false

function get_connection_sign(ex::GetField)
    if !isa(ex.base, This); error("Can only connect current model and submodels"); end
    true
end

function get_connection(ex::GetField, prefix::AbstractString)
    Connection(get_this_fieldname(flatten_this(ex, prefix)), get_connection_sign(ex.base))
end

function add_connection!(flat::Flat, prefix::AbstractString, instance::Instance, eq::Connect)
    # Check if connector exist in instance
    inst = lookup(instance, eq.a.base)
    if inst === nothing; return; end
    inst::Instance
    
    if eq.a.name in keys(vars_of(inst))
        atype = get_connector_type(lookup(instance, eq.a))
    else
        error("Connector $(eq.a.name) not found in: $eq in model $(instance.model_name)")
    end

    inst = lookup(instance, eq.b.base)
    if inst === nothing; return; end
    inst::Instance

    if eq.b.name in keys(vars_of(inst))
        btype = get_connector_type(lookup(instance, eq.b))
    else
        error("Connector $(eq.b.name) not found in: $eq in model $(instance.model_name)")
    end

    if atype != btype
        ModiaLogging.increaseLogCategory(:IncompatibleConnect)
    end

    @assert atype == btype "Ports in connect statement $eq are not of same type: $atype, $btype"

    aconn = get_connection(eq.a, prefix)
    bconn = get_connection(eq.b, prefix)

    connect!(flat, aconn, bconn, atype)
end

function flatten!(flat::Flat, prefix::AbstractString, name::Symbol, var)
    flat.vars[Symbol(prefix, name)] = var
    nothing
end

function flatten!(flat::Flat, prefix::AbstractString, name::Symbol, inst::Instance)
    flatten!(flat, string(prefix, name, "."), inst)
end

function flatten!(flat::Flat, prefix::AbstractString, name::Symbol, insts::Vector{Instance})
    for (k, inst) in enumerate(insts)
        flatten!(flat, string(prefix, name, "[$k]."), insts[k])
    end
end

function flatten!(flat::Flat, prefix::AbstractString, instance::Instance)
    for (name, var) in vars_of(instance)
        if logFlattening
            println(name)
        end
        flatten!(flat, prefix, name, var)
    end
    for eq in eqs_of(instance)
        if logFlattening
            println(prettyPrint(eq))
        end
        if isa(eq, Connect)
            add_connection!(flat, prefix, instance, eq)
        else
            push!(flat.eqs, flatten_this(eq, prefix))
        end
    end
end

to_access(connector::Symbol, field::Symbol) = GetField(This(), Symbol(connector, ".", field))
to_access(connector::Symbol, ::Nothing) = GetField(This(), Symbol(connector))

function to_access(connector::Connection, field::Union{Symbol,Nothing})
    ex = to_access(connector.first, field)
    connector.second ? (:( $(-)($ex) )) : ex
end

function flatten(instance::Instance)
    flat = Flat(VariableDict(), [], 
        Dict{Symbol,Symbol}(), Dict{Symbol,Dict}())
    flatten!(flat, "", instance)

    # Gather connection sets
    connection_sets = Dict{Connection,Set{Connection}}()
    for node in keys(flat.connections)
        rep = get_representative(flat, node)
        set = get!(() -> Set{Connection}([rep]), connection_sets, rep)
        push!(set, node)
    end

    @static if VERSION < v"0.7.0-DEV.2005"
        unconnected_flow_vars = VariableDict(filter((name, var) -> (isa(var, Variable) && var.flow), flat.vars))
    else
        unconnected_flow_vars = VariableDict(filter(p::Pair -> (isa(p.second, Variable) && p.second.flow), flat.vars))
    end
    # Create connection equations
    for (rep, set) in connection_sets
        ctype = flat.connection_types[rep]
        for (name, flow) in ctype
            if flow
                push!(flat.eqs, Expr(:(=), Expr(:call, +,
                    [to_access(node, name) for node in set]...), 0))
                for node in set
                    if node.second
                        delete!(unconnected_flow_vars, Symbol(node.first, ".", name))
                    end
                end
            else
                for node in set
                    if node === rep; continue; end
                    push!(flat.eqs, Expr(:(=), to_access(node.first, name), to_access(rep.first, name)))
                end
            end
        end
    end

    # Set flow variables that are unconnected from the outside to zero
    for (name, var) in unconnected_flow_vars
        dims = get_dims(var)
        z = dims == () ? 0 : zeros(dims)
        push!(flat.eqs, :($(GetField(This(), name)) = $z))
    end

    Instance(model_name_of(instance), instance.mod, instance.info, flat.vars, flat.eqs, instance.partial)
end


# -----------------------------------

addEquation!(M, e) = begin push!(M.initializers, Equations([e])) end

function deleteEquation!(M, e)
    M.initializers = filter!(i -> if typeof(i) == Equations; i.equations[1] != e else true end, M.initializers)
end


# --------------------------------- utilities ----------------------------------


# prettyfy(der::Der) = Symbol("der("*string(der.base.name)*")")
@static if VERSION < v"0.7.0-DEV.2005"
    prettyfy(der::Der) = Symbol("der(" * replace(string(der.base.name), ".", "_") * ")")
    prettyfy(get::GetField) = Symbol(replace(string(get.name), ".", "_")) # get.name # Handle dummy derivatives
else
    prettyfy(der::Der) = Symbol("der(" * replace(string(der.base.name), "." => "_") * ")")
    prettyfy(get::GetField) = Symbol(replace(string(get.name), "." => "_")) # get.name # Handle dummy derivatives
end
# prettyfy(get::GetField) = get.name
# prettyfy(s::Symbol) = s
prettyfy(ex) = ex

function prettyfy(ex::Expr)
    if isexpr(ex, :quote) || isexpr(ex, :line)
        nothing
    elseif isexpr(ex, :block)
#        if length(ex.args) >=2 # need to handle emtpy else 
          prettyfy(ex.args[2])
#        end
    else
        Expr(ex.head, [prettyfy(arg) for arg in ex.args]...)
    end
end

# Pretty printing of expressions
const oper = [:!, :(!=), :(!==), :%, :&, :*, :+, :-, :/, ://, :<, :<:, :<<, :(<=),
               :<|, :(==), :(===), :>, :>:, :(>=), :>>, :>>>, :\, :^, #= :colon, =#
               :ctranspose, :getindex, :hcat, :hvcat, :setindex!, :transpose, :vcat,
               :xor, :|, :|>, :~ #= , :× =# , :÷, :∈, :∉, :∋, :∌, :∘, :√, :∛, :∩, :∪, :≠, :≤,
               :≥ #=, :⊆, :⊈, :⊊, :⊻, :⋅=#]
               
const operator_table = Dict(getfield(Base,name) => name for name in
    filter(name->isdefined(Base,name), oper))

prettyPrint(ex) = get(operator_table, ex, ex)
Array{Any}

function prettyPrint(e::Expr)
    ex = prettyfy(e)
    if ex.head === :quote
        return ex
    elseif ex.head === :(:=)
        return string(prettyPrint(ex.args[1]), " := ", prettyPrint(ex.args[2]))
    end
    Expr(ex.head, [prettyPrint(arg) for arg in ex.args]...)
end


formatvar(v::Variable) = string("variable : ", v.info, " [", v.typ, "]")
formatvar(v::Instance) = string("model : ", v.model_name, v.info != "" ? string(" : ", v.info) : "")
formatvar(v) = string(typeof(v))

"Dynamic documentation for Models"
function Docs.getdoc(m::Model, args = [])
    docstr = try
        join(Docs.docstr(Docs.Binding(m.mod, m.name)).text, "\n")
    catch
        ""
    end
  
    #vars = instantiate(m, 0.0, args).variables
    vars = instantiate(m, 0.0).variables
    if length(vars) > 0
        docstr *= "\n#### Variables\n\n"
    end
    for (name, v) in vars
        docstr *= "* `$name` : "
        docstr *= formatvar(v)
        docstr *= "\n"
    end
    
    Markdown.parse(docstr)
end


end 

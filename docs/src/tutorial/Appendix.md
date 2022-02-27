# Appendix A

## A.1 Var constructor

The constructor `Var(..)` defines attributes of a variable with key/value pairs.
In column 1 the keys are shown. The default is that none of the keys are defined
(meaning `key = nothing`). Most of the keys are also provided as predefined constants as shown
in column 2 and 3. These constants can be used as shortcuts:

| Var key    | ShortCut  | Shortcut value        |  Description                                       |
|:---------- |:----------|:----------------------|:---------------------------------------------------|
| parameter  | parameter | Var(parameter = true) | If true, value is fixed during simulation          |
| input      | input     | Var(input = true)     | If true, input signal                              |
| output     | output    | Var(output = true)    | If true, output signal                             |
| potential  | potential | Var(potential = true) | If true, potential variable                        |
| flow       | flow      | Var(flow = true)      | If true, flow variable                             |
| init       | --        | --                    | Initial value of ODE state (defines unit and size) |
| start      | --        | --                    | Start value of variable (defines unit and size)    |
| hideResult | --        | --                    | If true, the variable is not stored in the result  |


Example:

```julia
v = output | Var(start = zeros(3)u"N*m")

# Same as: v = Var(output = true, start = zeros(3)u"N*m")
```

An attribute can be removed by using a value of `nothing`. Example:

```julia
System1 = Model(v = input | Var(init = 1.0), ...)

# System2 = Model(v = input, ...)
System2 = System1 | Map(v = Var(init = nothing), ...)
```

The following attributes are also defined for constructor `Var`,
but have **no effect yet**.
Using `min, max, info` already now, might be useful for model libraries:

| Var Key   | Shortcut          | Shortcut value        |  Description                     |
|:--------- |:------------------|:----------------------|:---------------------------------|
| constant  | constant          | Var(constant = true)  | If true, value cannot be changed |
| min, max  | interval(a,b)     | Var(min = a, max = b) | Allowed variable value range     |
| info      | info"..."         | Var(info="...")       | Description                      |

Example:
```julia
v = output | interval(0.0,1.0) | Var(start = zeros(3)u"N*m") | info"An output variable"

# Same as: v = Var(output = true, min = 0.0, max = 1.0,
#                  start = zeros(3)u"N*m", info = "An output variable")
```


## A.2 Dictionaries and quoted expressions

The fundamental mechanism for defining models, variables and parameter modifications in Modia are ordered dictionaries, i.e. a list of key/value pairs:

```julia
julia> using OrderedCollections

julia> S = OrderedDict(:p=>5, :q=>10)
OrderedDict{Symbol, Int64} with 2 entries:
  :p => 5
  :q => 10
```

It is also possible to define a constructor `Model` with keyword arguments which creates the ordered dictionary:

```julia
julia> Model(; kwargs...) = OrderedDict{Symbol, Any}(kwargs)
Model (generic function with 1 method)

julia> T=Model(q=100, r=200)
OrderedDict{Symbol, Any} with 2 entries:
  :q => 100
  :r => 200
```

The values can also be a quoted expression, i.e. an expression enclosed in `:( )`, an array of quoted expressions enclosed in `:[ ]` or just a quoted symbol, `:x`.
This mechanism is used to encode equations and expressions of the model which needs to be manipulated before the model can be simulated.

Julia defines a very useful merge operation between dictionaries:

```julia
julia> merge(S, T)
OrderedDict{Symbol, Any} with 3 entries:
  :p => 5
  :q => 100
  :r => 200
```

If a key already exists in the first dictionary (like `:q`), its value is overwritten (like `:r`) otherwise it's added (like `:p`).
Such a merge semantic allows for unification of parameter modifications and inheritance as will be demonstrated below.

## A.3 MergeModels algorithm

The basics of the `mergeModels` algorithm and the merge operator `|` are defined as follows (without logging):

```julia
function mergeModels(m1::AbstractDict, m2::AbstractDict, env=Symbol())
    result = deepcopy(m1)
    for (k,v) in m2)
        if typeof(v) <: AbstractDict
            if k in keys(result) && ! (:_redeclare in keys(v))
                if typeof(result[k]) <: AbstractDict
                    result[k] = mergeModels(result[k], v, k)
                end
            else
                result[k] = v
            end
        elseif v === nothing
            delete!(result, k)
        elseif k in keys(result) && k == :equations
            equa = copy(result[k])
            push!(equa.args, v.args...)
            result[k] = equa
        else
            result[k] = v
        end
    end
    return result
end

|(m::AbstractDict, n::AbstractDict) =  mergeModels(m, n)

```

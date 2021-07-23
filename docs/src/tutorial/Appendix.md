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


## A.2 Named tuples and quoted expressions

The fundamental mechanism for defining models in TinyModia are named tuples which is a list of key/value pairs enclosed in parentheses:

```julia
julia> S=(p=5, q=10)
(p = 5, q = 10)

julia> typeof(S)
NamedTuple{(:p, :q),Tuple{Int64,Int64}}
```

Named tuples are conceptually similar to dictionaries (Dict), but the constructor syntax is simpler. Note that if only one key/value pair is given, a comma must preceed the final parentheses: `(p = 5, )`. It is also possible to define named tuples using a keyword argument list, i.e. a list starting with a semi-colon: `z=(;p=5)`.

The values can also be a quoted expression, i.e. an expression enclosed in `:( )`, an array of quoted expressions encloded in `:[ ]` or just a quoted symbol, `:x`. This mechanism is used to encode equations and expressions of the model which needs to be manipulated before the model can be simulated.

Julia defines a very useful merge operation between named tuples (and dictionaries):

```julia
julia> T=(q=100, r=200)
(q = 100, r = 200)

julia> merge(S, T)
(p = 5, q = 100, r = 200)
```

If a key already exists `q` in the first named tuple, it's value is overwritten otherwise it's added, `r`. Such a merge semantics allows for unification of parameter modifications and inheritance as will be demonstrated below.

## A.3 MergeModels algorithm

The `mergeModels` algorithm is defined as follows (without logging):

```julia
function mergeModels(m1::NamedTuple, m2::NamedTuple, env=Symbol())
    mergedModels = OrderedDict{Symbol,Any}(pairs(m1)) # Convert the named tuple m1 to an OrderedDict
    for (k,v) in collect(pairs(m2))
        if typeof(v) <: NamedTuple
            if k in keys(mergedModels) && ! (:_redeclare in keys(v))
                mergedModels[k] = mergeModels(mergedModels[k], v, k)
            else
                mergedModels[k] = v
            end
        elseif v === nothing
            delete!(mergedModels, k)
        else
            mergedModels[k] = v
        end
    end
    return (; mergedModels...) # Transform OrderedDict to named tuple
end

|(m::NamedTuple, n::NamedTuple) =  mergeModels(m, n)

Redeclare = ( _redeclare = true, )
```

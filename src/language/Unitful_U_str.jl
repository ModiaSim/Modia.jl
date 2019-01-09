"""
    module Unitful_U_str

Temporary solution to fix issue with Unitful.
For details see: https://github.com/ModiaSim/Modia.jl/issues/49

This module exports string macro @U_str. It is then possible to use:

```julia
using Unitful
using Modia    # Unitful_U_str.jl is included in Modia.jl

const T1_type = u"N*m"    # gives warning
const T2_type = u"N*m"    # gives no warning
```
"""
module Unitful_U_str

export @U_str


using Unitful


macro U_str(unit)
    ex = Meta.parse(unit)
    esc(replace_value(ex))
end

const allowed_funcs = [:*, :/, :^, :sqrt, :√, :+, :-, ://]
function replace_value(ex::Expr)
    if ex.head == :call
        ex.args[1] in allowed_funcs ||
            error("""$(ex.args[1]) is not a valid function call when parsing a unit.
             Only the following functions are allowed: $allowed_funcs""")
        for i=2:length(ex.args)
            if typeof(ex.args[i])==Symbol || typeof(ex.args[i])==Expr
                ex.args[i]=replace_value(ex.args[i])
            end
        end
        return ex
    elseif ex.head == :tuple
        for i=1:length(ex.args)
            if typeof(ex.args[i])==Symbol
                ex.args[i]=replace_value(ex.args[i])
            else
                error("only use symbols inside the tuple.")
            end
        end
        return ex
    else
        error("Expr head $(ex.head) must equal :call or :tuple")
    end
end

replace_value(sym::Symbol) = Unitful.replace_value(sym)

replace_value(literal::Number) = literal

end
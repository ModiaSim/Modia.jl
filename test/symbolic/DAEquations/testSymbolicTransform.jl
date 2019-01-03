"""
Module with tests of SymbolicTransform.

* Author: Hilding Elmqvist, Mogram AB  
* Date: July-August 2016
* License: MIT

"""
module testSymbolicTransform

using Modia.SymbolicTransform
using Modia.Utilities

# Desired:
#   using Test
#   using LinearAlgebra
#
# In order that these packages need not to be defined in the user environment, they are included via Modia:
import Modia

@static if VERSION < v"0.7.0-DEV.2005"
  using Base.Test
else
  using Modia.Test
end

@static if VERSION < v"0.7.0-DEV.2005"
else
  using Modia.LinearAlgebra
end

using Base.Meta: isexpr

# Copied from Instantiation.jl
prettyfy(ex) = ex
#prettyfy(der::Der) = Symbol("der("*string(der.base.name)*")")
##prettyfy(der::Der) = Symbol("der("*replace(string(der.base.name), ".", "_")*")")
#prettyfy(get::GetField) = get.name
##prettyfy(get::GetField) = Symbol(replace(string(get.name), ".", "_")) # get.name # Handle dummy derivatives
#prettyfy(s::Symbol) = s
function prettyfy(ex::Expr)
  if isexpr(ex, :quote) || isexpr(ex, :line)
    nothing
  elseif isexpr(ex, :block)
    prettyfy(ex.args[2])
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
function prettyPrint(e::Expr)
    ex = prettyfy(e)
    if ex.head === :quote
      return ex
    elseif ex.head === :(:=)
      return string(prettyPrint(ex.args[1]), " := ", prettyPrint(ex.args[2]))
    end
    Expr(ex.head, [prettyPrint(arg) for arg in ex.args]...)
end

function showSolve(e, x)
  println("\nSolve: ", x, " from: ", prettyPrint(e))
  sol, solved = solve(e, x)
  if ! solved 
    println("NOT SOLVED")
  end
  println(prettyPrint(sol))
  return string(prettyPrint(sol))
end

function showDifferentiate(e)
  println("\nEquation: ", prettyPrint(e))
  der = differentiate(e)
  println("Differentiated: ", prettyPrint(der))
  return string(prettyPrint(der))
end

function testSolve()
  println("\nTest solve")

  sol = showSolve(:(y = x), :x)  
  @test sol == "x = y"

  sol = showSolve(:(y = x + z), :x)
  @test sol == "x = y - z"
  
  sol = showSolve(Expr(:(=), :y, Expr(:call, :+, :x, :z, :v, :w)), :x)
  @test prettyPrint(sol) == "x = y - (z + v + w)"

  sol = showSolve(Expr(:(=), :y, Expr(:call, :+, :x, :z, :v, :w)), :z)
  @test sol == "z = (y - x) - (v + w)"

  sol = showSolve(Expr(:(=), :y, Expr(:call, :+, :x, :z, :v, :w)), :v)
  @test sol == "v = ((y - x) - z) - w"

  sol = showSolve(Expr(:(=), :y, Expr(:call, :+, :x, :z, :v, :w)), :w)
  @test sol == "w = ((y - x) - z) - v"

  sol = showSolve(:(y = x - z), :x)
  @test sol == "x = y + z"

  sol = showSolve(:(y = x - z - w), :x)
  @test sol == "x = (y + w) + z"
    
  sol = showSolve(Expr(:(=), :y, Expr(:call, :-, :x, :z, :v, :w)), :x)
  @test sol == "x = y + (z + v + w)"

  sol = showSolve(Expr(:(=), :y, Expr(:call, :-, :x, :z, :v, :w)), :v) # Solve: v from: y = x - z - v - w
  @test sol == "v = ((x - y) - z) - w"

  sol = showSolve(:(y = z - x), :x)
  @test sol == "x = z - y"

  sol = showSolve(:(y = x*z), :x)
  @test sol == "x = y / z"
  
  sol = showSolve(:(y = x*z*z*z), :x)
  @test sol == "x = y / (z * z * z)"

  sol = showSolve(Expr(:(=), :y, Expr(:call, :/, :x, :z, :w)), :x)
  @test sol == "x = y * (z * w)"

  sol = showSolve(Expr(:(=), :y, Expr(:call, /, :x, :z, :w)), :z)
  @test sol == "z = (x / y) / w"

  sol = showSolve(:(y = x / z), :x)
  @test sol == "x = y * z"

  sol = showSolve(:(y = x / z), :z)
  @test sol == "z = x / y"

  sol = showSolve(:(y = x \ z), :x)
  @test sol == "x \\ z = y" 
  println("\n\n----------------------\n")

end

function testDifferentiate()
  println("\nTest differentiate")
    
  der = showDifferentiate(:(x + 5 + z = w))
  @test der == "der(x) + der(z) = der(w)"
     
  der = showDifferentiate(differentiate(:(x + 5 + z = w)))
  @test der == "der(der(x)) + der(der(z)) = der(der(w))"
  
  der = showDifferentiate(Expr(:(=), Expr(:call, :+, :x), :w))
  @test der == "der(x) = der(w)"
  
  der = showDifferentiate(:(2 + 3 = w))
  @test der == "0.0 = der(w)"
  
  der = showDifferentiate(Expr(:(=), Expr(:call, :-, :x), :w))
  @test der == "-(der(x)) = der(w)"
  
  der = showDifferentiate(:(x - 5 - z = w))
  @test der == "der(x) - der(z) = der(w)"
  
  der = showDifferentiate(:(5x = w))
  @test der == "5 * der(x) = der(w)"
    
  der = showDifferentiate(:(x * 5 * z = w))
  @test der == "der(x) * 5 * z + x * 5 * der(z) = der(w)"
  
  der = showDifferentiate(:(4 * 5 * 6 = w))
  @test der == "0.0 = der(w)"
  
  der = showDifferentiate(:(y = x/y))  
  @test der == "der(y) = der(x) / y + (x / y ^ 2) * der(y)"
  
  der = showDifferentiate(:(y = x/5))  
  @test der == "der(y) = der(x) / 5"
  
  der = showDifferentiate(:(y = 5/y))  
  @test der == "der(y) = (5 / y ^ 2) * der(y)"

  der = showDifferentiate(:(y = [1, x]))
  @test der == "der(y) = [0.0, der(x)]"
  
  der = showDifferentiate(:(y = [2x 3x; 4x 5x]))
  @test der == "der(y) = [2 * der(x) 3 * der(x); 4 * der(x) 5 * der(x)]"
  
  der = showDifferentiate(:(y = [2*x 3x; 4x 5x]*[1, x]))
  @test der == "der(y) = [2 * der(x) 3 * der(x); 4 * der(x) 5 * der(x)] * [1, x] + [2x 3x; 4x 5x] * [0.0, der(x)]"
  
  der = showDifferentiate(:(y = transpose(B) + B´))
  @test der == "der(y) = transpose(der(B)) + der(B´)"
   
  der = showDifferentiate(:(y = x[5, 6]))  
  @test der == "der(y) = (der(x))[5, 6]"
  
  der = showDifferentiate(:(y = x[5:7]))  
  @test der == "der(y) = (der(x))[5:7]"
  
#=    der = showDifferentiate(:(y = [x for x in z]))  
  @test der == "der(y) = [x for x = der(z)]"

  der = showDifferentiate(:(y = [x[i] for i in 1:5]))  
  @test der == "der(y) = [x[i] for i = nothing]"
=#      
  der = showDifferentiate(:(y = sin(x)))
  @test der == "der(y) = cos(x) * der(x)"
  
  der = showDifferentiate(:(y = cos(x)))
  @test der == "der(y) = -(sin(x)) * der(x)"
  
  der = showDifferentiate(:(y = tan(x)))
  @test der == "der(y) = (1 / cos(x) ^ 2) * der(x)"
  
  der = showDifferentiate(:(y = exp(x)))
  @test der == "der(y) = exp(x) * der(x)"

  der = showDifferentiate(:(y = x^y))
  @test der == "der(y) = y * x ^ (y - 1) * der(x) + x ^ y * log(x) * der(y)"
  
  der = showDifferentiate(:(y = log(x)))
  @test der == "der(y) = (1 / x) * der(x)"

  der = showDifferentiate(:(y = asin(x)))  
  @test der == "der(y) = (1 / sqrt(1 - x ^ 2)) * der(x)"

  der = showDifferentiate(:(y = acos(x)))  
  @test der == "der(y) = (-1 / sqrt(1 - x ^ 2)) * der(x)"

  der = showDifferentiate(:(y = atan(x)))  
  @test der == "der(y) = (1 / (1 + x ^ 2)) * der(x)"
  
  der = showDifferentiate(:(y = f(x, 5, z)))
  @test der == "der(y) = f_der_1(x, 5, z) * der(x) + f_der_3(x, 5, z) * der(z)"
  
  der = showDifferentiate(:(y = f(x, 5, g(z))))
  @test der == "der(y) = f_der_1(x, 5, g(z)) * der(x) + f_der_3(x, 5, g(z)) * (g_der_1(z) * der(z))"
  
  der = showDifferentiate(:(y = true ? x : y))  
@static if VERSION < v"0.7.0-DEV.2005"
    @test der == "der(y) = if true\n        der(x)\n    else \n        der(y)\n    end"
  else
    @test der == "der(y) = if true\n        der(x)\n    else\n        der(y)\n    end"
  end
#=
  der = showDifferentiate(:(y = if b; x elseif false; y else z end))  
  @test der == "der(y) = if b\n        der(x)\n    else\n        if false\n            der(y)\n        else \n            der(z)\n        end\n    end"
=#  
  der = showDifferentiate(:(y = time))  
  @test der == "der(y) = 1.0"
  
  setTimeInvariants([:a])
  der = showDifferentiate(:(y = a*x))  
  @test der == "der(y) = a * der(x)"
  
  println("\n\n----------------------\n")
  
  
  
  # check: time
end

@testset "Symbolic" begin

  @testset "Solve" begin
    testSolve()
  end # testset 

  @testset "Differentiate" begin
    testDifferentiate()
  end

end

end

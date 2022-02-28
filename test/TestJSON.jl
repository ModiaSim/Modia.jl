module TestJSON

using  Modia
import Modia.JSON

mutable struct S
	a
	b
end

mutable struct T
	a
	b
end

s = S(1,T(10, 25))

println("JSON.print(s):")
JSON.print(s)
println()

n = (a=1, b=(a=10, b=25))

println("JSON.print(n):")
JSON.print(n)
println()

println("JSON.print(s) == JSON.print(n):")
println(JSON.print(s) == JSON.print(n))

e = :(x = sin(time)*(a+1))

println()
JSON.print(e)
println()

println()
ej = JSON.json(e)
@show ej

println()
ee = JSON.parse(ej)
@show eval(ee)

end
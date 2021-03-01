module TestJSON
import JSON

mutable struct S
	a
	b
end

mutable struct T
	a
	b
end

s = S(1,T(10, 25))

JSON.print(s)

n = (a=1, b=(a=10, b=25))

JSON.print(n)

@show JSON.print(s) == JSON.print(n)


end
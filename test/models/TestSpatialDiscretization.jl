module TestSpatialDiscretization

using Modia
using ModiaMath.plot

const n = 5
@model SpatialDiscretization begin
    u = Var(start=zeros(n))
    c = 1
    @equations begin
        # der(u[2:n]) = ( u[2:n] - u[1:n-1] ) / c
        der(u) = [1; ( u[2:n] - u[1:n - 1] ) / c]
    end
end

result = simulate(SpatialDiscretization, 1)
plot(result, "u")

end


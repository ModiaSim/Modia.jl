module TestEquations

using Modia

# Desired:
#   using ModiaMath: plot
#   using Test
#   using LinearAlgebra
#
# In order that these packages need not to be defined in the user environment, they are included via Modia:
using Modia.ModiaMath: plot

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Modia.Test
    using Modia.LinearAlgebra
    eye(n) = Matrix{Float64}(I, n, n)
    diagm(v) = Matrix(Diagonal(v))
    q(qr) = qr.Q
    r(qr) = qr.R
end

# @testset "Equations" begin

@model Test begin
    a = -1; b = 1; c = 1
    n = 10; 
    M = Float()
    u = Float(); x = Float(start=0.0)
    on = Variable(); s = Variable()
    QR = Float(); Q = Float(); R = Float()
    dummy = Variable(size=())
    @equations begin
        u = if s == "on"; 10 else M end # Varying size of u
        der(x) - a * x = b * (s == "on" ? u : 10 * u[1,1])
        on = time > 1 && time < 3
        s := on ? "on" : "off"
        dummy = if time < 0.00001 || abs(time - 2.0) < 0.01; println(" Time=", time, ":  size of u: ", size(u)) else nothing end
        ((9 + 1) * diagm(1:n) + ones(n, n))^3 * M * diagm(1:n) = eye(n)
#        QR = qr(M)
#       Q = QR.Q; R = QR.R # ERROR: LoadError: No substitution provided for: this.QR.Q
#       Q = q(QR); R = r(QR) # ERROR: LoadError: MethodError: no method matching length(::LinearAlgebra.QRCompactWY{Float64,Array{Float64,2}})
        QR = Tuple{Array,Array}(qr(M)) # Tuple{Array,Array}() needed for Julia 0.7
        Q = QR[1]; R = QR[2]
    end
end;

result = simulate(Test, 5, storeEliminated=false, logTranslation=true, removeSingularities=false)
plot(result, ("x"), heading="TestEquations", figure=2)
@test isapprox(result["x"][end], 1.1776785083961023; atol=1e-8)


end

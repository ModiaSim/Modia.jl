
println("\nRunLargeTest: Demonstrates arrays of components and timing.")

using Modia
using Modia.Electric

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

@testset "RunLargeTests" begin

    @model LPfilter begin
        R = Resistor(R=100.0)
        C = Capacitor(C=2.5E-3, v=Float(start=0.0))
        V = ConstantVoltage(V=10.0)
        ground = Ground()
        @equations begin
            connect(V.n, ground.p)
            connect(V.p, R.p)
            connect(R.n, C.p)
            connect(C.n, V.n)
        end
    end 

    nFilters = 100
    @model ManyFilters begin
        F = [LPfilter() for i in 1:nFilters]
    end

    @time simulate(ManyFilters, 2, logTiming=true, storeEliminated=false, removeSingularities=false)

end # testset
nothing

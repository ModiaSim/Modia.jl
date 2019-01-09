module TestVariableTypes

println("\nTestVariableTypes: Demonstrating the handling of various variable types")

# Desired:
#   using Test
#   using LinearAlgebray
#   using ModiaMath.plot
#   using StaticArrays
#   using Unitful
#
# In order that these packages need not to be defined in the user environment, they are included via Modia:
using Modia

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
    identity(m,n) = eye(m,n)
else
    using Modia.Test
    using Modia.LinearAlgebra
    identity(m, n) = Matrix{Float64}(I, m, n)
end

using  Modia.ModiaMath: plot

using  Modia.StaticArrays
using  Modia.Unitful
using  Modia.Unitful.DefaultSymbols
import Modia.Unitful:
    nm, μm, mm, cm, m, km, inch, ft, mi,
    ac,
    mg, g, kg, A,
    Ra, °F, °C, K,
    rad, °,
    ms, s, minute, hr,
    J, A, N, mol, cd, V,
    mW, W, Hz


#const Vec3 = Vec{3,Float64}
const Vec3 = SVector{3,Float64}

@testset "VariableTypes" begin

    # Var(; args...) = Variable(; args...)
    # Float(; args...) = Variable(T=Float64; args...)

    FloatScalar(; args...) = Variable(T=Float64, size=(); args...)
    FloatArray(; args...) = Variable(T=Array{Float64,1}; args...)
    FloatMatrix(; args...) = Variable(T=Array{Float64,2}; args...)

    Float3(;args...) = Variable(T=Float64, size=(3,); args...)
    Float3x3(;args...) = Variable(T=Float64, size=(3, 3); args...)

    FixedFloat3(;args...) = Variable(T=Vec{3,Float64}; args...)
    FixedFloat3x3(;args...) = Variable(T=Mat{3,3,Float64}; args...)

    Voltage(; args...) = Variable(T=Unitful.V; args...)

    # Real = (;args...)->Variable(typ=Float64; args...)

    @model TestVariableTypes1 begin
        f = Float()
        x = Float(start=3.0E10, nominal=20.0E10)
        b = Boolean()
        i = Integ()
        s = Str()
        c = Var(T=Complex)
        re = Float()
        im = Float()
        c1 = Var(start=Complex(2.0, 0.0))
        # c2 = Var(start=Complex(2.0, 3.0))  # InexactError
        @equations begin
            f = 1
            der(x) = -x
            b = true
            i = 1
            s = "asdf"
            c = Complex(2.0, 3.0)
            re = real(c)
            im = imag(c)
            der(c1) = Complex(1.0, 0.0) # nominal no supported
            # der(c2) = Complex(1.5, 2.5)  # InexactError # Allocation of Complex in state vector not supported yet.
        end
    end 

    result = checkSimulation(TestVariableTypes1, 3, "i", 1)
    @show result["f"][end]
    @show result["b"][end]
    @show result["i"][end]
    @show result["s"][end]
    @show result["c"][end]
    @show result["re"][end]
    @show result["im"][end]
    # @show result["c2"][end]
    @test result["s"][end] == "asdf"

    # ----------------------

    @model TestArrays1 begin
        f = Float(start=[1,2,3], nominal=[1,2,3])
        b = Boolean()
        i = Integ()
        s = Str()
        c1 = Var(start=[Complex(1.0, 0.0), Complex(2.0)])

        @equations begin
            der(f) = [2,4,6]
            b = [false,true]
            i = [1,2]
            s = ["asdf", "qwerty"]
            der(c1) = [Complex(2.0, 0.0), Complex(4.0, 0.0)]
            # der(c2) = Complex(1.0, 2.0)  # InexactError
        end
    end 

    result = checkSimulation(TestArrays1, 1, "i", [1,2])
    # result = checkSimulation(TestArrays1, 1, "i", [1,2], expandArrayIncidence=true)
    if result != nothing
        @show keys(result)
        @show result["f"][end,:]
        @show result["der(f)"][end,:]
        @show result["b"][end]
        @show result["i"][end]
        @show result["s"][end]
        @show result["c1"][end,:]
        @show result["der(c1)"][end,:]
    end
    @test result["i"][end] == [1,2]
    println()

# ----------------------

    @model TestVariableTypes2 begin
        scalar = FloatScalar()
        vector = FloatArray()
        matrix = FloatMatrix()
        anyVector = Variable(T=Array{Float64,1})
        v1 = Variable(T=SVector{3,Float64}, start=ones(3), min=3)
        v2 = Variable(T=Vec3, start=ones(3), min=3)
        v3 = Variable(T=Vec3)

        v4 = Float3(min=3, nominal=[1,2,3])
        v5 = Var()
        v6 = Voltage(start=0.0)
        @equations begin
            scalar = 1
            vector = [1.0, 2.0, 3.0]
            matrix = identity(3, 3)
            anyVector = [1,2,3]
            v1 = ones(3)
            v2 = ones(3)
            v3 = 30.0

            v4 = [1,2,3]
            v5 = [1,2,3]
            #v6 = 10.0u"V"
            v6 = 10.0
        end
    end 

    result = simulate(TestVariableTypes2, 1, storeEliminated=false)

    println("Variable(T=Array{Float64,1}; args...) does not work with storeEliminated=true!")

# ----------------------

#=
    @model TestVariableUnits1 begin
        # Real3 := Variable{Float64,3}
        # typealias Voltage Variable{Volt}
        v1 = Variable()
        v2 = Var() # Voltage()  # Voltage should work!
        v3 = Var(T=typeof(1.0u"kg"), start=1.0u"kg")
        v4 = Var()
        m = 3.1u"kg" + 2.5kg + 3g
        F = 2.5u"kg*m/s^2" + 2u"N"
        a = Var() # (T=typeof(1.0u"m/s^2")) 
        l = 2u"m" + 30cm + 55.7mm  # + 2m not possible due to variable m
        @equations begin
            # SIunits
            v1 = 10.0u"m"
            v2 = 5.5u"V"
            # Unitful
            v3 = 1u"kg" + 2u"g" + 10 * sin(2 * time) * u"hg"  
            #  der(v3) = 1.0u"kg/s"  # Differential equations with units are not handled yet
            v4 = 3u"kg" * 0.5u"m/s^2" + time * u"N"
            a = F / m
        end
    end 

    result = simulate(TestVariableUnits1, 1)
    #result = checkSimulation(TestVariableUnits1, 1, "v1", 10u"m")
    
    @show result["v2"][end]
    @show result["v3"][end]
    @show result["v4"][end]
    @show result["m"][end]
    @show result["F"][end]
    @show result["a"][end]
    @test result["v3"][end] == 1.9112974268256817u"kg"
    @test result["v4"][end] == 2.5u"kg*m/s^2" == 2.5N
    @test result["a"][end] == 0.8031411743708727u"m/s^2"  # Unit deduction
    @show result["l"][end]

    plot(result, ("v3", "v4", "a"), figure=1)

=#

end

end

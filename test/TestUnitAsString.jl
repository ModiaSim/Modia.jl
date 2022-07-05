module TestUnitAsString

using  Modia
using  Test
using  Modia.Measurements
import Modia.MonteCarloMeasurements
import Modia.MonteCarloMeasurements.StaticParticles

v1 = 2.0u"m/s"
unit_v1 = unit(v1)
v1_unit = unitAsParseableString( unit_v1  )   # = "m*s^-1"
v2_withoutUnit = 2.0
code = :( $v2_withoutUnit*@u_str($v1_unit) )  # = 2.0u"m*s^-1"
v2 = eval(code)
#@show v1
#@show v1_unit
#@show v2

@test v1==v2


mutable struct Data{FloatType <: AbstractFloat}
    velocity::Modia.quantity(FloatType, u"m/s")
end

data = Data{Float64}(2.0u"mm/s")
v3   = data.velocity
#@show v3   # v3 = 0.002 m s^-1

@test unitAsParseableString( unit(v3) ) == "m*s^-1"


v = 2.0
@test typeof(v) <: AbstractFloat
@test Modia.baseType(typeof(v)) == Float64
@test Modia.isQuantity(             typeof(v)) == false
@test Modia.isMeasurements(          typeof(v)) == false
@test Modia.isMonteCarloMeasurements(typeof(v)) == false
Base.floatmax(typeof(v))

v = 2.0u"m/s"
@test typeof(v) <: Unitful.Quantity
@test Modia.baseType(typeof(v)) == Float64
@test Modia.isQuantity(             typeof(v)) == true
@test Modia.isMeasurements(          typeof(v)) == false
@test Modia.isMonteCarloMeasurements(typeof(v)) == false
Base.floatmax(typeof(v))


### Measurements
v = 2.0 ± 0.1
@test typeof(v) <: Measurements.Measurement
@test Modia.baseType(typeof(v)) == Float64
@test Modia.isQuantity(              typeof(v)) == false
@test Modia.isMeasurements(          typeof(v)) == true
@test Modia.isMonteCarloMeasurements(typeof(v)) == false
Base.floatmax(typeof(v))

v = (2.0 ± 0.1)u"m/s"
@test typeof(v) <: Unitful.Quantity
@test Modia.baseType(typeof(v)) <: Measurements.Measurement
@test Modia.baseType(Modia.baseType(typeof(v))) == Float64
@test Modia.isQuantity(              typeof(v)) == true
@test Modia.isMeasurements(          typeof(v)) == true
@test Modia.isMonteCarloMeasurements(typeof(v)) == false
Base.floatmax(typeof(v))


### MonteCarloMeasurements
const nparticles = 100
normal(mue,sigma) = StaticParticles(nparticles,MonteCarloMeasurements.Distributions.Normal(mue,sigma))

v = normal(2.0, 0.1) 
@test typeof(v) <: Modia.MonteCarloMeasurements.AbstractParticles
@test Modia.baseType(typeof(v)) == Float64
@test Modia.isQuantity(              typeof(v)) == false
@test Modia.isMeasurements(          typeof(v)) == false
@test Modia.isMonteCarloMeasurements(typeof(v)) == true
Base.floatmax(typeof(v))

v = normal(2.0, 0.1)u"m/s"
@test typeof(v) <: Modia.MonteCarloMeasurements.AbstractParticles
@test Modia.baseType(typeof(v)) <: Unitful.Quantity
@test Modia.baseType(Modia.baseType(typeof(v))) == Float64
@test Modia.isQuantity(              typeof(v)) == true
@test Modia.isMeasurements(          typeof(v)) == false
@test Modia.isMonteCarloMeasurements(typeof(v)) == true
Base.floatmax(typeof(v))


logCode = false

FirstOrder1 = Model(
    T = 0.2u"s",
    k = 2.0,
    x = Var(init=0.3),
    equations = :[T * der(x) + x = k]
)
FirstOrder2 = FirstOrder1 | Map(k = 2.0u"m/s", x = Var(init = 0.3u"m/s"))
FirstOrder3 = FirstOrder1 | Map(k = 2.0 ± 0.1, x = Var(init = 0.3 ± 0.1))
FirstOrder4 = FirstOrder1 | Map(k = (2.0 ± 0.1)u"m/s" , x = Var(init = (0.3 ± 0.1)u"m/s"))
FirstOrder5 = FirstOrder1 | Map(k = normal(2.0, 0.1/3), x = Var(init = normal(0.3, 0.1/3)))
FirstOrder6 = FirstOrder1 | Map(k = normal(2.0, 0.1/3)u"m/s", x = Var(init = normal(0.3, 0.1/3)u"m/s"))

if false
    # Remove tests, to drastically reduce the time for the test
    firstOrder1 = @instantiateModel(FirstOrder1, logCode=logCode)
    firstOrder2 = @instantiateModel(FirstOrder2, logCode=logCode)
    firstOrder3 = @instantiateModel(FirstOrder3, logCode=logCode, FloatType=Measurements.Measurement{Float64})
    firstOrder4 = @instantiateModel(FirstOrder4, logCode=logCode, FloatType=Measurements.Measurement{Float64})
    firstOrder5 = @instantiateModel(FirstOrder5, logCode=logCode, FloatType=StaticParticles{Float64,nparticles})
    firstOrder6 = @instantiateModel(FirstOrder6, logCode=logCode, FloatType=StaticParticles{Float64,nparticles})
    
    simulate!(firstOrder1)
    simulate!(firstOrder2)
    simulate!(firstOrder3)
    simulate!(firstOrder4)
    simulate!(firstOrder5)
    simulate!(firstOrder6)
end

end
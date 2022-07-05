module Tutorial

# 1 Getting started

using Modia

# Define model
SimpleModel = Model(
    T = 0.4,
    x = Var(init=0.2),
    equations = :[T * der(x) + x = 1],
)

# Transform to ODE form
simpleModel = @instantiateModel(SimpleModel)

# Simulate with a default integrator of DifferentialEquations
simulate!(simpleModel, stopTime = 1.2)

# Simulate with a specific integrator (Tsit5) and use a unit for stopTime
simulate!(simpleModel, Tsit5(), stopTime = 1.2u"s")

# Produce a line plot
@usingPlotPackage
plot(simpleModel, ("x", "der(x)"), figure=1)



# 2.1 - Equation oriented models

using Modia

LowPassFilter = Model(
    T = 0.2,
    u = input,
    y = output | Var(:x),
    x = Var(init=0),
    equations = :[T * der(x) + x = u],
)
    

   
# 2.2 Merging models  

HighPassFilter = LowPassFilter | Model( y = :(-x + u) )

setLogMerge(true)

LowAndHighPassFilter = LowPassFilter | Model(
        y = nothing,
        low = output | Var(:x),
        high = output | Var(:(-x + u)),
    )
    
setLogMerge(false)

@showModel LowAndHighPassFilter



# 2.3 Functions and tables

TestLowAndHighPassFilter = LowAndHighPassFilter | Model(
        u = :(sin( (time+1u"s")*u"1/s/s" * time)*u"V"),
        x = Var(init=0.2u"V")
    )
    
    
using Interpolations

table = CubicSplineInterpolation(0:0.5:2.0, [0.0, 0.7, 2.0, 1.8, 1.2])
TestLowAndHighPassFilter2 = TestLowAndHighPassFilter | Map(u = :(table(time*u"1/s")*u"V")) 
    
    
function ref(time)
    y1 = sin(time)
    y2 = cos(time)
    return (y1,y2)
end

using Modia
TestMultiReturningFunction1 = Model(
    equations = :[
        (y1,y2) = ref(time)
        y3 = y1+y2
    ]
)
testMultiReturningFunction1 = @instantiateModel(TestMultiReturningFunction1, unitless=true)
simulate!(testMultiReturningFunction1, stopTime=10.0)
plot(testMultiReturningFunction1, ["y1", "y2", "y3"], figure=2)


# 2.4 Hierarchical modeling

TwoFilters = (
    high = HighPassFilter,
    low = LowPassFilter,
)

BandPassFilter = (
    u = input,
    y = output,
    high = HighPassFilter | Map(T=0.5, x=Var(init=0.1u"V")),
    low = LowPassFilter | Map(x=Var(init=0.2u"V")),
    equations = :[
        high.u = u,
        low.u = high.y,
        y = low.y]
)

TestBandPassFilter = BandPassFilter | Map(
        u = :(sin( (time+1u"s")*u"1/s/s" * time)*u"V")
    )

#bandPassFilter = @instantiateModel(TestBandPassFilter, logStateSelection=true)
#simulate!(bandPassFilter, Tsit5(), stopTime = 50u"s")
#plot(bandPassFilter, ["u", "y"], figure=3)


StateSpace = Model(
    A = fill(0.0, 0, 0),
    B = fill(0.0, 0, 0),
    C = fill(0.0, 0, 0),
    D = fill(0.0, 0, 0),
    u = input,
    y = output,
    x = Var(init = zeros(0)),
    equations = :[
        der(x) = A*x + B*u
             y = C*x + D*u
    ]
)

SecondOrder = Model(
    w = 20.0,
    D =  0.1,
    k =  2.0,
    sys = StateSpace | Map(A = :([  0        1;
                                 -w^2  -2*D*w]),
                           B = :([0; w^2;;]),    # Julia 1.7: Trailing ";" defines a column matrix
                           C = :([k 0]),
                           D = :(zeros(1,1)),
                           x = Var(init = zeros(2)) ),
    equations = :[sys.u = [1.0]]
)
secondOrder = @instantiateModel(SecondOrder)
simulate!(secondOrder, stopTime=2.0)
plot(secondOrder, ("sys.u", "sys.x", "sys.y"), figure=3)


using Modia.StaticArrays
TestArray1 = Model(
    v = Var(init=SVector{3,Float64}(1.0, 2.0, 3.0)),
    equations = :[der(v) = -v]
)
testArray1 = @instantiateModel(TestArray1, logCode=true)
simulate!(testArray1, stopTime=2.0)
plot(testArray1, "v", figure=4)

TestArray2 = Model(
    v = Var(init=[1.0, 2.0, 3.0]),
    equations = :[der(v) = -v]
)
testArray2 = @instantiateModel(TestArray2)
simulate!(testArray2, stopTime=2.0, merge=Map(v = [4.0, 3.0, 2.0, 1.0]))
plot(testArray2, "v", figure=5)

TwoInertiasAndIdealGearTooManyInits = Model(
    J1    = 50.0,
    J2    = 100.0,
    ratio = 2.0,
    f     = 3.0, # Hz

    phi1 = Var(init = 0.0), # Absolute angle of inertia1
    w1   = Var(init = 0.0), # Absolute angular velocity of inertia1
    phi2 = Var(init = 0.0), # Absolute angle of inertia2
    w2   = Var(init = 0.0), # Absolute angular velocity of inertia2

    equations = :[
        tau = 2.0*sin(2*3.14*f*time/u"s")

        # inertia1
        w1 = der(phi1)
        J1*der(w1) = tau - tau1

        # ideal gear
        phi1 = ratio*phi2
        ratio*tau1 = tau2

        # inertia2
        w2 = der(phi2)
        J2*der(w2) = tau2
    ]
)

drive1 = @instantiateModel(TwoInertiasAndIdealGearTooManyInits, logStateSelection=true)

end
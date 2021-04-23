module Tutorial

using Modia, ModiaPlot

# 1 Getting started

using Modia, ModiaPlot

# Define model
SimpleModel = Model(
    T = 0.4,
    x = Var(init=0.2),
    equation = :[T * der(x) + x = 1],
)

# Transform to ODE form
simpleModel = @instantiateModel(SimpleModel)

# Simulate with a default integrator of DifferentialEquations
simulate!(simpleModel, stopTime = 1.2)

# Simulate with a specific integrator (Tsit5) and use a unit for stopTime
simulate!(simpleModel, Tsit5(), stopTime = 1.2u"s")

# Produce a line plot with GLMakie
plot(simpleModel, ("x", "der(x)"), figure=1)



# 2.1 - Equation oriented models

using Modia

LowPassFilter = Model(
    T = 0.2,
    u = input,
    y = output | Var(:x),
    x = Var(init=0),
    equation = :[T * der(x) + x = u],
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

#bandPassFilter = @instantiateModel(TestBandPassFilter)
#simulate!(bandPassFilter, Tsit5(), stopTime = 50u"s")
#plot(bandPassFilter, ["u", "y"], figure=2)


# 2.5 Physically oriented modeling

# 2.5.1 Connectors

Pin = Model( v = potential, i = flow )


# 2.5.2 Components

Resistor = Model(
    R = 1.0u"Ω",
    p = Pin,
    n = Pin,
    equations = :[
        0 = p.i + n.i
        v = p.v - n.v
        i = p.i
        R*i = v ]
    )
    
    
# 2.5.3 Inheritance

OnePort = Model(
    p = Pin,
    n = Pin,
    partialEquations = :[
        0 = p.i + n.i
        v = p.v - n.v
        i = p.i ] )
        
Resistor = OnePort | Model( R = 1.0u"Ω", equation = :[ R*i = v ], )

Capacitor = OnePort | Model( C = 1.0u"F", v=Map(init=0.0u"V"), equation = :[ C*der(v) = i ] )

Inductor = OnePort | Model( L = 1.0u"H", i=Map(init=0.0u"A"), equation = :[ L*der(i) = v ] )

ConstantVoltage = OnePort | Model( V = 1.0u"V", equation = :[ v = V ] )        



# 2.5.5 Connected models

Filter = (
    R = Resistor | Map(R=0.5u"Ω"),
    C = Capacitor | Map(C=2.0u"F"),
    V = ConstantVoltage | Map(V=10.0u"V"),
    connect = :[
      (V.p, R.p)
      (R.n, C.p)
      (C.n, V.n)
    ]
)


# 2.5.6 Parameter propagation

Filter2 = Model(
    r = 2.0u"Ω",
    c = 1.0u"F",
    v = 10u"V",
    R = Resistor | Map(R=:r),
    C = Capacitor | Map(C=:c),
    V = ConstantVoltage | Map(V=:v),
    connect = :[
      (V.p, R.p)
      (R.n, C.p)
      (C.n, V.n)
    ]
)

TwoFilters = Model( f1 = Filter | Map( r = 10.0, c = 2.0), f2 = Filter )


# 2.5.7 Redeclarations

VoltageDividerAndFilter = TwoFilters | Map(f1 = Map(C = Redeclare | Resistor | Map(R = 20.0)))


# 2.6 Arrays


StateSpace = Model(
    A = fill(0.0,0,0),
    B = fill(0.0,0,0),
    C = fill(0.0,0,0),
    D = fill(0.0,0,0),
    u = input,
    y = output,
    x = Var(init = zeros(0)),
    equations = :[
        der(x) = A*x + B*u
             y = C*x + D*u
    ]
)

col(args...) = hvcat(1, args...)  # Construct a column matrix from a vector

SecondOrder = Model(
    w = 20.0,
    D =  0.1,
    k =  2.0,
    sys = StateSpace | Map(A = :([  0        1;
                                 -w^2  -2*D*w]),
                           B = :(col([0; w^2])),
                           C = :([k 0]),
                           D = :(zeros(1,1)),
                           x = Var(init = zeros(2)) ),
    equations = :[sys.u = [1.0]]
)


# 2.7 Model libraries

FilterCircuit = Model(
    R = Modia.Resistor  | Map(R=0.5u"Ω"),
    C = Modia.Capacitor | Map(C=2.0u"F", v=Var(init=0.1u"V")),
    V = Modia.ConstantVoltage | Map(V=10.0u"V"),
    ground = Modia.Ground,
    connect = :[
      (V.p, R.p)
      (R.n, C.p)
      (C.n, V.n, ground.p)
    ]
)

filterCircuit = @instantiateModel(FilterCircuit)
simulate!(filterCircuit, Tsit5(), stopTime=10.0)
plot(filterCircuit, ["C.v", "C.i"], figure=3)

end
module TestLinearEquationSystemWithUnitsAndMonteCarlo

using ModiaLang
using ModiaLang.Unitful
using ModiaLang.DifferentialEquations
using ModiaLang.MonteCarloMeasurements
@usingModiaPlot


include("../models/Electric.jl")

### With units

FilterCircuit = Model(
    R1 = Resistor  | Map(R=100u"Ω"),
    R2 = Resistor  | Map(R=200u"Ω", i=Var(start=0.0u"A")),    
    Ri = Resistor  | Map(R=10u"Ω"),
    C  = Capacitor | Map(C=2.5e-3u"F", v=Var(init=0.0u"V"), i=Var(start=0.0u"A")),
    V  = ConstantVoltage | Map(V=10.0u"V"),
    ground = Ground,     
    connect = :[
      (V.p, Ri.n)
      (Ri.p, R1.p)
      (Ri.p, R2.p)
      (R1.n, R2.n, C.p)
      (C.n, V.n, ground.p)
    ]
)


filterCircuit1 = @instantiateModel(FilterCircuit, unitless=false)
simulate!(filterCircuit1, QBDF(autodiff=false), stopTime = 1.0, log=true) 
plot(filterCircuit1, ("V.v", "Ri.v", "R1.v", "C.v"), figure=1)


filterCircuit2 = @instantiateModel(FilterCircuit, unitless=false, FloatType=StaticParticles{Float64,100})
simulate!(filterCircuit2, QBDF(autodiff=false), stopTime = 1.0, log=true, 
          merge = Map(R1=Map(R=(100∓10)u"Ω"),
                      R2=Map(R=(200∓20)u"Ω"),
                      C =Map(C=(2.5e-3∓1e-4)u"F")
                     )
         ) 
plot(filterCircuit2, ("V.v", "Ri.v", "R1.v", "C.v"), figure=2)


filterCircuit3 = @instantiateModel(FilterCircuit, unitless=false, FloatType=Particles{Float64,2000})
simulate!(filterCircuit3, QBDF(autodiff=false), stopTime = 1.0, log=true, 
          merge = Map(R1=Map(R=(100±10)u"Ω"),
                      R2=Map(R=(200±20)u"Ω"),
                      C =Map(C=(2.5e-3±1e-4)u"F")
                     )
         ) 
plot(filterCircuit3, ("V.v", "Ri.v", "R1.v", "C.v"), figure=3, MonteCarloAsArea=true)


### Without units

filterCircuit4 = @instantiateModel(FilterCircuit, unitless=true)
simulate!(filterCircuit4, QBDF(autodiff=false), stopTime = 1.0, log=true) 
plot(filterCircuit4, ("V.v", "Ri.v", "R1.v", "C.v"), figure=4)


filterCircuit5 = @instantiateModel(FilterCircuit, unitless=true, FloatType=StaticParticles{Float64,100})
simulate!(filterCircuit5, QBDF(autodiff=false), stopTime = 1.0, log=true, 
          merge = Map(R1=Map(R=100∓10),
                      R2=Map(R=200∓20),
                      C =Map(C=2.5e-3∓1e-4)
                     )
         ) 
plot(filterCircuit5, ("V.v", "Ri.v", "R1.v", "C.v"), figure=5)


filterCircuit6 = @instantiateModel(FilterCircuit, unitless=true, FloatType=Particles{Float64,2000})
simulate!(filterCircuit6, QBDF(autodiff=false), stopTime = 1.0, log=true, 
          merge = Map(R1=Map(R=100±10),
                      R2=Map(R=200±20),
                      C =Map(C=2.5e-3±1e-4)
                     )
         ) 
plot(filterCircuit6, ("V.v", "Ri.v", "R1.v", "C.v"), figure=6, MonteCarloAsArea=true)


end
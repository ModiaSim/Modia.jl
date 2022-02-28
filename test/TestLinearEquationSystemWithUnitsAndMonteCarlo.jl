module TestLinearEquationSystemWithUnitsAndMonteCarlo

using Modia
using Modia.MonteCarloMeasurements
@usingModiaPlot


include("../models/Electric.jl")

### With units

FilterCircuit = Model(
    R1 = Resistor  | Map(R=100.0u"Ω"),
    R2 = Resistor  | Map(R=200.0u"Ω", i=Var(start=0.0u"A")),    
    Ri = Resistor  | Map(R=10.0u"Ω"),
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
plotVariables = ("V.v", "Ri.v", "R1.v", "C.v")
logCode = false


# Units defined (and unitless=false and true)
filterCircuit1 = @instantiateModel(FilterCircuit, unitless=true, logCode=logCode)
simulate!(filterCircuit1, QBDF(autodiff=false), stopTime = 1.0, log=true) 
plot(filterCircuit1, plotVariables, figure=1)

filterCircuit2 = @instantiateModel(FilterCircuit, unitless=false, logCode=logCode)
simulate!(filterCircuit2, QBDF(autodiff=false), stopTime = 1.0, log=true) 
plot(filterCircuit2, plotVariables, figure=2)


# Units and StaticParticles defined (and unitless=false and true)
FilterCircuitStaticParticles = FilterCircuit | Map(R1=Map(R=(100.0∓10.0)u"Ω"),
                                                   R2=Map(R=(200.0∓20.0)u"Ω"),
                                                   C =Map(C=(2.5e-3∓1e-4)u"F")
                                                  )
filterCircuit3 = @instantiateModel(FilterCircuitStaticParticles, unitless=true, FloatType=StaticParticles{Float64,100}, logCode=logCode)
simulate!(filterCircuit3, QBDF(autodiff=false), stopTime = 1.0, log=true, 
          merge = Map(R1=Map(R=(110.0∓10.0)u"Ω"),
                      R2=Map(R=(220.0∓20.0)u"Ω"),
                      C =Map(C=(2.6e-3∓1e-4)u"F")
                     )
         ) 
plot(filterCircuit3, plotVariables, figure=3)

filterCircuit4 = @instantiateModel(FilterCircuitStaticParticles, unitless=false, FloatType=StaticParticles{Float64,100}, logCode=logCode)
simulate!(filterCircuit4, QBDF(autodiff=false), stopTime = 1.0, log=true, 
          merge = Map(R1=Map(R=(110.0∓10.0)u"Ω"),
                      R2=Map(R=(220.0∓20.0)u"Ω"),
                      C =Map(C=(2.6e-3∓1e-4)u"F")
                     )
         ) 
plot(filterCircuit4, plotVariables, figure=4)



# Units and Particles defined (and unitless=false and true)
FilterCircuitParticles = FilterCircuit | Map(R1=Map(R=(100.0±0.0)u"Ω"),
                                             R2=Map(R=(200.0±20.0)u"Ω"),
                                             C =Map(C=(2.5e-3±1e-4)u"F")
                                            )             
filterCircuit5 = @instantiateModel(FilterCircuitParticles, unitless=true, FloatType=Particles{Float64,2000}, logCode=logCode)
simulate!(filterCircuit5, QBDF(autodiff=false), stopTime = 1.0, log=true, 
          merge = Map(R1=Map(R=(110.0±10.0)u"Ω"),
                      R2=Map(R=(220.0±20.0)u"Ω"),
                      C =Map(C=(2.6e-3±1e-4)u"F")
                     )
         ) 
plot(filterCircuit5, plotVariables, figure=5, MonteCarloAsArea=true)

filterCircuit6 = @instantiateModel(FilterCircuitParticles, unitless=false, FloatType=Particles{Float64,2000}, logCode=logCode)
simulate!(filterCircuit6, QBDF(autodiff=false), stopTime = 1.0, log=true, 
          merge = Map(R1=Map(R=(110.0±10.0)u"Ω"),
                      R2=Map(R=(220.0±20.0)u"Ω"),
                      C =Map(C=(2.6e-3±1e-4)u"F")
                     )
         ) 
plot(filterCircuit6, plotVariables, figure=6, MonteCarloAsArea=true)


end
module BouncingBallSimulation

using TinyModia
using DifferentialEquations
using ModiaPlot

BouncingBall = Model(
    e = 0.7,
    g = 9.81,
    h = Var(init = 1.0),
    v = Var(init = 0.0),
    flying = Var(start = true),   # due to pre(flying) -> SimulationModel(....; pre_startValues = [true], ....)
    equations = :[
        # desired: 
        #   flying = edge(-h) ? reinit(v, -e*v) : pre(flying)
        flying = edge(instantiatedModel, 1, -h, "-h", _leq_mode) ?
                    !reinit(instantiatedModel, _x, 2, -e*v, _leq_mode) : pre(instantiatedModel, 1),
        der(h) = v,
        der(v) = flying ? -g : 0.0
    ]
)


model = @instantiateModel(BouncingBall, logCode=true)

@show model.equationInfo.x_info

model.pre = [true]

simulate!(model, Tsit5(), stopTime = 3.0, nz=1, log=true, logEvents=true)   # requiredFinalStates = [-0.3617373025974107]

plot(model, ["h", "v", "flying"])

end